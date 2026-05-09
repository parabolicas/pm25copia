# =============================================================================
# matching_cohort.R — Pareamento 1:2 (Casos vs Controles)
# Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP
#
# Objetivo:
#   Parear cada caso de CPNPC (nunca/leves tabagistas) com 2 controles
#   oncológicos não-pulmonares, conforme Seção 4.1 do protocolo:
#     - Pareamento exato por: sexo + mesorregiao
#     - Pareamento por proximidade: idade (caliper ±5 anos)
#     - Razão 1:2 (um caso para dois controles, sem reposição)
#
# Entradas (CSV):
#   base_clinica.csv com colunas:
#     id_paciente   : identificador único (string ou int)
#     is_case       : 1 (caso CPNPC) | 0 (controle não-pulmonar)
#     idade         : numérico (anos completos no diagnóstico)
#     sexo          : "M" | "F" (será padronizado)
#     mesorregiao   : nome da mesorregiao SP (será padronizado)
#
# Saídas (em ./matching_output/):
#   base_clinica_pareada_1_2.csv     — base pareada (com subclass = id do par)
#   casos_nao_pareados.csv           — casos descartados pelo caliper
#   balanceamento_pre.csv            — SMD/médias antes do matching
#   balanceamento_pos.csv            — SMD/médias depois do matching
#   sumario_consort.txt              — fluxograma CONSORT-like
#   diagnostico_pareamento.pdf       — Love plot (cobalt) + propensity overlap
#   sensibilidade_nearest.csv        — pareamento alternativo (Sensibilidade)
#
# Referências metodológicas:
#   - Stuart EA. Matching methods for causal inference: a review. Stat Sci 2010.
#   - Austin PC. Balance diagnostics for comparing the distribution of baseline
#     covariates between treatment groups in propensity-score matched samples.
#     Stat Med 2009;28:3083-3107. (origem do threshold SMD ≤ 0.10)
#   - Greifer N. cobalt: Covariate Balance Tables and Plots. CRAN.
#
# Requisitos de pacotes:
#   install.packages(c("MatchIt", "dplyr", "cobalt", "optmatch"))
# =============================================================================

# ---------------------------------------------------------------------
# 1. SETUP
# ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(MatchIt)
  library(dplyr)
  library(cobalt)
})

# Reprodutibilidade — qualquer ordem aleatória de empate fica determinística
set.seed(20260507)

# Pasta de saída
out_dir <- "matching_output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("=============================================\n")
cat("Pareamento 1:2 — PM2.5 e CPNPC FMUSP\n")
cat("Seed: 20260507  |  Saída: ", out_dir, "\n")
cat("=============================================\n\n")


# ---------------------------------------------------------------------
# 2. CARREGAMENTO E VALIDAÇÃO DA BASE
# ---------------------------------------------------------------------
input_csv <- "base_clinica.csv"
if (!file.exists(input_csv)) {
  stop("Arquivo de entrada nao encontrado: ", input_csv,
       "\n  Gere a partir do pipeline Python (run_pipeline.py).")
}

df <- read.csv(input_csv, stringsAsFactors = FALSE)
cat("Carregado:", nrow(df), "registros de", input_csv, "\n")

# --- Asserções de integridade ---
required_cols <- c("id_paciente", "is_case", "idade", "sexo", "mesorregiao")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Colunas obrigatorias ausentes: ",
       paste(missing_cols, collapse = ", "))
}

# is_case deve ser estritamente 0 ou 1
df <- df %>% filter(!is.na(is_case))
if (!all(df$is_case %in% c(0, 1))) {
  bad_vals <- unique(df$is_case[!df$is_case %in% c(0, 1)])
  stop("is_case contem valores invalidos: ",
       paste(bad_vals, collapse = ", "), ". Deve ser 0 ou 1.")
}

# Idade plausivel (pacientes oncologicos)
df <- df %>% filter(!is.na(idade))
df$idade <- as.numeric(df$idade)
if (any(df$idade < 18 | df$idade > 110, na.rm = TRUE)) {
  bad_n <- sum(df$idade < 18 | df$idade > 110, na.rm = TRUE)
  warning("Idades fora do intervalo [18, 110]: ", bad_n,
          " registros. Removendo.")
  df <- df %>% filter(idade >= 18, idade <= 110)
}

# Padronizar sexo e mesorregiao (case + trimming)
df$sexo <- as.factor(toupper(trimws(df$sexo)))
df$mesorregiao <- as.factor(toupper(trimws(df$mesorregiao)))

# Remover NAs nas variaveis-chave
df <- df %>% filter(!is.na(sexo), !is.na(mesorregiao))

n_pre <- nrow(df)
n_cases_pre <- sum(df$is_case == 1)
n_ctrl_pre <- sum(df$is_case == 0)
cat("Apos validacao:", n_pre, "registros (", n_cases_pre, "casos +",
    n_ctrl_pre, "controles disponiveis)\n")

# Garantir que ha controles suficientes para ratio 1:2
if (n_ctrl_pre < 2 * n_cases_pre) {
  warning("Controles disponiveis (", n_ctrl_pre,
          ") < 2x casos (", 2 * n_cases_pre,
          "). Pareamento ratio=2 pode descartar casos.")
}


# ---------------------------------------------------------------------
# 3. PAREAMENTO PRINCIPAL — OPTIMAL + PROPENSITY SCORE LOGISTICO
# ---------------------------------------------------------------------
# Justificativas das escolhas:
#   - method = "optimal": minimiza distancia total da rede de pares
#     (mais robusto que "nearest" greedy, especialmente quando casos
#     concorrem por controles raros). Requer pacote optmatch.
#   - distance = "glm": propensity score via logit. Mais defensavel que
#     mahalanobis quando ha apenas 1 covariavel continua na formula.
#     Permite incluir mais covariaveis no futuro sem mudar arquitetura.
#   - exact = sexo + mesorregiao: pareamento estrito, conforme protocolo.
#   - ratio = 2, replace = FALSE: 1:2 sem reposicao.
#   - caliper = c(idade = 5), std.caliper = FALSE: ±5 anos brutos.
# ---------------------------------------------------------------------

cat("\n--- PAREAMENTO PRINCIPAL (optimal + glm) ---\n")

# Tentar optimal; se optmatch nao instalado, cair para nearest com aviso
m.out <- tryCatch({
  if (!requireNamespace("optmatch", quietly = TRUE)) {
    stop("Pacote optmatch nao instalado.")
  }
  matchit(
    formula     = is_case ~ idade,
    data        = df,
    method      = "optimal",
    distance    = "glm",
    exact       = ~ sexo + mesorregiao,
    ratio       = 2,
    caliper     = c(idade = 5),
    std.caliper = FALSE,
    replace     = FALSE
  )
}, error = function(e) {
  cat("  AVISO:", conditionMessage(e), "\n")
  cat("  Caindo para method='nearest'. Para usar 'optimal':\n")
  cat("    install.packages('optmatch')\n\n")
  matchit(
    formula     = is_case ~ idade,
    data        = df,
    method      = "nearest",
    distance    = "glm",
    exact       = ~ sexo + mesorregiao,
    ratio       = 2,
    caliper     = c(idade = 5),
    std.caliper = FALSE,
    replace     = FALSE
  )
})

method_used <- m.out$call$method
cat("  Metodo aplicado:", as.character(method_used), "\n")


# ---------------------------------------------------------------------
# 4. DIAGNOSTICOS NUMERICOS DE BALANCEAMENTO
# ---------------------------------------------------------------------
cat("\n--- DIAGNOSTICOS NUMERICOS ---\n")

bal_summary <- summary(m.out, un = TRUE, standardize = TRUE)
print(bal_summary)

# Salvar tabelas
if (!is.null(bal_summary$sum.all)) {
  write.csv(bal_summary$sum.all,
            file.path(out_dir, "balanceamento_pre.csv"),
            row.names = TRUE)
}
if (!is.null(bal_summary$sum.matched)) {
  write.csv(bal_summary$sum.matched,
            file.path(out_dir, "balanceamento_pos.csv"),
            row.names = TRUE)
}

# --- Quantificar descartes (CONSORT-like) ---
matched_data <- match.data(m.out)
ids_pareados <- matched_data$id_paciente

casos_pareados <- matched_data %>% filter(is_case == 1)
controles_pareados <- matched_data %>% filter(is_case == 0)

casos_originais <- df %>% filter(is_case == 1)
casos_nao_pareados <- casos_originais %>%
  filter(!id_paciente %in% ids_pareados)

n_cases_post <- nrow(casos_pareados)
n_cases_lost <- nrow(casos_nao_pareados)
n_ctrl_post <- nrow(controles_pareados)

cat("\n=== Fluxograma CONSORT-like ===\n")
cat("Casos elegiveis pre-matching:    ", n_cases_pre, "\n")
cat("Casos pareados (com 2 controles):", n_cases_post, "\n")
cat("Casos descartados (sem par valido):", n_cases_lost, "\n")
cat("Controles utilizados:            ", n_ctrl_post, "\n")
cat("Controles disponiveis nao usados:", n_ctrl_pre - n_ctrl_post, "\n")
cat("Razao alcancada:", round(n_ctrl_post / max(n_cases_post, 1), 2), ":1\n\n")

# --- Detectar pares incompletos (algum subclass com #ctrls != 2) ---
pair_table <- table(matched_data$subclass, matched_data$is_case)
incomplete_pairs <- which(
  pair_table[, "1"] != 1 | pair_table[, "0"] != 2
)
if (length(incomplete_pairs) > 0) {
  cat("AVISO: pares incompletos detectados:", length(incomplete_pairs),
      "(idealmente 0)\n")
  print(pair_table[incomplete_pairs, , drop = FALSE])
} else {
  cat("Todos os pares estao completos (1 caso + 2 controles cada).\n")
}

# --- Salvar casos descartados para auditoria ---
if (n_cases_lost > 0) {
  write.csv(casos_nao_pareados,
            file.path(out_dir, "casos_nao_pareados.csv"),
            row.names = FALSE)
  cat("\nCasos descartados salvos em:",
      file.path(out_dir, "casos_nao_pareados.csv"), "\n")
}


# ---------------------------------------------------------------------
# 5. ANALISE DE SENSIBILIDADE — NEAREST vs OPTIMAL
# ---------------------------------------------------------------------
# Reportar diferencas: e comum que optimal pareia 1-3% mais casos.
# ---------------------------------------------------------------------
cat("\n--- SENSIBILIDADE: nearest (greedy) ---\n")

m.nearest <- matchit(
  formula     = is_case ~ idade,
  data        = df,
  method      = "nearest",
  distance    = "glm",
  exact       = ~ sexo + mesorregiao,
  ratio       = 2,
  caliper     = c(idade = 5),
  std.caliper = FALSE,
  replace     = FALSE
)
md_nearest <- match.data(m.nearest)
n_cases_nearest <- sum(md_nearest$is_case == 1)
cat("  Nearest pareia:", n_cases_nearest, "casos\n")
cat("  Optimal pareia:", n_cases_post, "casos\n")
cat("  Diferenca:     ", n_cases_post - n_cases_nearest, "casos\n")

write.csv(md_nearest,
          file.path(out_dir, "sensibilidade_nearest.csv"),
          row.names = FALSE)


# ---------------------------------------------------------------------
# 6. DIAGNOSTICOS GRAFICOS (PDF unico, multi-paginas)
# ---------------------------------------------------------------------
cat("\n--- DIAGNOSTICOS GRAFICOS ---\n")

pdf(file.path(out_dir, "diagnostico_pareamento.pdf"),
    width = 9, height = 6.5)

# Pagina 1: Love plot (Standardized Mean Differences)
print(love.plot(
  m.out,
  binary       = "std",
  abs          = TRUE,
  thresholds   = c(m = 0.10),
  var.order    = "unadjusted",
  title        = paste0(
    "Balanceamento de Covariaveis: Casos vs Controles (1:2)\n",
    "Threshold SMD = 0.10 (Austin 2009)  |  Metodo: ",
    as.character(method_used)
  ),
  sample.names = c("Pre-matching", "Pos-matching")
))

# Pagina 2: Propensity score overlap (densidade)
print(bal.plot(
  m.out,
  var.name = "distance",
  which    = "both",
  type     = "density",
  mirror   = TRUE,
  title    = "Distribuicao do Propensity Score: Pre vs Pos-matching"
))

# Pagina 3: Distribuicao de idade pos-matching
print(bal.plot(
  m.out,
  var.name = "idade",
  which    = "both",
  type     = "density",
  title    = "Distribuicao de Idade: Pre vs Pos-matching (caliper +/-5 anos)"
))

invisible(dev.off())
cat("  PDF salvo em:", file.path(out_dir, "diagnostico_pareamento.pdf"),
    "\n")


# ---------------------------------------------------------------------
# 7. SALVAR BASE PAREADA (PRINCIPAL) PARA O PIPELINE PYTHON
# ---------------------------------------------------------------------
output_main <- file.path(out_dir, "base_clinica_pareada_1_2.csv")
write.csv(matched_data, output_main, row.names = FALSE)
cat("\nBase pareada principal:", output_main, "\n")
cat("Colunas adicionadas pelo MatchIt:\n")
cat("  - distance: propensity score logistico\n")
cat("  - weights:  peso da observacao (=1 para 1:k sem reposicao)\n")
cat("  - subclass: ID do par/estrato (USAR como cluster na regressao\n")
cat("              logistica condicional do pipeline Python).\n")


# ---------------------------------------------------------------------
# 8. SUMARIO TEXTO PARA O ARTIGO 1
# ---------------------------------------------------------------------
sumario_path <- file.path(out_dir, "sumario_consort.txt")
sink(sumario_path)
cat("=== Pareamento 1:2 - Sumario CONSORT-like ===\n\n")
cat("Metodologia:\n")
cat("  Pacote: MatchIt (Stuart 2011, R) + cobalt (Greifer)\n")
cat("  Metodo:", as.character(method_used), "\n")
cat("  Distance: propensity score via GLM logistico\n")
cat("  Pareamento exato: sexo + mesorregiao\n")
cat("  Caliper: idade +/- 5 anos (escala bruta)\n")
cat("  Razao: 1:2, sem reposicao\n")
cat("  Seed RNG: 20260507\n\n")

cat("Fluxograma:\n")
cat("  Casos elegiveis:                    ", n_cases_pre, "\n")
cat("  Controles disponiveis:              ", n_ctrl_pre, "\n")
cat("  Casos pareados:                     ", n_cases_post, "\n")
cat("  Casos descartados (caliper/exact):  ", n_cases_lost,
    " (", round(100 * n_cases_lost / n_cases_pre, 1), "%)\n")
cat("  Controles utilizados:               ", n_ctrl_post, "\n")
cat("  Razao alcancada:                    ",
    round(n_ctrl_post / max(n_cases_post, 1), 2), ":1\n\n")

cat("Sensibilidade (nearest vs optimal):\n")
cat("  Nearest pareia:", n_cases_nearest, "casos\n")
cat("  Optimal pareia:", n_cases_post, "casos\n")
cat("  Vantagem do optimal: +", n_cases_post - n_cases_nearest,
    " casos preservados\n\n", sep = "")

cat("Balanceamento (ver balanceamento_pos.csv):\n")
cat("  Threshold de aceitacao SMD <= 0.10 (Austin 2009)\n")

sink()
cat("\nSumario CONSORT salvo em:", sumario_path, "\n")

cat("\n=============================================\n")
cat("PAREAMENTO CONCLUIDO\n")
cat("=============================================\n")
