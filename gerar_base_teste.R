# =============================================================================
# gerar_base_teste.R — Base sintética para validar matching_cohort.R
# Projeto: PM2.5 e CPNPC — Doutorado Direto FMUSP
#
# Cria uma base hipotética com 500 pacientes (~165 casos + ~335 controles)
# distribuídos em 5 mesorregiões, com sexo e idade plausíveis.
# Permite testar o pareamento ANTES de ter dados reais da coleta.
#
# Uso:
#   Rscript gerar_base_teste.R
# =============================================================================

set.seed(42)
n <- 500

test_df <- data.frame(
  id_paciente = paste0("P", sprintf("%04d", 1:n)),
  is_case     = rbinom(n, 1, 0.33),                # ~33% casos
  idade       = round(runif(n, 35, 80)),
  sexo        = sample(c("M", "F"), n, replace = TRUE, prob = c(0.5, 0.5)),
  mesorregiao = sample(
    c("RMSP", "Campinas", "Ribeirao Preto", "Bauru", "Aracatuba"),
    n, replace = TRUE
  )
)

write.csv(test_df, "base_clinica.csv", row.names = FALSE)

cat("Base sintetica criada: base_clinica.csv (",
    nrow(test_df), "registros)\n")
cat("  Casos:    ", sum(test_df$is_case == 1), "\n")
cat("  Controles:", sum(test_df$is_case == 0), "\n")
cat("\nProximo passo: Rscript matching_cohort.R\n")
