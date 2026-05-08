# Script para Pareamento 1:2 (Casos vs Controles) usando R e MatchIt
# Requisitos: install.packages(c("MatchIt", "dplyr", "cobalt"))

library(MatchIt)
library(dplyr)
library(cobalt) # Excelente para plots de balanceamento (Love plots)

# 1. Carregar os dados
# Substitua 'base_clinica.csv' pelo caminho real do seu arquivo gerado pelo pipeline Python
# A base deve conter as colunas: id_paciente, is_case (1 para caso, 0 para controle),
# idade, sexo, mesorregiao.
df <- read.csv("base_clinica.csv")

# Tratamento de NAs (MatchIt não aceita NAs nas variáveis de pareamento)
df <- df %>% filter(!is.na(idade), !is.na(sexo), !is.na(mesorregiao))

# 2. Executar o Pareamento 1:2
# Método: Nearest Neighbor (pode ser trocado por 'optimal' se o dataset não for gigante)
# Exact: pareamento estrito por sexo e mesorregião
# Caliper na idade: garante que a diferença de idade seja no máximo 5 anos (±5)
# std.caliper = FALSE garante que o caliper seja em unidades brutas (anos) e não em desvios-padrão

cat("Iniciando o pareamento 1:2...\n")

m.out <- matchit(
  formula = is_case ~ idade,
  data = df,
  method = "nearest",       # Use "optimal" para pareamento ótimo global (requer pacote 'optmatch')
  distance = "mahalanobis", # A distância será baseada puramente na idade
  exact = ~ sexo + mesorregiao,
  ratio = 2,                # Pareamento 1 caso para 2 controles
  caliper = c(idade = 5),   # Tolerância máxima de ±5 anos
  std.caliper = FALSE,      # Usa anos brutos, não desvios-padrão
  replace = FALSE           # Sem reposição
)

# 3. Diagnósticos de Balanceamento
cat("\n--- Sumário do Pareamento ---\n")
print(summary(m.out, un = FALSE))

# Salvar gráficos de diagnóstico em PDF
cat("\nGerando gráficos de diagnóstico (Love Plot)...\n")
pdf("diagnostico_pareamento.pdf", width = 8, height = 6)

# Love plot usando o pacote cobalt (padrão ouro em epidemiologia)
love.plot(m.out,
          binary = "std",
          thresholds = c(m = .1), # Threshold aceitável de Standardized Mean Difference (SMD)
          title = "Balanceamento de Covariáveis: Casos vs Controles (1:2)")

dev.off()

# 4. Extrair o Dataset Pareado
matched_data <- match.data(m.out)

# O MatchIt adiciona três colunas úteis:
# distance: a distância usada no pareamento
# weights: peso da observação (1 para todos no matching 1:1 ou 1:k sem reposição)
# subclass: ID do par/estrato (útil para regressão logística condicional!)

# 5. Salvar a base pareada para continuar o pipeline em Python
write.csv(matched_data, "base_clinica_pareada_1_2.csv", row.names = FALSE)

cat("\nPareamento concluído com sucesso. Base salva como 'base_clinica_pareada_1_2.csv'.\n")
cat("Nota para a análise em Python: use a coluna 'subclass' como o ID do grupo na Regressão Logística Condicional.\n")
