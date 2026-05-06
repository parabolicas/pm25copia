# ============================================================
# SCRIPT: Download PM2.5 de TODAS as estações CETESB (2012-2025)
# Projeto: PM2.5 e CPNPC - Doutorado FMUSP
# Método: Web scraping direto do Qualar V3.83
# ============================================================

library(httr)
library(qualR)

# ------------------------------------------------------------
# 1. CREDENCIAIS E CONFIGURAÇÃO
# ------------------------------------------------------------
meu_usuario <- "vitor.felixx@gmail.com"
minha_senha  <- "jmc7qwB7NY9TZ9."
param_pm25   <- 57
anos         <- 2012:2025
estacoes     <- cetesb_aqs[, c("name", "code")]

dir.create("dados_pm25_cetesb", showWarnings = FALSE)
write.csv(cetesb_aqs, "dados_pm25_cetesb/cetesb_estacoes_coordenadas.csv",
          row.names = FALSE)

cat("============================================\n")
cat("Download de PM2.5 - CETESB Qualar\n")
cat("Estações:", nrow(estacoes), "| Anos:", min(anos), "a", max(anos), "\n")
cat("============================================\n\n")

# ------------------------------------------------------------
# 2. LOGIN
# ------------------------------------------------------------
cat("Fazendo login...\n")
res <- GET("https://qualar.cetesb.sp.gov.br/qualar/home.do")
my_cookie <- setNames(cookies(res)$value, cookies(res)$name)

login_res <- POST(
  "https://qualar.cetesb.sp.gov.br/qualar/autenticador",
  body = list(
    cetesb_login    = meu_usuario,
    cetesb_password = minha_senha,
    enviar          = "Enviar"
  ),
  encode = "form",
  set_cookies(my_cookie)
)
all_cookies <- c(my_cookie, setNames(cookies(login_res)$value, cookies(login_res)$name))
cat("Login status:", status_code(login_res), "\n\n")

# ------------------------------------------------------------
# 3. DOWNLOAD PRINCIPAL
# ------------------------------------------------------------
todos_dados <- list()
contador <- 1
estacoes_com_dados <- c()

for (ano in anos) {
  
  data_inicio <- paste0("01/01/", ano)
  data_fim    <- paste0("31/12/", ano)
  
  cat("\n======== ANO:", ano, "========\n")
  
  for (i in 1:nrow(estacoes)) {
    
    nome <- estacoes$name[i]
    cod  <- estacoes$code[i]
    cat("  ", nome, "(", cod, ") ... ")
    
    # --- Baixar CSV ---
    csv_res <- tryCatch(
      POST(
        "https://qualar.cetesb.sp.gov.br/qualar/exportaDadosAvanc.do",
        query = list(method = "exportar"),
        body = list(
          dataInicialStr         = data_inicio,
          dataFinalStr           = data_fim,
          estacaoVO.nestcaMonto  = as.character(cod),
          nparmtsSelecionados    = as.character(param_pm25)
        ),
        encode = "form",
        set_cookies(all_cookies)
      ),
      error = function(e) { cat("ERRO CONEXAO\n"); NULL }
    )
    
    if (is.null(csv_res) || status_code(csv_res) != 200) {
      cat("FALHA (status:", ifelse(is.null(csv_res), "NULL", status_code(csv_res)), ")\n")
      next
    }
    
    # Verificar se retornou CSV
    ct <- headers(csv_res)$`content-type`
    if (is.null(ct) || !grepl("csv", ct)) {
      cat("NAO CSV (content-type:", ct, ")\n")
      next
    }
    
    # --- Ler e processar ---
    csv_raw <- content(csv_res, "raw")
    csv_texto <- tryCatch(
      iconv(rawToChar(csv_raw), from = "CP1252", to = "UTF-8"),
      error = function(e) NULL
    )
    
    if (is.null(csv_texto)) { cat("ERRO ENCODING\n"); next }
    
    linhas <- strsplit(csv_texto, "\n")[[1]]
    idx_header <- grep("^Data;Hora", linhas)
    
    if (length(idx_header) == 0) { cat("SEM HEADER\n"); next }
    
    idx_dados <- (idx_header + 2):length(linhas)
    dados_linhas <- linhas[idx_dados]
    dados_linhas <- dados_linhas[nchar(trimws(dados_linhas)) > 0]
    
    if (length(dados_linhas) == 0) { cat("SEM DADOS\n"); next }
    
    df <- tryCatch({
      dados_texto <- paste(dados_linhas, collapse = "\n")
      read.csv(text = dados_texto, header = FALSE, sep = ";",
               stringsAsFactors = FALSE, na.strings = c("", " "))
    }, error = function(e) NULL)
    
    if (is.null(df) || nrow(df) == 0) { cat("PARSE ERRO\n"); next }
    
    names(df) <- c("data", "hora", "pm25")
    df$pm25 <- as.numeric(gsub(",", ".", df$pm25))
    
    n_validos <- sum(!is.na(df$pm25))
    if (n_validos == 0) { cat("APENAS NA\n"); next }
    
    # Adicionar metadados
    df$estacao_codigo <- cod
    df$estacao_nome   <- nome
    df$ano            <- ano
    idx_aqs <- which(cetesb_aqs$code == cod)
    df$lat <- ifelse(length(idx_aqs) > 0, cetesb_aqs$lat[idx_aqs], NA)
    df$lon <- ifelse(length(idx_aqs) > 0, cetesb_aqs$lon[idx_aqs], NA)
    
    media <- round(mean(df$pm25, na.rm = TRUE), 1)
    cat("OK -", n_validos, "medições, média:", media, "ug/m3\n")
    
    todos_dados[[contador]] <- df
    contador <- contador + 1
    estacoes_com_dados <- unique(c(estacoes_com_dados, nome))
    
    Sys.sleep(1)
  }
  
  # Re-login a cada ano
  cat("\n>> Re-autenticando... ")
  res <- GET("https://qualar.cetesb.sp.gov.br/qualar/home.do")
  my_cookie <- setNames(cookies(res)$value, cookies(res)$name)
  login_res <- POST(
    "https://qualar.cetesb.sp.gov.br/qualar/autenticador",
    body = list(
      cetesb_login    = meu_usuario,
      cetesb_password = minha_senha,
      enviar          = "Enviar"
    ),
    encode = "form",
    set_cookies(my_cookie)
  )
  all_cookies <- c(my_cookie, setNames(cookies(login_res)$value, cookies(login_res)$name))
  cat("OK\n")
  
  # Salvar parcial
  if (length(todos_dados) > 0) {
    parcial <- do.call(rbind, todos_dados)
    write.csv(parcial,
              paste0("dados_pm25_cetesb/pm25_parcial_ate_", ano, ".csv"),
              row.names = FALSE)
    cat(">> Parcial salvo:", nrow(parcial), "registros\n")
  }
}

# ------------------------------------------------------------
# 4. CONSOLIDAR
# ------------------------------------------------------------
if (length(todos_dados) > 0) {
  
  df_completo <- do.call(rbind, todos_dados)
  
  write.csv(df_completo,
            "dados_pm25_cetesb/PM25_TODAS_ESTACOES_2012_2025.csv",
            row.names = FALSE)
  
  resumo <- aggregate(
    df_completo$pm25,
    by = list(
      estacao_nome   = df_completo$estacao_nome,
      estacao_codigo = df_completo$estacao_codigo,
      ano            = df_completo$ano,
      lat            = df_completo$lat,
      lon            = df_completo$lon
    ),
    FUN = function(x) round(mean(x, na.rm = TRUE), 1)
  )
  names(resumo)[6] <- "media_anual_pm25_ugm3"
  resumo <- resumo[order(resumo$estacao_nome, resumo$ano), ]
  
  write.csv(resumo,
            "dados_pm25_cetesb/RESUMO_MEDIAS_ANUAIS_PM25.csv",
            row.names = FALSE)
  
  cat("\n============================================\n")
  cat("DOWNLOAD COMPLETO!\n")
  cat("Total registros:", nrow(df_completo), "\n")
  cat("Estações com dados:", length(unique(df_completo$estacao_nome)), "\n")
  cat("Período:", min(df_completo$ano), "-", max(df_completo$ano), "\n")
  cat("Arquivos em: dados_pm25_cetesb/\n")
  cat("============================================\n")
} else {
  cat("\nNenhum dado baixado.\n")
}
