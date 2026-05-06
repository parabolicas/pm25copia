# ============================================================
# SCRIPT: Validar arquivo CAMS EAC4 PM2.5 (NetCDF)
# Confirma que a região corresponde ao Estado de São Paulo
# ============================================================

# Instalar pacotes necessários (só precisa rodar uma vez)
if (!require("ncdf4"))   install.packages("ncdf4")
if (!require("fields"))  install.packages("fields")
if (!require("maps"))    install.packages("maps")

library(ncdf4)
library(fields)
library(maps)

# ------------------------------------------------------------
# 1. ABRIR O ARQUIVO NetCDF
# ------------------------------------------------------------
# ALTERE O CAMINHO abaixo para onde você salvou o arquivo
arquivo <- "~/Downloads/data_sfc.nc"

nc <- nc_open(arquivo)

cat("============================================\n")
cat("VALIDAÇÃO DO ARQUIVO CAMS EAC4\n")
cat("============================================\n\n")

# ------------------------------------------------------------
# 2. INFORMAÇÕES GERAIS
# ------------------------------------------------------------
cat("--- INFORMAÇÕES GERAIS ---\n")
print(nc)

# ------------------------------------------------------------
# 3. COORDENADAS
# ------------------------------------------------------------
lat <- ncvar_get(nc, "latitude")
lon <- ncvar_get(nc, "longitude")
time <- ncvar_get(nc, "time")

cat("\n--- COORDENADAS ---\n")
cat("Latitudes:", min(lat), "a", max(lat), "(", length(lat), "pontos )\n")
cat("Longitudes:", min(lon), "a", max(lon), "(", length(lon), "pontos )\n")
cat("Resolução lat:", abs(diff(lat[1:2])), "graus\n")
cat("Resolução lon:", abs(diff(lon[1:2])), "graus\n")
cat("Pontos temporais:", length(time), "\n")

# Verificar se cobre SP
cat("\n--- VERIFICAÇÃO REGIÃO SP ---\n")
sp_lat_ok <- min(lat) <= -26 & max(lat) >= -19
sp_lon_ok <- min(lon) <= -54 & max(lon) >= -44
if (sp_lat_ok & sp_lon_ok) {
  cat("✓ REGIÃO CORRETA - Cobre o Estado de São Paulo\n")
} else {
  cat("✗ ATENÇÃO - Região pode não cobrir SP adequadamente\n")
  cat("  Esperado: lat -26 a -19, lon -54 a -44\n")
}

# ------------------------------------------------------------
# 4. DADOS DE PM2.5
# ------------------------------------------------------------
pm25 <- ncvar_get(nc, "pm2p5")
cat("\n--- DADOS PM2.5 ---\n")
cat("Dimensões:", paste(dim(pm25), collapse=" x "), "\n")
cat("Unidade original: kg/m³\n")

# Converter para µg/m³
pm25_ugm3 <- pm25 * 1e9

cat("Faixa de valores (µg/m³):", round(min(pm25_ugm3, na.rm=T), 2),
    "a", round(max(pm25_ugm3, na.rm=T), 2), "\n")
cat("Média geral (µg/m³):", round(mean(pm25_ugm3, na.rm=T), 2), "\n")

# ------------------------------------------------------------
# 5. CALCULAR MÉDIA ANUAL
# ------------------------------------------------------------
# Média ao longo de todos os tempos para cada ponto do grid
pm25_media <- apply(pm25_ugm3, c(1,2), mean, na.rm=TRUE)

cat("\nMédia anual PM2.5 (µg/m³):\n")
cat("  Mínima:", round(min(pm25_media, na.rm=T), 1), "\n")
cat("  Máxima:", round(max(pm25_media, na.rm=T), 1), "\n")
cat("  Média:", round(mean(pm25_media, na.rm=T), 1), "\n")

# ------------------------------------------------------------
# 6. GERAR MAPA DE VALIDAÇÃO
# ------------------------------------------------------------
cat("\nGerando mapa...\n")

# Salvar como PNG
png("CAMS_PM25_validacao_mapa.png", width=1000, height=800, res=120)

# Criar grid de coordenadas
image.plot(lon, lat, t(pm25_media),
           main = "CAMS EAC4 - PM2.5 Média Anual 2008 (µg/m³)",
           xlab = "Longitude", ylab = "Latitude",
           col = rev(heat.colors(20)),
           zlim = c(0, max(pm25_media, na.rm=T)))

# Adicionar contorno do Brasil/SP
map("world", regions = "Brazil", add = TRUE, col = "black", lwd = 2)

# Marcar algumas cidades de referência
points(-46.63, -23.55, pch=19, cex=1.5, col="blue")  # São Paulo
text(-46.63, -23.55, "São Paulo", pos=4, col="blue", font=2, cex=0.8)

points(-47.06, -22.91, pch=19, cex=1.2, col="blue")  # Campinas
text(-47.06, -22.91, "Campinas", pos=4, col="blue", cex=0.7)

points(-49.38, -20.82, pch=19, cex=1.2, col="blue")  # S.J.Rio Preto
text(-49.38, -20.82, "S.J.Rio Preto", pos=4, col="blue", cex=0.7)

points(-47.81, -21.18, pch=19, cex=1.2, col="blue")  # Ribeirão Preto
text(-47.81, -21.18, "Rib. Preto", pos=4, col="blue", cex=0.7)

points(-46.32, -23.96, pch=19, cex=1.2, col="blue")  # Santos
text(-46.32, -23.96, "Santos", pos=4, col="blue", cex=0.7)

points(-51.39, -22.12, pch=19, cex=1.2, col="blue")  # Pres. Prudente
text(-51.39, -22.12, "P.Prudente", pos=4, col="blue", cex=0.7)

# Linha da diretriz OMS (15 µg/m³)
contour(lon, lat, t(pm25_media), levels = 15,
        add = TRUE, col = "red", lwd = 2, labcex = 0.8)

# Legenda
legend("bottomleft",
       legend = c("Cidades SP", "Contorno OMS (15 µg/m³)"),
       pch = c(19, NA), lty = c(NA, 1),
       col = c("blue", "red"), lwd = c(NA, 2),
       bg = "white", cex = 0.8)

dev.off()

cat("Mapa salvo: CAMS_PM25_validacao_mapa.png\n")

# ------------------------------------------------------------
# 7. EXTRAIR PM2.5 PARA CIDADES-CHAVE
# ------------------------------------------------------------
cat("\n--- PM2.5 MÉDIO ANUAL POR CIDADE (2008) ---\n")

cidades <- data.frame(
  nome = c("São Paulo (centro)", "Campinas", "S.J.Rio Preto",
           "Ribeirão Preto", "Santos", "Sorocaba",
           "Piracicaba", "Pres. Prudente", "S.J.Campos"),
  lat = c(-23.55, -22.91, -20.82, -21.18, -23.96, -23.50,
          -22.73, -22.12, -23.18),
  lon = c(-46.63, -47.06, -49.38, -47.81, -46.32, -47.46,
          -47.65, -51.39, -45.88)
)

for (i in 1:nrow(cidades)) {
  idx_lat <- which.min(abs(lat - cidades$lat[i]))
  idx_lon <- which.min(abs(lon - cidades$lon[i]))
  valor <- round(pm25_media[idx_lon, idx_lat], 1)
  cat(sprintf("  %-22s (%.1f, %.1f): %5.1f µg/m³\n",
              cidades$nome[i], cidades$lat[i], cidades$lon[i], valor))
}

nc_close(nc)

cat("\n============================================\n")
cat("VALIDAÇÃO COMPLETA\n")
cat("============================================\n")

# # Inverter latitude para ordem crescente
lat_ord <- rev(lat)
pm25_media_ord <- pm25_media[, rev(seq_along(lat))]

# Mapa
png("CAMS_PM25_validacao_SP.png", width = 1000, height = 800, res = 120)
image.plot(lon, lat_ord, pm25_media_ord,
           main = "CAMS EAC4 - PM2.5 Média Anual 2008 (µg/m³)",
           xlab = "Longitude", ylab = "Latitude",
           col = rev(heat.colors(20)))
map("world", regions = "Brazil", add = TRUE, col = "black", lwd = 2)
points(-46.63, -23.55, pch = 19, cex = 1.5, col = "blue")
text(-46.63, -23.55, "São Paulo", pos = 4, col = "blue", font = 2)
points(-47.06, -22.91, pch = 19, col = "blue")
text(-47.06, -22.91, "Campinas", pos = 4, col = "blue", cex = 0.8)
points(-49.38, -20.82, pch = 19, col = "blue")
text(-49.38, -20.82, "S.J.Rio Preto", pos = 4, col = "blue", cex = 0.8)
dev.off()

cat("Mapa salvo: CAMS_PM25_validacao_SP.png\n")

# Valores por cidade (isso não muda)
cat("\n--- PM2.5 por cidade (média 2008, µg/m³) ---\n")
cidades <- list(
  c("São Paulo", -23.55, -46.63),
  c("Campinas", -22.91, -47.06),
  c("Ribeirão Preto", -21.18, -47.81),
  c("S.J.Rio Preto", -20.82, -49.38),
  c("Santos", -23.96, -46.32),
  c("Sorocaba", -23.50, -47.46),
  c("Piracicaba", -22.73, -47.65)
)

for (c in cidades) {
  idx_lat <- which.min(abs(lat - as.numeric(c[2])))
  idx_lon <- which.min(abs(lon - as.numeric(c[3])))
  val <- round(pm25_media[idx_lon, idx_lat], 1)
  cat(sprintf("  %-18s: %.1f µg/m³\n", c[1], val))
}

nc_close(nc) 