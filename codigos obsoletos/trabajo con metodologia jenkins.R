# ==============================================================================
# PROYECTO: APTITUD AGROCLIMÁTICA - REGIÓN PUNO
# OBJETIVO: Dos Planificaciones → (1) Espacial General | (2) Producto: Quinua
# Metodología: Jenks Natural Breaks | NASA POWER | SAL250 | ENA 2025
# VERSIÓN FINAL - Con todos los fixes aplicados
# ==============================================================================

library(terra)
library(sf)
library(data.table)
library(classInt)

setwd("E:/estadistica espacial")

# ==============================================================================
# BLOQUE 1: MARCO GEOGRÁFICO
# ==============================================================================
cat("--- Cargando límites de Puno ---\n")

shp_peru_vect  <- vect("gadm41_PER_3.shp")
crs(shp_peru_vect) <- "EPSG:4326"
shp_puno_wgs84 <- shp_peru_vect[shp_peru_vect$NAME_1 == "Puno", ]

crs_wgs84  <- "+proj=longlat +datum=WGS84 +no_defs"
crs_utm19s <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs"
crs(shp_puno_wgs84) <- crs_wgs84

shp_puno_utm <- project(shp_puno_wgs84, crs_utm19s)
writeVector(shp_puno_utm, "distritos_puno_utm.shp", overwrite = TRUE)
cat("Distritos de Puno:", nrow(shp_puno_utm), "\n")

# ==============================================================================
# BLOQUE 2: DEM - UNIFICAR Y REPROYECTAR
# FIX: dem_puno_final ahora se define correctamente
# ==============================================================================
cat("--- Procesando DEM ---\n")

if (file.exists("dem_puno_unificado.tif")) {
  dem_puno_wgs84 <- rast("dem_puno_unificado.tif")
  cat("DEM cargado desde disco.\n")
} else {
  zips_dem     <- list.files("dem_zips/", pattern = "\\.zip$", full.names = TRUE)
  carpeta_temp <- "dem_descomprimido"
  dir.create(carpeta_temp, showWarnings = FALSE)
  for (zip in zips_dem) unzip(zip, exdir = carpeta_temp)
  
  archivos_dem <- list.files(carpeta_temp, pattern = "\\.(hgt|tif)$", full.names = TRUE)
  dem_completo <- do.call(mosaic, lapply(archivos_dem, rast))
  crs(dem_completo) <- "EPSG:4326"
  
  dem_puno_wgs84 <- crop(dem_completo, shp_puno_wgs84)
  dem_puno_wgs84 <- mask(dem_puno_wgs84, shp_puno_wgs84)
  writeRaster(dem_puno_wgs84, "dem_puno_unificado.tif",
              gdal = c("COMPRESS=LZW"), overwrite = TRUE)
}

# FIX CRÍTICO: Reproyectar a UTM 19S → esto es dem_puno_final
dem_puno_final <- project(dem_puno_wgs84, crs_utm19s)
names(dem_puno_final) <- "Elevacion_m"

# ==============================================================================
# BLOQUE 3: NASA POWER → TEMPERATURA Y PRECIPITACIÓN ESPACIAL
# FIX: Nombres de columna reales del CSV (temp_C, precip_mmday, humidity_pct)
# ==============================================================================
cat("--- Integrando NASA POWER ---\n")

nasa_raw <- fread("nasa_power_tidy.csv")

# Promedios con los nombres reales de tu CSV
temp_base_nasa    <- mean(nasa_raw$temp_C,       na.rm = TRUE)
precip_base_nasa  <- mean(nasa_raw$precip_mmday, na.rm = TRUE) * 365  # mm/día → mm/año
altitud_base_nasa <- 4105.46

cat(sprintf("Temperatura base: %.2f °C | Precipitación media: %.1f mm/año\n",
            temp_base_nasa, precip_base_nasa))

# --- Raster de Temperatura (gradiente adiabático) ---
# FIX: Escribir directamente sin aritmética intermedia en RAM
v_dem  <- values(dem_puno_final, mat = FALSE)
v_temp <- temp_base_nasa - 0.0065 * (v_dem - altitud_base_nasa)
v_temp <- pmax(pmin(v_temp, 15), -2)  # Clamp [-2, 15] °C en RAM

raster_temperatura        <- dem_puno_final
values(raster_temperatura) <- v_temp
names(raster_temperatura)  <- "Temperatura_Media_C"
writeRaster(raster_temperatura, "raster_temperatura_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)

# --- Raster de Precipitación (modelo orográfico) ---
delta_alt <- v_dem - altitud_base_nasa
v_prec    <- precip_base_nasa + 0.05 * delta_alt - 0.000005 * delta_alt^2
v_prec    <- pmax(pmin(v_prec, 900), 200)  # Clamp [200, 900] mm en RAM

raster_precip        <- dem_puno_final
values(raster_precip) <- v_prec
names(raster_precip)  <- "Precipitacion_mm"
writeRaster(raster_precip, "raster_precipitacion_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)

cat("Rasters climáticos guardados.\n")
rm(v_dem, v_temp, v_prec, delta_alt); gc()

# ==============================================================================
# BLOQUE 4: SUELOS SAL250
# FIX: id_suelo se recrea en RAM (no depende del .shp anterior)
# ==============================================================================
cat("--- Procesando suelos SAL250 ---\n")

if (file.exists("suelos_puno_utm.shp")) {
  suelos_puno_utm <- vect("suelos_puno_utm.shp")
} else {
  suelos_peru     <- vect("Capacidad_de_Uso_Mayor_ONERN_CUM_geogpsperu_SuyoPomalia_931381206.shp")
  suelos_puno_utm <- project(mask(crop(suelos_peru, shp_puno_wgs84), shp_puno_wgs84), crs_utm19s)
  writeVector(suelos_puno_utm, "suelos_puno_utm.shp", overwrite = TRUE)
}

# FIX: Recrear id_suelo en el objeto (puede no estar en el .shp guardado)
columna_suelo           <- names(suelos_puno_utm)[1]
suelos_puno_utm$id_suelo <- as.numeric(as.factor(
  values(suelos_puno_utm)[, columna_suelo]))

cat("Columna de suelo usada:", columna_suelo,
    "| Clases:", length(unique(suelos_puno_utm$id_suelo)), "\n")

# ==============================================================================
# BLOQUE 5: ENA 2025 - LIMPIEZA Y QUINUA CON JENKS
# ==============================================================================
cat("--- Procesando ENA 2025 ---\n")

normalizar_texto <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  gsub("[^A-Z0-9]", "", x)
}

data_agro  <- fread("Data_Agro_puno.csv")
data_limpia <- data_agro[!is.na(P204_COD) & (!is.na(P210_SUP_1) | !is.na(P219_CANT_1))]
data_limpia <- data_limpia[CCDD == 21]

data_limpia[, Afectado_Helada := ifelse(P223B_3 == 1, 1L, 0L)]
data_limpia[, Afectado_Sequia := ifelse(P223B_6 == 1, 1L, 0L)]
data_limpia[, P204_NOM        := normalizar_texto(P204_NOM)]

distritos_rendimiento <- data_limpia[, .(
  Total_Productores    = .N,
  Sup_Cosechada_Ha     = sum(P217_SUP_1,  na.rm = TRUE),
  Produccion_Total_T   = sum(P219_CANT_1, na.rm = TRUE),
  Rendimiento_Promedio = mean(P219_CANT_1 / (P217_SUP_1 + 0.001), na.rm = TRUE),
  Porcentaje_Heladas   = mean(Afectado_Helada, na.rm = TRUE) * 100,
  Porcentaje_Sequia    = mean(Afectado_Sequia, na.rm = TRUE) * 100
), by = .(CCDI, P204_NOM)]

# Enlace ENA ↔ mapa distrital
shp_puno_utm <- vect("distritos_puno_utm.shp")
shp_puno_utm$NAME_CLEAN <- normalizar_texto(shp_puno_utm$NAME_3)

llaves_mapa <- data.table(
  NAME_CLEAN = unique(shp_puno_utm$NAME_CLEAN),
  CCDI       = seq_along(unique(shp_puno_utm$NAME_CLEAN))
)
tabla_ENA <- merge(distritos_rendimiento, llaves_mapa, by = "CCDI", all.x = TRUE)

# --- Mapa de Quinua ---
tabla_quinua <- tabla_ENA[grepl("QUINUA", P204_NOM)]
mapa_quinua  <- merge(shp_puno_utm, tabla_quinua, by = "NAME_CLEAN", all.x = FALSE)

if (nrow(mapa_quinua) > 0) {
  writeVector(mapa_quinua, "mapa_distrital_quinua_puno.shp", overwrite = TRUE)
  cat("Distritos con datos de quinua:", nrow(mapa_quinua), "\n")
} else {
  cat("Sin coincidencias para QUINUA. Verifica P204_NOM en la consola:\n")
  print(unique(tabla_ENA$P204_NOM)[1:20])
  stop("Revisar el cruce de nombres antes de continuar.")
}

# --- JENKS para Quinua ---
vals_rend   <- mapa_quinua$Rendimiento_Promedio
vals_rend   <- vals_rend[!is.na(vals_rend) & is.finite(vals_rend)]
n_cl        <- min(5, length(unique(vals_rend)) - 1)
jenks_obj   <- classIntervals(vals_rend, n = n_cl, style = "jenks")
cortes_jenks <- unique(jenks_obj$brks)

cat("\nCortes Jenks - Rendimiento QUINUA (T/Ha):\n")
etiq <- c("Muy Bajo","Bajo","Medio","Alto","Muy Alto")
for (i in seq_len(length(cortes_jenks) - 1))
  cat(sprintf("  Clase %d (%s): %.2f - %.2f T/Ha\n",
              i, etiq[i], cortes_jenks[i], cortes_jenks[i + 1]))

# ==============================================================================
# OUTPUT 1: PLANIFICACIÓN ESPACIAL GENERAL (4 paneles)
# Capas: DEM | Temperatura | Precipitación | Suelos SAL250
# ==============================================================================
cat("\n--- Generando Output 1: Planificación Espacial General ---\n")

png("planificacion_espacial_general.png", width = 3200, height = 2400, res = 250)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 5))

# Panel 1 - Elevación
plot(dem_puno_final,
     main = "1. Elevación (m.s.n.m.) - DEM NASA",
     col  = colorRampPalette(c("#3288bd","#99d594","#ffffbf","#fc8d59","#d53e4f"))(100))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

# Panel 2 - Temperatura
plot(raster_temperatura,
     main = "2. Temperatura Media Estimada (°C)\nNASA POWER + Gradiente Adiabático",
     col  = colorRampPalette(c("#313695","#74add1","#fee090","#f46d43","#d73027"))(100))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

# Panel 3 - Precipitación
plot(raster_precip,
     main = "3. Precipitación Estimada (mm/año)\nNASA POWER + Modelo Orográfico",
     col  = colorRampPalette(c("#ffffd9","#41b6c4","#225ea8","#081d58"))(100))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

# Panel 4 - Suelos SAL250 (vectorial directo, sin rasterizar)
n_clases_suelo <- length(unique(suelos_puno_utm$id_suelo))
plot(shp_puno_utm,
     main = "4. Capacidad de Uso Mayor de Suelos\n(SAL250 - ONERN)",
     col  = "white", border = "gray70")
plot(suelos_puno_utm, "id_suelo", add = TRUE,
     col = rainbow(n_clases_suelo, alpha = 0.75))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

dev.off()
cat("Output 1 guardado: planificacion_espacial_general.png\n")

# ==============================================================================
# OUTPUT 2: PLANIFICACIÓN DE PRODUCTO - QUINUA (4 paneles)
# Capas: Rendimiento Jenks | Zona térmica | Zona precipitación | Riesgo heladas
# ==============================================================================
cat("--- Generando Output 2: Planificación de Producto QUINUA ---\n")

paleta_q <- colorRampPalette(c("#FEF0D9","#FDCC8A","#FC8D59",
                               "#E34A33","#B30000"))(length(cortes_jenks) - 1)

png("planificacion_quinua_puno.png", width = 3200, height = 2400, res = 250)
par(mfrow = c(2, 2), mar = c(2, 2, 3, 5))

# Panel 1 - Rendimiento por distrito con Jenks
plot(mapa_quinua, "Rendimiento_Promedio",
     type   = "interval",
     breaks = cortes_jenks,
     col    = paleta_q,
     main   = "1. Rendimiento QUINUA por Distrito (T/Ha)\nClasificación Jenks · ENA 2025")
plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)

# Panel 2 - Zona térmica óptima para quinua
# FIX: Operación en RAM con pmax/pmin, sin aritmética raster
v_t      <- values(raster_temperatura, mat = FALSE)
v_zona_t <- ifelse(v_t >= 8  & v_t <= 12, 1L,
                   ifelse(v_t >= 6  & v_t <  8,  2L,
                          ifelse(v_t > 12  & v_t <= 14, 3L, NA_integer_)))
zona_termica        <- raster_temperatura
values(zona_termica) <- v_zona_t
names(zona_termica)  <- "Zona_Termica"
rm(v_t, v_zona_t); gc()

plot(zona_termica,
     main   = "2. Zonas Térmicas para QUINUA",
     col    = c("#2c7bb6","#abd9e9","#fdae61"),
     type   = "classes",
     levels = c("Óptima 8-12°C","Sub-óptima fría 6-8°C","Sub-óptima cálida 12-14°C"))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

# Panel 3 - Zona de precipitación para quinua
v_p      <- values(raster_precip, mat = FALSE)
v_zona_p <- ifelse(v_p >= 400 & v_p <= 700, 1L,
                   ifelse(v_p <  400,              2L,
                          ifelse(v_p >  700,              3L, NA_integer_)))
zona_precip        <- raster_precip
values(zona_precip) <- v_zona_p
names(zona_precip)  <- "Zona_Precip"
rm(v_p, v_zona_p); gc()

plot(zona_precip,
     main   = "3. Zonas de Precipitación para QUINUA",
     col    = c("#1a9641","#fdae61","#d7191c"),
     type   = "classes",
     levels = c("Óptima 400-700mm","Déficit <400mm","Exceso >700mm"))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

# Panel 4 - Riesgo por heladas con Jenks (desde ENA 2025)
vals_h <- mapa_quinua$Porcentaje_Heladas
vals_h <- vals_h[!is.na(vals_h)]

if (length(unique(vals_h)) >= 3) {
  cortes_h <- unique(classIntervals(vals_h,
                                    n = min(5, length(unique(vals_h)) - 1),
                                    style = "jenks")$brks)
  paleta_h <- colorRampPalette(c("#f7fbff","#6baed6","#08306b"))(length(cortes_h) - 1)
  plot(mapa_quinua, "Porcentaje_Heladas",
       type   = "interval",
       breaks = cortes_h,
       col    = paleta_h,
       main   = "4. Riesgo por Heladas por Distrito (%)\nClasificación Jenks · ENA 2025")
  cat("\nCortes Jenks - Riesgo Heladas (%):\n")
  for (i in seq_len(length(cortes_h) - 1))
    cat(sprintf("  Clase %d: %.1f%% - %.1f%%\n", i, cortes_h[i], cortes_h[i + 1]))
} else {
  plot(mapa_quinua, "Porcentaje_Sequia",
       col  = colorRampPalette(c("#fff5eb","#fd8d3c","#7f2704"))(100),
       main = "4. Riesgo por Sequía por Distrito (%)\nENA 2025")
}
plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)

dev.off()
cat("Output 2 guardado: planificacion_quinua_puno.png\n")

# ==============================================================================
# RESUMEN FINAL
# ==============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════╗\n")
cat("║     PLANIFICACIÓN AGROCLIMÁTICA PUNO - COMPLETA      ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Temp. base NASA POWER   : %5.2f °C               ║\n", temp_base_nasa))
cat(sprintf("║  Precipitación media     : %5.1f mm/año           ║\n", precip_base_nasa))
cat(sprintf("║  Distritos quinueros     : %3d                    ║\n", nrow(mapa_quinua)))
cat(sprintf("║  Clases Jenks aplicadas  : %3d                    ║\n", length(cortes_jenks)-1))
cat("╠══════════════════════════════════════════════════════╣\n")
cat("║  OUTPUT 1 → planificacion_espacial_general.png      ║\n")
cat("║    Panel 1: Elevación DEM                           ║\n")
cat("║    Panel 2: Temperatura (NASA POWER + gradiente)    ║\n")
cat("║    Panel 3: Precipitación (NASA POWER + orográfico) ║\n")
cat("║    Panel 4: Suelos SAL250 (ONERN)                   ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
cat("║  OUTPUT 2 → planificacion_quinua_puno.png           ║\n")
cat("║    Panel 1: Rendimiento Jenks (ENA 2025)            ║\n")
cat("║    Panel 2: Zonas térmicas óptimas                  ║\n")
cat("║    Panel 3: Zonas de precipitación                  ║\n")
cat("║    Panel 4: Riesgo por heladas Jenks (ENA 2025)     ║\n")
cat("╚══════════════════════════════════════════════════════╝\n")