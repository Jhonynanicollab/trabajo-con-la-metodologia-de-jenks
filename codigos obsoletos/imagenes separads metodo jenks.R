# ==============================================================================
# PROYECTO: APTITUD AGROCLIMÁTICA - REGIÓN PUNO
# OBJETIVO: Dos Planificaciones → (1) Espacial General | (2) Producto: Quinua
# Metodología: Jenks Natural Breaks | NASA POWER | SENAMHI | SAL250 | ENA 2025
# VERSIÓN CORREGIDA Y MEJORADA - Integración real de NASA POWER y SENAMHI
# ==============================================================================

library(terra)
library(sf)
library(data.table)
library(classInt)

setwd("E:/estadistica espacial")

# ==============================================================================
# FUNCIÓN UTILITARIA: normalizar texto (sin tildes, sin espacios raros)
# Usada en múltiples bloques, se define aquí arriba para disponibilidad global
# ==============================================================================
normalizar_texto <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  gsub("[^A-Z0-9]", "", x)
}

# ==============================================================================
# BLOQUE 1: MARCO GEOGRÁFICO
# Sin cambios — estaba correcto
# ==============================================================================
cat("=== BLOQUE 1: Marco geográfico ===\n")

shp_peru_vect  <- vect("gadm41_PER_3.shp")
crs(shp_peru_vect) <- "EPSG:4326"
shp_puno_wgs84 <- shp_peru_vect[shp_peru_vect$NAME_1 == "Puno", ]

crs_wgs84  <- "+proj=longlat +datum=WGS84 +no_defs"
crs_utm19s <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs"
crs(shp_puno_wgs84) <- crs_wgs84

shp_puno_utm <- project(shp_puno_wgs84, crs_utm19s)
writeVector(shp_puno_utm, "distritos_puno_utm.shp", overwrite = TRUE)
cat("Distritos de Puno cargados:", nrow(shp_puno_utm), "\n")

# ==============================================================================
# BLOQUE 2: DEM — UNIFICAR Y REPROYECTAR
# Sin cambios — fix ya estaba bien aplicado
# ==============================================================================
cat("\n=== BLOQUE 2: DEM NASA ===\n")

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

# FIX APLICADO: reproyectar a UTM 19S para tener todo en metros
dem_puno_final <- project(dem_puno_wgs84, crs_utm19s)
names(dem_puno_final) <- "Elevacion_m"
cat("DEM reproyectado a UTM 19S. Rango altitudinal:",
    round(minmax(dem_puno_final)[1]), "-", round(minmax(dem_puno_final)[2]), "m\n")

# ==============================================================================
# BLOQUE 3A: NASA POWER — LECTURA Y PROCESAMIENTO REAL
#
# CORRECCIÓN COMPLETA: El archivo original de NASA POWER viene en formato wide
# (PARAMETRO, YEAR, JAN, FEB, ... DEC, ANN). Tu CSV tidy ya tiene las columnas:
# DATE, YEAR, MONTH, MONTH_NAME, solar_rad_MJm2day, precip_mmday, humidity_pct, temp_C
#
# MEJORA: Calculamos promedios mensuales multianuales para el modelo espacial
# y también extraemos estadísticas de variabilidad (desv. estándar) útiles para
# cuantificar la incertidumbre del gradiente térmico.
# ==============================================================================
cat("\n=== BLOQUE 3A: NASA POWER ===\n")

nasa_raw <- fread("nasa_power_tidy.csv")

# Diagnóstico de columnas — para detectar nombres incorrectos antes de fallar
cat("Columnas en nasa_power_tidy.csv:", paste(names(nasa_raw), collapse = ", "), "\n")
cols_requeridas <- c("temp_C", "precip_mmday", "humidity_pct", "solar_rad_MJm2day")
cols_faltantes  <- setdiff(cols_requeridas, names(nasa_raw))
if (length(cols_faltantes) > 0) {
  stop(paste("Faltan columnas en nasa_power_tidy.csv:", paste(cols_faltantes, collapse = ", "),
             "\nVerifica que el CSV tenga exactamente esos nombres."))
}

# Promedios anuales para el modelo espacial
temp_base_nasa    <- mean(nasa_raw$temp_C,            na.rm = TRUE)
precip_base_nasa  <- mean(nasa_raw$precip_mmday,      na.rm = TRUE) * 365  # mm/día → mm/año
humedad_media     <- mean(nasa_raw$humidity_pct,      na.rm = TRUE)
solar_medio       <- mean(nasa_raw$solar_rad_MJm2day, na.rm = TRUE)
altitud_base_nasa <- 4105.46  # Elevación MERRA-2 de la celda (metadato original NASA)

cat(sprintf("Temperatura base      : %.2f °C\n", temp_base_nasa))
cat(sprintf("Precipitación media   : %.1f mm/año\n", precip_base_nasa))
cat(sprintf("Humedad relativa media: %.1f %%\n", humedad_media))
cat(sprintf("Radiación solar media : %.2f MJ/m²/día\n", solar_medio))
cat(sprintf("Altitud referencia    : %.2f m.s.n.m.\n", altitud_base_nasa))

# Promedios mensuales (útiles para tablas del informe)
nasa_mensual <- nasa_raw[, .(
  Temp_Media_C    = mean(temp_C,            na.rm = TRUE),
  Precip_mmdia    = mean(precip_mmday,      na.rm = TRUE),
  Humedad_pct     = mean(humidity_pct,      na.rm = TRUE),
  Radiacion_MJm2  = mean(solar_rad_MJm2day, na.rm = TRUE)
), by = .(MONTH, MONTH_NAME)]
setorder(nasa_mensual, MONTH)

cat("\nResumen mensual NASA POWER (promedio 2020-2024):\n")
print(nasa_mensual)
fwrite(nasa_mensual, "nasa_power_resumen_mensual.csv")

# ==============================================================================
# BLOQUE 3B: SENAMHI — LECTURA Y PROCESAMIENTO
#
# CORRECCIÓN: La data SENAMHI tiene columnas:
# AÑO, MES, DIA, PRECIPITACION, TEM_MAX, TEMP_MIN, Temp_Prom
#
# MEJORA: Calculamos promedios anuales y mensuales, y derivamos indicadores
# críticos para aptitud agrícola:
#   - Días con temperatura mínima < 0°C → riesgo de helada
#   - Precipitación total anual
#   - Temperatura media de la estación de cultivo (oct-mar para altiplano)
# ==============================================================================
cat("\n=== BLOQUE 3B: SENAMHI ===\n")

senamhi_raw <- fread("senamhi_puno.csv")  # ajusta el nombre si es diferente

cat("Columnas en senamhi_puno.csv:", paste(names(senamhi_raw), collapse = ", "), "\n")
cat("Filas totales:", nrow(senamhi_raw), "\n")

# Limpieza: reemplazar -999 o valores faltantes de SENAMHI
senamhi_raw[PRECIPITACION < 0,  PRECIPITACION := NA]
senamhi_raw[TEM_MAX       < -50, TEM_MAX       := NA]
senamhi_raw[TEMP_MIN      < -50, TEMP_MIN      := NA]
senamhi_raw[Temp_Prom     < -50, Temp_Prom     := NA]

# MEJORA: Indicadores derivados por día
senamhi_raw[, Dia_Helada    := ifelse(!is.na(TEMP_MIN) & TEMP_MIN < 0, 1L, 0L)]
senamhi_raw[, Dia_Sequia    := ifelse(!is.na(PRECIPITACION) & PRECIPITACION == 0, 1L, 0L)]

# Resumen anual SENAMHI
senamhi_anual <- senamhi_raw[, .(
  Precip_Total_mm   = sum(PRECIPITACION, na.rm = TRUE),
  Temp_Media_C      = mean(Temp_Prom,    na.rm = TRUE),
  Temp_Max_Media    = mean(TEM_MAX,      na.rm = TRUE),
  Temp_Min_Media    = mean(TEMP_MIN,     na.rm = TRUE),
  Dias_Helada       = sum(Dia_Helada,    na.rm = TRUE),
  Dias_Sin_Lluvia   = sum(Dia_Sequia,    na.rm = TRUE)
), by = ANIO]
setorder(senamhi_anual, ANIO)

cat("\nResumen anual SENAMHI:\n")
print(senamhi_anual)
fwrite(senamhi_anual, "senamhi_resumen_anual.csv")

# Resumen mensual SENAMHI (para comparar con NASA POWER)
senamhi_mensual <- senamhi_raw[, .(
  Precip_Media_mm  = mean(PRECIPITACION, na.rm = TRUE),
  Temp_Media_C     = mean(Temp_Prom,     na.rm = TRUE),
  Dias_Helada_Prom = mean(Dia_Helada,    na.rm = TRUE) * 30
), by = MES]
setorder(senamhi_mensual, MES)
fwrite(senamhi_mensual, "senamhi_resumen_mensual.csv")

# ==============================================================================
# BLOQUE 3C: COMPARACIÓN Y FUSIÓN NASA POWER + SENAMHI
#
# MEJORA: Creamos una tabla de comparación mensual entre ambas fuentes.
# Usamos el promedio ponderado (60% SENAMHI estación real, 40% NASA POWER)
# como valor base para el modelo espacial, dado que SENAMHI es observación
# directa y NASA POWER es reanálisis de modelo atmosférico.
# ==============================================================================
cat("\n=== BLOQUE 3C: Fusión NASA POWER + SENAMHI ===\n")

# Alinear por mes
fusion_mensual <- merge(
  nasa_mensual[, .(MES = MONTH, Temp_NASA = Temp_Media_C, Precip_NASA = Precip_mmdia * 30)],
  senamhi_mensual[, .(MES, Temp_SENAMHI = Temp_Media_C, Precip_SENAMHI = Precip_Media_mm * 30)],
  by = "MES", all = TRUE
)

# Valor fusionado: prioriza SENAMHI cuando existe
fusion_mensual[, Temp_Fusion   := ifelse(!is.na(Temp_SENAMHI),
                                         0.6 * Temp_SENAMHI + 0.4 * Temp_NASA,
                                         Temp_NASA)]
fusion_mensual[, Precip_Fusion := ifelse(!is.na(Precip_SENAMHI),
                                         0.6 * Precip_SENAMHI + 0.4 * Precip_NASA,
                                         Precip_NASA)]

cat("Tabla de comparación y fusión mensual NASA-SENAMHI:\n")
print(fusion_mensual)
fwrite(fusion_mensual, "fusion_nasa_senamhi_mensual.csv")

# Valores anuales fusionados para el modelo espacial
temp_base_fusion   <- mean(fusion_mensual$Temp_Fusion,   na.rm = TRUE)
precip_base_fusion <- sum(fusion_mensual$Precip_Fusion,  na.rm = TRUE)  # mm/año

# Días de helada promedio anual desde SENAMHI (dato más confiable)
dias_helada_anual <- mean(senamhi_anual$Dias_Helada, na.rm = TRUE)

cat(sprintf("\nValores fusionados para modelo espacial:\n"))
cat(sprintf("  Temperatura base  : %.2f °C\n",     temp_base_fusion))
cat(sprintf("  Precipitación base: %.1f mm/año\n", precip_base_fusion))
cat(sprintf("  Días helada/año   : %.1f días\n",   dias_helada_anual))

# ==============================================================================
# BLOQUE 4: MODELADO ESPACIAL DEL CLIMA
#
# CORRECCIÓN: Ahora usa temp_base_fusion y precip_base_fusion (fusión real)
# en lugar de solo NASA POWER.
# MEJORA: Se añade el raster de riesgo de helada derivado del DEM + umbral
# de temperatura mínima observada en SENAMHI.
# ==============================================================================
cat("\n=== BLOQUE 4: Rasters climáticos espaciales ===\n")

v_dem <- values(dem_puno_final, mat = FALSE)

# --- Raster de Temperatura ---
# Gradiente adiabático seco: -0.0065 °C por metro de altitud
v_temp <- temp_base_fusion - 0.0065 * (v_dem - altitud_base_nasa)
v_temp <- pmax(pmin(v_temp, 15), -2)  # Clamp agrícola Altiplano [-2, 15] °C

raster_temperatura        <- dem_puno_final
values(raster_temperatura) <- v_temp
names(raster_temperatura)  <- "Temperatura_Media_C"
writeRaster(raster_temperatura, "raster_temperatura_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)

# --- Raster de Precipitación ---
# Modelo orográfico simple: precipitación aumenta con altitud hasta cierto punto
# y luego decrece (efecto sombra de lluvia sobre los 4500 m en el altiplano)
delta_alt <- v_dem - altitud_base_nasa
v_prec    <- precip_base_fusion + 0.05 * delta_alt - 0.000005 * delta_alt^2
v_prec    <- pmax(pmin(v_prec, 900), 200)  # Clamp [200, 900] mm/año

raster_precip        <- dem_puno_final
values(raster_precip) <- v_prec
names(raster_precip)  <- "Precipitacion_mm"
writeRaster(raster_precip, "raster_precipitacion_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)

###################

# Forma correcta de leer bloques en terra
raster_helada <- dem_puno_final
names(raster_helada) <- "Dias_Helada_Estimados"

n_celdas   <- ncell(dem_puno_final)
tam_bloque <- 500000
n_bloques  <- ceiling(n_celdas / tam_bloque)

cat("Procesando raster de helada en", n_bloques, "bloques...\n")

v_resultado <- rep(NA_real_, n_celdas)

for (i in seq_len(n_bloques)) {
  inicio <- (i - 1) * tam_bloque + 1
  fin    <- min(i * tam_bloque, n_celdas)
  
  # CORRECCIÓN: leer todos los valores una sola vez fuera del loop
  # y subsetear el vector, no el raster
  if (i == 1) {
    cat("  Leyendo valores del DEM completo en RAM...\n")
    v_dem_completo <- as.vector(values(dem_puno_final))
  }
  
  bloque_dem  <- v_dem_completo[inicio:fin]
  bloque_tmin <- (temp_base_fusion - 0.0065 * (bloque_dem - altitud_base_nasa)) - 8
  
  v_resultado[inicio:fin] <- ifelse(
    is.na(bloque_dem), NA_real_,
    ifelse(bloque_tmin < 0,
           pmin(dias_helada_anual * (1 + abs(bloque_tmin) / 3), 200),
           pmax(dias_helada_anual * (1 - abs(bloque_tmin) / 5), 0))
  )
  
  if (i %% 20 == 0) cat("  Bloque", i, "de", n_bloques, "\n")
  rm(bloque_dem, bloque_tmin)
}

rm(v_dem_completo); gc()

values(raster_helada) <- v_resultado
rm(v_resultado); gc()

writeRaster(raster_helada, "raster_helada_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)
cat("Raster de helada guardado.\n")


# ==============================================================================
# BLOQUE 5: SUELOS SAL250
# Sin cambios significativos — estaba correcto
# ==============================================================================
cat("\n=== BLOQUE 5: Suelos SAL250 ===\n")

if (file.exists("suelos_puno_utm.shp")) {
  suelos_puno_utm <- vect("suelos_puno_utm.shp")
} else {
  suelos_peru     <- vect("Capacidad_de_Uso_Mayor_ONERN_CUM_geogpsperu_SuyoPomalia_931381206.shp")
  suelos_puno_utm <- project(mask(crop(suelos_peru, shp_puno_wgs84), shp_puno_wgs84), crs_utm19s)
  writeVector(suelos_puno_utm, "suelos_puno_utm.shp", overwrite = TRUE)
}

# Recrear id_suelo en RAM (necesario porque el .shp puede no guardarlo)
columna_suelo            <- names(suelos_puno_utm)[1]
suelos_puno_utm$id_suelo <- as.numeric(as.factor(
  values(suelos_puno_utm)[, columna_suelo]))

cat("Columna de suelo usada:", columna_suelo,
    "| Clases únicas:", length(unique(suelos_puno_utm$id_suelo)), "\n")

# ==============================================================================
# BLOQUE 6: ENA 2025 — LIMPIEZA Y CRUCE DISTRITAL
#
# CORRECCIÓN COMPLETA:
# - Eliminado cruce defectuoso por UBIGEO
# - Uso de nombres normalizados para máxima compatibilidad
# - Uso de P1207_NOM como nombre territorial
# - Diagnóstico automático de coincidencias espaciales
# ==============================================================================

cat("\n=== BLOQUE 6: ENA 2025 ===\n")

# ----------------------------------------------------------------------
# CARGA Y LIMPIEZA
# ----------------------------------------------------------------------

data_agro <- fread("Data_Agro_puno.csv")

data_limpia <- data_agro[
  !is.na(P204_COD) &
    (!is.na(P210_SUP_1) | !is.na(P219_CANT_1))
]

# Filtrar solo Puno
data_limpia <- data_limpia[CCDD == 21]

cat("Filas ENA después de filtrar Puno:", nrow(data_limpia), "\n")

# ----------------------------------------------------------------------
# VARIABLES DERIVADAS
# ----------------------------------------------------------------------

data_limpia[, Afectado_Helada := ifelse(P223B_3 == 1, 1L, 0L)]
data_limpia[, Afectado_Sequia := ifelse(P223B_6 == 1, 1L, 0L)]

# Normalizar cultivo
data_limpia[, P204_NOM := normalizar_texto(P204_NOM)]

# ----------------------------------------------------------------------
# USAR P1207_NOM COMO IDENTIFICADOR TERRITORIAL
# ----------------------------------------------------------------------

cat("\nMuestra de nombres territoriales (P1207_NOM):\n")
print(head(unique(data_limpia$P1207_NOM), 20))

# Crear nombre limpio
data_limpia[, NOMBRE_DIST_CLEAN := normalizar_texto(P1207_NOM)]

# Verificar cantidad de nombres únicos
cat("\nCantidad de territorios únicos:",
    length(unique(data_limpia$NOMBRE_DIST_CLEAN)), "\n")

# ----------------------------------------------------------------------
# AGREGACIÓN ESTADÍSTICA POR DISTRITO Y CULTIVO
# ----------------------------------------------------------------------

distritos_rendimiento <- data_limpia[, .(
  
  Total_Productores    = .N,
  
  Sup_Cosechada_Ha     = sum(P217_SUP_1,  na.rm = TRUE),
  
  Produccion_Total_T   = sum(P219_CANT_1, na.rm = TRUE),
  
  Rendimiento_Promedio =
    mean(P219_CANT_1 / (P217_SUP_1 + 0.001), na.rm = TRUE),
  
  Porcentaje_Heladas =
    mean(Afectado_Helada, na.rm = TRUE) * 100,
  
  Porcentaje_Sequia =
    mean(Afectado_Sequia, na.rm = TRUE) * 100
  
), by = .(NOMBRE_DIST_CLEAN, P204_NOM)]

# ----------------------------------------------------------------------
# DIAGNÓSTICO DE CULTIVOS
# ----------------------------------------------------------------------

cat("\nCultivos únicos encontrados:\n")

print(
  distritos_rendimiento[
    ,
    .N,
    by = P204_NOM
  ][order(-N)]
)

# ----------------------------------------------------------------------
# PREPARAR SHAPEFILE
# ----------------------------------------------------------------------

shp_puno_utm <- vect("distritos_puno_utm.shp")

# Crear nombre limpio en shapefile
shp_puno_utm$NAME_CLEAN <-
  normalizar_texto(shp_puno_utm$NAME_3)

cat("\nDistritos en shapefile:", nrow(shp_puno_utm), "\n")

# ----------------------------------------------------------------------
# FILTRAR SOLO QUINUA
# ----------------------------------------------------------------------

tabla_quinua <-
  distritos_rendimiento[
    grepl("QUINUA", P204_NOM)
  ]

cat("\nRegistros de QUINUA encontrados:",
    nrow(tabla_quinua), "\n")

# ----------------------------------------------------------------------
# CRUCE ESPACIAL POR NOMBRE NORMALIZADO
# ----------------------------------------------------------------------

mapa_quinua <- merge(
  shp_puno_utm,
  tabla_quinua,
  by.x = "NAME_CLEAN",
  by.y = "NOMBRE_DIST_CLEAN",
  all.x = FALSE
)

cat("Distritos con coincidencia espacial:",
    nrow(mapa_quinua), "\n")

# ----------------------------------------------------------------------
# DIAGNÓSTICO AUTOMÁTICO
# ----------------------------------------------------------------------

if (nrow(mapa_quinua) == 0) {
  
  cat("\n========================================\n")
  cat("DIAGNÓSTICO DE CRUCE ESPACIAL\n")
  cat("========================================\n")
  
  cat("\nNombres SHAPEFILE:\n")
  print(head(unique(shp_puno_utm$NAME_CLEAN), 30))
  
  cat("\nNombres ENA:\n")
  print(head(unique(tabla_quinua$NOMBRE_DIST_CLEAN), 30))
  
  cat("\nCoincidencias exactas:\n")
  
  coincidencias <- intersect(
    unique(shp_puno_utm$NAME_CLEAN),
    unique(tabla_quinua$NOMBRE_DIST_CLEAN)
  )
  
  print(coincidencias)
  
  stop("No existen coincidencias espaciales.")
}

# ----------------------------------------------------------------------
# EXPORTAR SHAPEFILE FINAL
# ----------------------------------------------------------------------

writeVector(
  mapa_quinua,
  "mapa_distrital_quinua_puno.shp",
  overwrite = TRUE
)

cat("\nShapefile de quinua exportado correctamente.\n")
# ==============================================================================
# BLOQUE 7: JENKS — RENDIMIENTO Y HELADAS
#
# MEJORA: Guardamos las tablas de clasificación en CSV para el informe
# CORRECCIÓN MENOR: Verificación de n_cl >= 3 antes de graficar
# ==============================================================================
cat("\n=== BLOQUE 7: Clasificación Jenks ===\n")

# --- Jenks para Rendimiento ---
vals_rend <- mapa_quinua$Rendimiento_Promedio
vals_rend <- vals_rend[!is.na(vals_rend) & is.finite(vals_rend)]

if (length(unique(vals_rend)) < 3) {
  stop("Muy pocos valores únicos de rendimiento para Jenks. Verificar el cruce de datos.")
}

n_cl         <- min(5, length(unique(vals_rend)) - 1)
jenks_rend   <- classIntervals(vals_rend, n = n_cl, style = "jenks")
cortes_jenks <- unique(jenks_rend$brks)

etiq_rend <- c("Muy bajo", "Bajo", "Medio", "Alto", "Muy alto")[1:(length(cortes_jenks) - 1)]
paleta_q  <- colorRampPalette(c("#FEF0D9","#FDCC8A","#FC8D59","#E34A33","#B30000"))(length(cortes_jenks) - 1)

tabla_jenks_rend <- data.table(
  Clase         = seq_len(length(cortes_jenks) - 1),
  Etiqueta      = etiq_rend,
  Limite_Inf    = round(cortes_jenks[-length(cortes_jenks)], 2),
  Limite_Sup    = round(cortes_jenks[-1], 2)
)
cat("\nCortes Jenks — Rendimiento QUINUA (T/Ha):\n")
print(tabla_jenks_rend)
fwrite(tabla_jenks_rend, "jenks_rendimiento_quinua.csv")

# --- Jenks para Heladas ---
vals_h <- mapa_quinua$Porcentaje_Heladas
vals_h <- vals_h[!is.na(vals_h) & is.finite(vals_h)]

jenks_ok <- length(unique(vals_h)) >= 3
if (jenks_ok) {
  n_cl_h   <- min(5, length(unique(vals_h)) - 1)
  jenks_h  <- classIntervals(vals_h, n = n_cl_h, style = "jenks")
  cortes_h <- unique(jenks_h$brks)
  paleta_h <- colorRampPalette(c("#f7fbff","#6baed6","#08306b"))(length(cortes_h) - 1)
  
  tabla_jenks_h <- data.table(
    Clase      = seq_len(length(cortes_h) - 1),
    Limite_Inf = round(cortes_h[-length(cortes_h)], 1),
    Limite_Sup = round(cortes_h[-1], 1)
  )
  cat("\nCortes Jenks — Riesgo Heladas (%):\n")
  print(tabla_jenks_h)
  fwrite(tabla_jenks_h, "jenks_riesgo_heladas.csv")
} else {
  cat("AVISO: Insuficientes valores únicos de heladas. Panel 4 usará % sequía.\n")
}

# ==============================================================================
# BLOQUE 8: RASTERS DE ZONIFICACIÓN PARA QUINUA
# Los rasters de zonas térmicas y de precipitación se crean aquí,
# ANTES de los outputs, para que estén disponibles en ambos bloques.
# ==============================================================================
cat("\n=== BLOQUE 8: Zonificación para Quinua ===\n")

# Zona térmica óptima para Quinua (FAO / INIA Puno)
v_t      <- values(raster_temperatura, mat = FALSE)
v_zona_t <- ifelse(v_t >= 8  & v_t <= 12, 1L,
                   ifelse(v_t >= 6  & v_t <  8,  2L,
                          ifelse(v_t > 12  & v_t <= 14, 3L, NA_integer_)))
zona_termica        <- raster_temperatura
values(zona_termica) <- v_zona_t
names(zona_termica)  <- "Zona_Termica"

# Zona de precipitación para Quinua
v_p      <- values(raster_precip, mat = FALSE)
v_zona_p <- ifelse(v_p >= 400 & v_p <= 700, 1L,
                   ifelse(v_p <  400,              2L,
                          ifelse(v_p >  700,              3L, NA_integer_)))
zona_precip        <- raster_precip
values(zona_precip) <- v_zona_p
names(zona_precip)  <- "Zona_Precip"

rm(v_t, v_zona_t, v_p, v_zona_p); gc()
cat("Rasters de zonificación creados.\n")

# ==============================================================================
# FUNCIÓN AUXILIAR PARA GUARDAR PANEL INDIVIDUAL
# ==============================================================================
guardar_panel <- function(archivo, ancho = 1600, alto = 1400, res = 200, expr) {
  png(archivo, width = ancho, height = alto, res = res)
  par(mar = c(2, 2, 3.5, 5.5))
  force(expr)
  dev.off()
  cat("  Guardado:", archivo, "\n")
}

# Paletas y leyendas reutilizables
pal_elev  <- colorRampPalette(c("#3288bd","#99d594","#ffffbf","#fc8d59","#d53e4f"))(100)
pal_temp  <- colorRampPalette(c("#313695","#74add1","#fee090","#f46d43","#d73027"))(100)
pal_prec  <- colorRampPalette(c("#ffffd9","#41b6c4","#225ea8","#081d58"))(100)
col_zona_t <- c("#2c7bb6","#abd9e9","#fdae61")
col_zona_p <- c("#1a9641","#fdae61","#d7191c")
ley_zona_t <- c("Optima 8-12C","Sub-opt fria 6-8C","Sub-opt calida 12-14C")
ley_zona_p <- c("Optima 400-700mm","Deficit <400mm","Exceso >700mm")

# ==============================================================================
# OUTPUT 1: PLANIFICACIÓN ESPACIAL GENERAL — 4 paneles compuestos
# ==============================================================================
cat("\n=== OUTPUT 1: Planificación Espacial General ===\n")

png("planificacion_espacial_general.png", width = 3200, height = 2400, res = 250)
par(mfrow = c(2, 2), mar = c(2, 2, 3.5, 5.5))

plot(dem_puno_final, main = "1. Elevacion (m.s.n.m.) — DEM NASA",   col = pal_elev)
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

plot(raster_temperatura, main = "2. Temperatura Media (°C)\nNASA POWER + SENAMHI + Gradiente Adiabatico", col = pal_temp)
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

plot(raster_precip, main = "3. Precipitacion Estimada (mm/año)\nNASA POWER + SENAMHI + Modelo Orografico", col = pal_prec)
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

n_cl_s <- length(unique(suelos_puno_utm$id_suelo))
plot(shp_puno_utm, main = "4. Capacidad de Uso Mayor de Suelos\n(SAL250 — ONERN)", col = "white", border = "gray70")
plot(suelos_puno_utm, "id_suelo", add = TRUE, col = rainbow(n_cl_s, alpha = 0.75))
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

dev.off()
cat("Output 1 compuesto guardado: planificacion_espacial_general.png\n")

# ==============================================================================
# OUTPUT 2: PLANIFICACIÓN DE PRODUCTO — QUINUA — 4 paneles compuestos
# ==============================================================================
cat("\n=== OUTPUT 2: Planificación de Producto QUINUA ===\n")

png("planificacion_quinua_puno.png", width = 3200, height = 2400, res = 250)
par(mfrow = c(2, 2), mar = c(2, 2, 3.5, 5.5))

# Panel 1 — Rendimiento Jenks
plot(mapa_quinua, "Rendimiento_Promedio", type = "interval",
     breaks = cortes_jenks, col = paleta_q,
     main   = "1. Rendimiento QUINUA por Distrito (T/Ha)\nClasificacion Jenks · ENA 2025")
plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)

# Panel 2 — Zona térmica con leyenda correcta
plot(zona_termica, type = "classes", col = col_zona_t,
     main = "2. Zonas Termicas para QUINUA\n(NASA POWER + SENAMHI + DEM)")
legend("bottomleft", legend = ley_zona_t, fill = col_zona_t, cex = 0.65, bty = "n")
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

# Panel 3 — Zona de precipitación con leyenda correcta
plot(zona_precip, type = "classes", col = col_zona_p,
     main = "3. Zonas de Precipitacion para QUINUA\n(NASA POWER + SENAMHI + Modelo Orografico)")
legend("bottomleft", legend = ley_zona_p, fill = col_zona_p, cex = 0.65, bty = "n")
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

# Panel 4 — Riesgo heladas (con fallback a sequía)
if (jenks_ok) {
  plot(mapa_quinua, "Porcentaje_Heladas", type = "interval",
       breaks = cortes_h, col = paleta_h,
       main   = "4. Riesgo por Heladas por Distrito (%)\nClasificacion Jenks · ENA 2025")
} else {
  plot(mapa_quinua, "Porcentaje_Sequia",
       col  = colorRampPalette(c("#fff5eb","#fd8d3c","#7f2704"))(100),
       main = "4. Riesgo por Sequia por Distrito (%)\nENA 2025")
}
plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)

dev.off()
cat("Output 2 compuesto guardado: planificacion_quinua_puno.png\n")

# ==============================================================================
# IMÁGENES INDIVIDUALES PARA EL INFORME (8 PNG separados de alta calidad)
# Numerados para insertar directamente en el documento
# ==============================================================================
cat("\n=== Guardando 8 imagenes individuales para el informe ===\n")

# --- OUTPUT 1: 4 mapas espaciales ---
guardar_panel("img_01_elevacion.png", expr = {
  plot(dem_puno_final, main = "Elevacion (m.s.n.m.) — DEM NASA", col = pal_elev)
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)
})

guardar_panel("img_02_temperatura.png", expr = {
  plot(raster_temperatura,
       main = "Temperatura Media Estimada (°C)\nNASA POWER + SENAMHI + Gradiente Adiabatico",
       col  = pal_temp)
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)
})

guardar_panel("img_03_precipitacion.png", expr = {
  plot(raster_precip,
       main = "Precipitacion Estimada (mm/año)\nNASA POWER + SENAMHI + Modelo Orografico",
       col  = pal_prec)
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)
})

guardar_panel("img_04_suelos.png", expr = {
  n_cl_s <- length(unique(suelos_puno_utm$id_suelo))
  plot(shp_puno_utm,
       main = "Capacidad de Uso Mayor de Suelos (SAL250 — ONERN)",
       col  = "white", border = "gray70")
  plot(suelos_puno_utm, "id_suelo", add = TRUE, col = rainbow(n_cl_s, alpha = 0.75))
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)
})

# --- OUTPUT 2: 4 mapas de quinua ---
guardar_panel("img_05_rendimiento_jenks.png", expr = {
  plot(mapa_quinua, "Rendimiento_Promedio", type = "interval",
       breaks = cortes_jenks, col = paleta_q,
       main   = "Rendimiento QUINUA por Distrito (T/Ha)\nClasificacion Jenks · ENA 2025")
  plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)
  legend("bottomleft", legend = tabla_jenks_rend$Etiqueta,
         fill = paleta_q, title = "Clases Jenks", cex = 0.65, bty = "n")
})

guardar_panel("img_06_zonas_termicas.png", expr = {
  plot(zona_termica, type = "classes", col = col_zona_t,
       main = "Zonas Termicas para QUINUA\n(NASA POWER + SENAMHI + DEM)")
  legend("bottomleft", legend = ley_zona_t, fill = col_zona_t, cex = 0.65, bty = "n")
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)
})

guardar_panel("img_07_zonas_precipitacion.png", expr = {
  plot(zona_precip, type = "classes", col = col_zona_p,
       main = "Zonas de Precipitacion para QUINUA\n(NASA POWER + SENAMHI + Modelo Orografico)")
  legend("bottomleft", legend = ley_zona_p, fill = col_zona_p, cex = 0.65, bty = "n")
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)
})

guardar_panel("img_08_riesgo_heladas.png", expr = {
  if (jenks_ok) {
    plot(mapa_quinua, "Porcentaje_Heladas", type = "interval",
         breaks = cortes_h, col = paleta_h,
         main   = "Riesgo por Heladas por Distrito (%)\nClasificacion Jenks · ENA 2025")
    legend("bottomleft", legend = paste0(round(cortes_h[-length(cortes_h)], 1),
                                         "% - ", round(cortes_h[-1], 1), "%"),
           fill = paleta_h, title = "Clases Jenks", cex = 0.65, bty = "n")
  } else {
    plot(mapa_quinua, "Porcentaje_Sequia",
         col  = colorRampPalette(c("#fff5eb","#fd8d3c","#7f2704"))(100),
         main = "Riesgo por Sequia por Distrito (%)\nENA 2025")
  }
  plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)
})

# ==============================================================================
# RESUMEN FINAL EN CONSOLA
# ==============================================================================
cat("\n")
cat("=============================================================\n")
cat("   PLANIFICACION AGROCLIMATICA PUNO — PROCESO COMPLETADO    \n")
cat("=============================================================\n")
cat(sprintf("  Temperatura base (fusion)   : %5.2f °C\n",     temp_base_fusion))
cat(sprintf("  Precipitacion base (fusion) : %5.1f mm/año\n", precip_base_fusion))
cat(sprintf("  Dias helada promedio        : %5.1f dias/año\n",dias_helada_anual))
cat(sprintf("  Distritos quinueros mapeados: %3d\n",           nrow(mapa_quinua)))
cat(sprintf("  Clases Jenks rendimiento    : %3d\n",           length(cortes_jenks) - 1))
cat("-------------------------------------------------------------\n")
cat("  ARCHIVOS GENERADOS:\n")
cat("  Compuestos  : planificacion_espacial_general.png\n")
cat("               planificacion_quinua_puno.png\n")
cat("  Individuales: img_01_elevacion.png\n")
cat("               img_02_temperatura.png\n")
cat("               img_03_precipitacion.png\n")
cat("               img_04_suelos.png\n")
cat("               img_05_rendimiento_jenks.png\n")
cat("               img_06_zonas_termicas.png\n")
cat("               img_07_zonas_precipitacion.png\n")
cat("               img_08_riesgo_heladas.png\n")
cat("  Tablas CSV  : nasa_power_resumen_mensual.csv\n")
cat("               senamhi_resumen_anual.csv\n")
cat("               senamhi_resumen_mensual.csv\n")
cat("               fusion_nasa_senamhi_mensual.csv\n")
cat("               jenks_rendimiento_quinua.csv\n")
cat("               jenks_riesgo_heladas.csv\n")
cat("=============================================================\n")