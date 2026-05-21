
library(terra)
library(sf)
library(data.table)
library(classInt)

setwd("E:/estadistica espacial")

# ==============================================================================
# FUNCIÓN UTILITARIA: normalizar texto (sin tildes, sin espacios raros)
# ==============================================================================
normalizar_texto <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  gsub("[^A-Z0-9]", "", x)
}

# ==============================================================================
# BLOQUE 1: MARCO GEOGRÁFICO
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
# [C1] CORRECCIÓN: mosaic() en terra requiere argumento fun="mean" y
#      se pasa como lista con c(..., fun="mean") para que do.call funcione.
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
  
  # [C1] FIX: do.call con mosaic necesita fun= explícito en terra
  lista_dem    <- lapply(archivos_dem, rast)
  dem_completo <- do.call(mosaic, c(lista_dem, fun = "mean"))
  crs(dem_completo) <- "EPSG:4326"
  
  dem_puno_wgs84 <- crop(dem_completo, shp_puno_wgs84)
  dem_puno_wgs84 <- mask(dem_puno_wgs84, shp_puno_wgs84)
  writeRaster(dem_puno_wgs84, "dem_puno_unificado.tif",
              gdal = c("COMPRESS=LZW"), overwrite = TRUE)
}

dem_puno_final <- project(dem_puno_wgs84, crs_utm19s)
names(dem_puno_final) <- "Elevacion_m"
cat("DEM reproyectado a UTM 19S. Rango altitudinal:",
    round(minmax(dem_puno_final)[1]), "-", round(minmax(dem_puno_final)[2]), "m\n")

# ==============================================================================
# BLOQUE 3A: NASA POWER — LECTURA Y PROCESAMIENTO REAL
# ==============================================================================
cat("\n=== BLOQUE 3A: NASA POWER ===\n")

nasa_raw <- fread("nasa_power_tidy.csv")

cat("Columnas en nasa_power_tidy.csv:", paste(names(nasa_raw), collapse = ", "), "\n")
cols_requeridas <- c("temp_C", "precip_mmday", "humidity_pct", "solar_rad_MJm2day")
cols_faltantes  <- setdiff(cols_requeridas, names(nasa_raw))
if (length(cols_faltantes) > 0) {
  stop(paste("Faltan columnas en nasa_power_tidy.csv:", paste(cols_faltantes, collapse = ", "),
             "\nVerifica que el CSV tenga exactamente esos nombres."))
}

temp_base_nasa    <- mean(nasa_raw$temp_C,            na.rm = TRUE)
precip_base_nasa  <- mean(nasa_raw$precip_mmday,      na.rm = TRUE) * 365
humedad_media     <- mean(nasa_raw$humidity_pct,      na.rm = TRUE)
solar_medio       <- mean(nasa_raw$solar_rad_MJm2day, na.rm = TRUE)
altitud_base_nasa <- 4105.46

cat(sprintf("Temperatura base      : %.2f °C\n",       temp_base_nasa))
cat(sprintf("Precipitación media   : %.1f mm/año\n",   precip_base_nasa))
cat(sprintf("Humedad relativa media: %.1f %%\n",       humedad_media))
cat(sprintf("Radiación solar media : %.2f MJ/m²/día\n",solar_medio))
cat(sprintf("Altitud referencia    : %.2f m.s.n.m.\n", altitud_base_nasa))

nasa_mensual <- nasa_raw[, .(
  Temp_Media_C   = mean(temp_C,            na.rm = TRUE),
  Precip_mmdia   = mean(precip_mmday,      na.rm = TRUE),
  Humedad_pct    = mean(humidity_pct,      na.rm = TRUE),
  Radiacion_MJm2 = mean(solar_rad_MJm2day, na.rm = TRUE)
), by = .(MONTH, MONTH_NAME)]
setorder(nasa_mensual, MONTH)

cat("\nResumen mensual NASA POWER (promedio 2020-2024):\n")
print(nasa_mensual)
fwrite(nasa_mensual, "nasa_power_resumen_mensual.csv")

# ==============================================================================
# BLOQUE 3B: SENAMHI — LECTURA Y PROCESAMIENTO
# [C2] CORRECCIÓN: renombrar "AÑO" → "ANIO" para evitar problemas con tilde,
#      y asegurar que MES sea integer para el merge posterior.
# ==============================================================================
cat("\n=== BLOQUE 3B: SENAMHI ===\n")

senamhi_raw <- fread("senamhi_puno.csv")

cat("Columnas originales en senamhi_puno.csv:", paste(names(senamhi_raw), collapse = ", "), "\n")

# [C2] FIX: renombrar columna con tilde si existe
if ("AÑO" %in% names(senamhi_raw)) {
  setnames(senamhi_raw, "AÑO", "ANIO")
  cat("Columna 'AÑO' renombrada a 'ANIO'.\n")
} else if (!"ANIO" %in% names(senamhi_raw)) {
  stop("No se encontró columna de año ('AÑO' o 'ANIO') en senamhi_puno.csv")
}

# [C2] FIX: asegurar que MES es integer
senamhi_raw[, MES := as.integer(MES)]
cat("Tipo de MES después del fix:", class(senamhi_raw$MES), "\n")

cat("Filas totales:", nrow(senamhi_raw), "\n")

# Limpieza de valores inválidos
senamhi_raw[PRECIPITACION < 0,   PRECIPITACION := NA]
senamhi_raw[TEM_MAX       < -50, TEM_MAX       := NA]
senamhi_raw[TEMP_MIN      < -50, TEMP_MIN      := NA]
senamhi_raw[Temp_Prom     < -50, Temp_Prom     := NA]

senamhi_raw[, Dia_Helada := ifelse(!is.na(TEMP_MIN) & TEMP_MIN < 0, 1L, 0L)]
senamhi_raw[, Dia_Sequia := ifelse(!is.na(PRECIPITACION) & PRECIPITACION == 0, 1L, 0L)]

# Resumen anual
senamhi_anual <- senamhi_raw[, .(
  Precip_Total_mm = sum(PRECIPITACION, na.rm = TRUE),
  Temp_Media_C    = mean(Temp_Prom,    na.rm = TRUE),
  Temp_Max_Media  = mean(TEM_MAX,      na.rm = TRUE),
  Temp_Min_Media  = mean(TEMP_MIN,     na.rm = TRUE),
  Dias_Helada     = sum(Dia_Helada,    na.rm = TRUE),
  Dias_Sin_Lluvia = sum(Dia_Sequia,    na.rm = TRUE)
), by = ANIO]
setorder(senamhi_anual, ANIO)

cat("\nResumen anual SENAMHI:\n")
print(senamhi_anual)
fwrite(senamhi_anual, "senamhi_resumen_anual.csv")

# Resumen mensual
senamhi_mensual <- senamhi_raw[, .(
  Precip_Media_mm  = mean(PRECIPITACION, na.rm = TRUE),
  Temp_Media_C     = mean(Temp_Prom,     na.rm = TRUE),
  Dias_Helada_Prom = mean(Dia_Helada,    na.rm = TRUE) * 30
), by = MES]
setorder(senamhi_mensual, MES)
fwrite(senamhi_mensual, "senamhi_resumen_mensual.csv")

# ==============================================================================
# BLOQUE 3C: FUSIÓN NASA POWER + SENAMHI
# ==============================================================================
cat("\n=== BLOQUE 3C: Fusión NASA POWER + SENAMHI ===\n")

fusion_mensual <- merge(
  nasa_mensual[, .(MES = MONTH, Temp_NASA = Temp_Media_C, Precip_NASA = Precip_mmdia * 30)],
  senamhi_mensual[, .(MES, Temp_SENAMHI = Temp_Media_C, Precip_SENAMHI = Precip_Media_mm * 30)],
  by = "MES", all = TRUE
)

fusion_mensual[, Temp_Fusion   := ifelse(!is.na(Temp_SENAMHI),
                                         0.6 * Temp_SENAMHI + 0.4 * Temp_NASA,
                                         Temp_NASA)]
fusion_mensual[, Precip_Fusion := ifelse(!is.na(Precip_SENAMHI),
                                         0.6 * Precip_SENAMHI + 0.4 * Precip_NASA,
                                         Precip_NASA)]

cat("Tabla de comparación y fusión mensual NASA-SENAMHI:\n")
print(fusion_mensual)
fwrite(fusion_mensual, "fusion_nasa_senamhi_mensual.csv")

temp_base_fusion   <- mean(fusion_mensual$Temp_Fusion,  na.rm = TRUE)
precip_base_fusion <- sum(fusion_mensual$Precip_Fusion, na.rm = TRUE)
dias_helada_anual  <- mean(senamhi_anual$Dias_Helada,   na.rm = TRUE)

cat(sprintf("\nValores fusionados para modelo espacial:\n"))
cat(sprintf("  Temperatura base  : %.2f °C\n",     temp_base_fusion))
cat(sprintf("  Precipitación base: %.1f mm/año\n", precip_base_fusion))
cat(sprintf("  Días helada/año   : %.1f días\n",   dias_helada_anual))

# ==============================================================================
# BLOQUE 4: MODELADO ESPACIAL DEL CLIMA
# [C4] CORRECCIÓN: rm(v_dem) después de crear rasters de temperatura y
#      precipitación para liberar RAM antes del loop de helada.
# ==============================================================================
cat("\n=== BLOQUE 4: Rasters climáticos espaciales ===\n")

v_dem <- values(dem_puno_final, mat = FALSE)

# --- Raster de Temperatura ---
v_temp <- temp_base_fusion - 0.0065 * (v_dem - altitud_base_nasa)
v_temp <- pmax(pmin(v_temp, 15), -2)

raster_temperatura        <- dem_puno_final
values(raster_temperatura) <- v_temp
names(raster_temperatura)  <- "Temperatura_Media_C"
writeRaster(raster_temperatura, "raster_temperatura_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)

# --- Raster de Precipitación ---
delta_alt <- v_dem - altitud_base_nasa
v_prec    <- precip_base_fusion + 0.05 * delta_alt - 0.000005 * delta_alt^2
v_prec    <- pmax(pmin(v_prec, 900), 200)

raster_precip        <- dem_puno_final
values(raster_precip) <- v_prec
names(raster_precip)  <- "Precipitacion_mm"
writeRaster(raster_precip, "raster_precipitacion_puno.tif",
            gdal = c("COMPRESS=LZW"), overwrite = TRUE)

# [C4] FIX: liberar v_dem y vectores intermedios antes del loop de helada
rm(v_dem, v_temp, delta_alt, v_prec); gc()

# --- Raster de Helada (procesamiento por bloques) ---
raster_helada        <- dem_puno_final
names(raster_helada) <- "Dias_Helada_Estimados"

n_celdas   <- ncell(dem_puno_final)
tam_bloque <- 500000
n_bloques  <- ceiling(n_celdas / tam_bloque)

cat("Procesando raster de helada en", n_bloques, "bloques...\n")
v_resultado    <- rep(NA_real_, n_celdas)
v_dem_completo <- as.vector(values(dem_puno_final))  # lectura única fuera del loop

for (i in seq_len(n_bloques)) {
  inicio <- (i - 1) * tam_bloque + 1
  fin    <- min(i * tam_bloque, n_celdas)
  
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
# ==============================================================================
cat("\n=== BLOQUE 5: Suelos SAL250 ===\n")

if (file.exists("suelos_puno_utm.shp")) {
  suelos_puno_utm <- vect("suelos_puno_utm.shp")
} else {
  suelos_peru     <- vect("Capacidad_de_Uso_Mayor_ONERN_CUM_geogpsperu_SuyoPomalia_931381206.shp")
  suelos_puno_utm <- project(mask(crop(suelos_peru, shp_puno_wgs84), shp_puno_wgs84), crs_utm19s)
  writeVector(suelos_puno_utm, "suelos_puno_utm.shp", overwrite = TRUE)
}

columna_suelo            <- names(suelos_puno_utm)[1]
suelos_puno_utm$id_suelo <- as.numeric(as.factor(
  values(suelos_puno_utm)[, columna_suelo]))

cat("Columna de suelo usada:", columna_suelo,
    "| Clases únicas:", length(unique(suelos_puno_utm$id_suelo)), "\n")

# ==============================================================================
# BLOQUE 6: ENA 2025 — LIMPIEZA Y CRUCE DISTRITAL
# [C3] CORRECCIÓN: diagnóstico de CC_2/CC_3 antes de construir UBIGEO
# [C5] CORRECCIÓN: verificar columnas P223B_3 y P223B_6 antes de usarlas
# ==============================================================================
cat("\n=== BLOQUE 6: ENA 2025 ===\n")

# 1. CARGAR DATA
data_agro <- fread("Data_Agro_puno.csv")
cat("Filas originales:", nrow(data_agro), "\n")

# 2. FILTRAR PUNO Y REGISTROS VÁLIDOS
data_limpia <- data_agro[
  !is.na(P204_COD) & (!is.na(P210_SUP_1) | !is.na(P219_CANT_1))
]
data_limpia <- data_limpia[CCDD == 21]
cat("Filas ENA después de filtrar Puno:", nrow(data_limpia), "\n")

# 3. VARIABLES DERIVADAS
# [C5] FIX: verificar que las columnas de afectación existen antes de usarlas
if (!"P223B_3" %in% names(data_limpia)) {
  warning("Columna P223B_3 no encontrada. Afectado_Helada se asignará como 0.")
  data_limpia[, Afectado_Helada := 0L]
} else {
  data_limpia[, Afectado_Helada := ifelse(P223B_3 == 1, 1L, 0L)]
}

if (!"P223B_6" %in% names(data_limpia)) {
  warning("Columna P223B_6 no encontrada. Afectado_Sequia se asignará como 0.")
  data_limpia[, Afectado_Sequia := 0L]
} else {
  data_limpia[, Afectado_Sequia := ifelse(P223B_6 == 1, 1L, 0L)]
}

# 4. NORMALIZAR NOMBRE DE CULTIVO
data_limpia[, P204_NOM := normalizar_texto(P204_NOM)]

# 5. CREAR UBIGEO
data_limpia[, UBIGEO_DIST := paste0(
  sprintf("%02d", as.integer(CCDD)),
  sprintf("%02d", as.integer(CCPP)),
  sprintf("%02d", as.integer(CCDI))
)]
cat("\nEjemplo UBIGEO ENA:\n")
print(head(unique(data_limpia$UBIGEO_DIST), 10))

# 6. AGREGACIÓN ESTADÍSTICA
distritos_rendimiento <- data_limpia[, .(
  Total_Productores   = .N,
  Sup_Cosechada_Ha    = sum(P217_SUP_1,                           na.rm = TRUE),
  Produccion_Total_T  = sum(P219_CANT_1,                          na.rm = TRUE),
  Rendimiento_Promedio= mean(P219_CANT_1 / (P217_SUP_1 + 0.001), na.rm = TRUE),
  Porcentaje_Heladas  = mean(Afectado_Helada,                     na.rm = TRUE) * 100,
  Porcentaje_Sequia   = mean(Afectado_Sequia,                     na.rm = TRUE) * 100
), by = .(CCDD, CCPP, CCDI, UBIGEO_DIST, P204_NOM)]

setnames(distritos_rendimiento, "UBIGEO_DIST", "UBIGEO")

cat("\nColumnas de distritos_rendimiento:\n"); print(names(distritos_rendimiento))
cat("\nPrimeras filas:\n");                    print(head(distritos_rendimiento))
cat("\nCultivos detectados:\n")
print(distritos_rendimiento[, .N, by = P204_NOM][order(-N)])

# 7. FILTRAR QUINUA
tabla_quinua <- copy(distritos_rendimiento[grepl("QUINUA", P204_NOM)])
cat("\nFilas de QUINUA:", nrow(tabla_quinua), "\n")

# 8. PREPARAR SHAPEFILE
shp_puno_utm <- vect("distritos_puno_utm.shp")
cat("\nColumnas shapefile:\n"); print(names(shp_puno_utm))


# Construir UBIGEO desde GID_3: formato "PER.22.PP.DI_1"
gid3   <- as.character(shp_puno_utm$GID_3)
partes <- regmatches(gid3, gregexpr("[0-9]+", gid3))

# CORRECCIÓN: reemplazar el primer número (22 de GADM) por 21 (código INEI de Puno)
shp_puno_utm$UBIGEO <- sapply(partes, function(x) {
  if (length(x) >= 3) {
    paste0(
      "21",                                # forzar código INEI de Puno
      sprintf("%02d", as.integer(x[2])),   # provincia
      sprintf("%02d", as.integer(x[3]))    # distrito
    )
  } else NA_character_
})

cat("NAs en UBIGEO shapefile:", sum(is.na(shp_puno_utm$UBIGEO)), "\n")
cat("Muestra UBIGEO shapefile:\n")
print(head(unique(shp_puno_utm$UBIGEO), 20))

# 9. DIAGNÓSTICO DE UBIGEO
cat("\nUBIGEO ENA (muestra):\n");       print(head(unique(tabla_quinua$UBIGEO), 20))
cat("\nUBIGEO SHAPEFILE (muestra):\n"); print(head(unique(shp_puno_utm$UBIGEO), 20))

# 10. VERIFICAR COINCIDENCIAS
coincidencias <- intersect(unique(tabla_quinua$UBIGEO), unique(shp_puno_utm$UBIGEO))
cat("\nCantidad de coincidencias:", length(coincidencias), "\n")
print(head(coincidencias, 20))

if (length(coincidencias) == 0) {
  cat("\n========================================\n")
  cat("NO EXISTEN COINCIDENCIAS UBIGEO\n")
  cat("========================================\n")
  cat("Posibles causas:\n")
  cat("- CC_2 y CC_3 no son provincia/distrito\n")
  cat("- El shapefile usa otra codificación\n")
  cat("- CCDI no coincide con GADM\n")
  stop("Cruce espacial fallido.")
}

# 11. CRUCE ESPACIAL FINAL
mapa_quinua <- merge(shp_puno_utm, tabla_quinua, by = "UBIGEO", all.x = FALSE)
cat("\nDistritos con datos de QUINUA:", nrow(mapa_quinua), "\n")

if (nrow(mapa_quinua) == 0) stop("El merge final produjo 0 filas.")

writeVector(mapa_quinua, "mapa_distrital_quinua_puno.shp", overwrite = TRUE)
cat("\nShapefile final guardado correctamente.\n")
cat("\n========================================\n")
cat("CRUCE ESPACIAL COMPLETADO\n")
cat("========================================\n")
cat("Distritos mapeados :", nrow(mapa_quinua), "\n")
cat("Coincidencias      :", length(coincidencias), "\n")
cat("========================================\n")

# ==============================================================================
# BLOQUE 7: JENKS — RENDIMIENTO Y HELADAS
# ==============================================================================
cat("\n=== BLOQUE 7: Clasificación Jenks ===\n")

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
  Clase      = seq_len(length(cortes_jenks) - 1),
  Etiqueta   = etiq_rend,
  Limite_Inf = round(cortes_jenks[-length(cortes_jenks)], 2),
  Limite_Sup = round(cortes_jenks[-1], 2)
)
cat("\nCortes Jenks — Rendimiento QUINUA (T/Ha):\n"); print(tabla_jenks_rend)
fwrite(tabla_jenks_rend, "jenks_rendimiento_quinua.csv")

vals_h  <- mapa_quinua$Porcentaje_Heladas
vals_h  <- vals_h[!is.na(vals_h) & is.finite(vals_h)]
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
  cat("\nCortes Jenks — Riesgo Heladas (%):\n"); print(tabla_jenks_h)
  fwrite(tabla_jenks_h, "jenks_riesgo_heladas.csv")
} else {
  cat("AVISO: Insuficientes valores únicos de heladas. Panel 4 usará % sequía.\n")
}

# ==============================================================================
# BLOQUE 8: RASTERS DE ZONIFICACIÓN PARA QUINUA
# ==============================================================================
cat("\n=== BLOQUE 8: Zonificación para Quinua ===\n")

v_t      <- values(raster_temperatura, mat = FALSE)
v_zona_t <- ifelse(v_t >= 8  & v_t <= 12, 1L,
                   ifelse(v_t >= 6  & v_t <  8,  2L,
                          ifelse(v_t > 12  & v_t <= 14, 3L, NA_integer_)))
zona_termica        <- raster_temperatura
values(zona_termica) <- v_zona_t
names(zona_termica)  <- "Zona_Termica"

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
# [C6] CORRECCIÓN: tryCatch() para capturar errores por panel sin detener todo
# ==============================================================================
guardar_panel <- function(archivo, ancho = 1600, alto = 1400, res = 200, expr) {
  tryCatch({
    png(archivo, width = ancho, height = alto, res = res)
    par(mar = c(2, 2, 3.5, 5.5))
    force(expr)
    dev.off()
    cat("  Guardado:", archivo, "\n")
  }, error = function(e) {
    # Cerrar el dispositivo gráfico si quedó abierto
    if (dev.cur() > 1) dev.off()
    cat("  ERROR al guardar", archivo, ":", conditionMessage(e), "\n")
  })
}

# Paletas y leyendas reutilizables
pal_elev   <- colorRampPalette(c("#3288bd","#99d594","#ffffbf","#fc8d59","#d53e4f"))(100)
pal_temp   <- colorRampPalette(c("#313695","#74add1","#fee090","#f46d43","#d73027"))(100)
pal_prec   <- colorRampPalette(c("#ffffd9","#41b6c4","#225ea8","#081d58"))(100)
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

plot(dem_puno_final,    main = "1. Elevacion (m.s.n.m.) — DEM NASA",                              col = pal_elev)
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

plot(raster_temperatura, main = "2. Temperatura Media (°C)\nNASA POWER + SENAMHI + Gradiente",   col = pal_temp)
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

plot(raster_precip,      main = "3. Precipitacion Estimada (mm/año)\nNASA POWER + Modelo Orogr.", col = pal_prec)
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.4)

n_cl_s <- length(unique(suelos_puno_utm$id_suelo))
plot(shp_puno_utm, main = "4. Capacidad de Uso Mayor de Suelos\n(SAL250 — ONERN)",
     col = "white", border = "gray70")
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

plot(mapa_quinua, "Rendimiento_Promedio", type = "interval",
     breaks = cortes_jenks, col = paleta_q,
     main   = "1. Rendimiento QUINUA por Distrito (T/Ha)\nClasificacion Jenks · ENA 2025")
plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)

plot(zona_termica, type = "classes", col = col_zona_t,
     main = "2. Zonas Termicas para QUINUA\n(NASA POWER + SENAMHI + DEM)")
legend("bottomleft", legend = ley_zona_t, fill = col_zona_t, cex = 0.65, bty = "n")
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

plot(zona_precip, type = "classes", col = col_zona_p,
     main = "3. Zonas de Precipitacion para QUINUA\n(NASA POWER + SENAMHI + Modelo Orografico)")
legend("bottomleft", legend = ley_zona_p, fill = col_zona_p, cex = 0.65, bty = "n")
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)

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
# IMÁGENES INDIVIDUALES PARA EL INFORME (8 PNG separados)
# ==============================================================================
cat("\n=== Guardando 8 imágenes individuales para el informe ===\n")

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
  plot(shp_puno_utm, main = "Capacidad de Uso Mayor de Suelos (SAL250 — ONERN)",
       col = "white", border = "gray70")
  plot(suelos_puno_utm, "id_suelo", add = TRUE, col = rainbow(n_cl_s, alpha = 0.75))
  plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)
})

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
    legend("bottomleft",
           legend = paste0(round(cortes_h[-length(cortes_h)], 1),
                           "% - ", round(cortes_h[-1], 1), "%"),
           fill = paleta_h, title = "Clases Jenks", cex = 0.65, bty = "n")
  } else {
    plot(mapa_quinua, "Porcentaje_Sequia",
         col  = colorRampPalette(c("#fff5eb","#fd8d3c","#7f2704"))(100),
         main = "Riesgo por Sequia por Distrito (%)\nENA 2025")
  }
  plot(shp_puno_utm, add = TRUE, border = "gray30", lwd = 0.5)
})

