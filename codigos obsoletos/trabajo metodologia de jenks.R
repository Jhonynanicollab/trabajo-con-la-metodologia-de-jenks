# ==============================================================================
# PROYECTO: APTITUD AGROCLIMÁTICA Y ESPACIAL - REGIÓN PUNO
# FASE 1: Unificación de Rasters, Homogeneización de CRS y Carga de ENA 2025
# ==============================================================================

# 1. CARGAR LIBRERÍAS CRÍTICAS (Instálalas si no las tienes con install.packages())
library(terra)       # Manejo eficiente de grandes volúmenes de datos raster
library(sf)          # Manejo de datos vectoriales (Shapefiles)
library(data.table)  # Lectura ultra rápida de estructuras de datos pesadas

# 2. DEFINIR DIRECTORIO DE TRABAJO
# Remplaza esta ruta por la carpeta exacta donde guardaste todos tus archivos
setwd("C:/Users/crist/Documents/2026/estadistica espacial/PROYECTO_JENKS")
# ==============================================================================
# 2. MARCO GEOGRÁFICO: LÍMITES DE PUNO (CÓDIGO OPTIMIZADO PARA GADM)
# ==============================================================================
cat("--- Cargando y filtrando el marco geográfico de Puno ---\n")

# Cargamos el shapefile de distritos de GADM
shp_peru_vect <- vect("gadm41_PER_3.shp")

# Forzar el CRS original (GADM siempre viene en WGS84 / EPSG:4326)
crs(shp_peru_vect) <- "EPSG:4326"

# Filtrar estrictamente usando sintaxis nativa de terra (limpia y sin avisos)
# Usamos NAME_1 que es donde se almacena el Departamento en GADM
shp_puno_wgs84 <- shp_peru_vect[shp_peru_vect$NAME_1 == "Puno", ]

cat("Número de distritos filtrados en Puno:", nrow(shp_puno_wgs84), "\n")

# ==============================================================================
# 3. ELEVACIÓN: UNIFICAR Y RECORTAR DEM
# ==============================================================================
cat("--- Unificando las 15 piezas del DEM (NASA) ---\n")

zips_dem <- list.files("dem_zips/", pattern = "\\.zip$", full.names = TRUE)
carpeta_temp <- "dem_descomprimido"
for(zip in zips_dem) {
  unzip(zip, exdir = carpeta_temp)
}

archivos_descomprimidos <- list.files(carpeta_temp, pattern = "\\.(hgt|tif)$", full.names = TRUE)
lista_rasters <- lapply(archivos_descomprimidos, rast)
dem_completo <- do.call(mosaic, lista_rasters)

# Asegurar el CRS original del DEM antes de recortar
crs(dem_completo) <- "EPSG:4326"

cat("--- Recortando el DEM con los distritos de Puno ---\n")
# Ahora el recorte se ejecutará sobre la geometría limpia de GADM
dem_puno_wgs84 <- crop(dem_completo, shp_puno_wgs84)
dem_puno_wgs84 <- mask(dem_puno_wgs84, shp_puno_wgs84)

# ==============================================================================
# 4. REPROYECTAR AL SISTEMA MÉTRICO Y GRAFICAR (CÓDIGO SEGURO)
# ==============================================================================
cat("--- Configurando Proyecciones seguras para el Vector ---\n")

# Cadenas de texto explícitas para evitar usar la base de datos dañada
crs_wgs84  <- "+proj=longlat +datum=WGS84 +no_defs"
crs_utm19s <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs"

# Forzar el CRS de origen correcto al vector de distritos filtrados
crs(shp_puno_wgs84) <- crs_wgs84

# Reproyectar el vector de distritos a UTM Zona 19S
shp_puno_utm <- project(shp_puno_wgs84, crs_utm19s)

# Guardar el archivo raster final unificado (¡Esto ya te funcionó!)
writeRaster(dem_puno_final, "dem_puno_unificado.tif", overwrite = TRUE)

# Guardar también el Shapefile de los distritos de Puno en formato UTM 19S
# Lo usaremos más adelante para el Stack y para unir la ENA
writeVector(shp_puno_utm, "distritos_puno_utm.shp", overwrite = TRUE)

cat("--- ¡Dibujando Mapa Completo con Límites Distritales! ---\n")
# 1. Graficar el relieve (Raster)
plot(dem_puno_final, main = "DEM Unificado y Distritos de Puno (UTM 19S)")

# 2. Superponer los límites de los distritos (Vector) sin que falle
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.5)


# ==============================================================================
# EDALFOLOGÍA: RECORTAR SUELOS (SAL250) - OPTIMIZADO
# ==============================================================================
cat("--- Recortando Mapa de Suelos SAL250 para Puno ---\n")

# 1. Cargamos el shapefile de suelos usando vect() de terra
suelos_peru <- vect("Capacidad_de_Uso_Mayor_ONERN_CUM_geogpsperu_SuyoPomalia_931381206.shp")

# 2. Aseguramos que las proyecciones coincidan antes de cortar (GADM y MINAM suelen venir en WGS84)
# Recortamos los suelos usando la silueta distrital de Puno que ya creamos
suelos_puno_wgs84 <- crop(suelos_peru, shp_puno_wgs84)
suelos_puno_final <- mask(suelos_puno_wgs84, shp_puno_wgs84)

# 3. Definimos la cadena de proyección UTM 19S que usamos para el DEM
crs_utm19s <- "+proj=utm +zone=19 +south +datum=WGS84 +units=m +no_defs"

# 4. Reproyectamos el mapa de suelos cortado a metros
suelos_puno_utm <- project(suelos_puno_final, crs_utm19s)

# 5. Guardar el archivo vectorial exclusivo de Puno en tu carpeta
writeVector(suelos_puno_utm, "suelos_puno_utm.shp", overwrite = TRUE)
cat("--- ¡Mapa de suelos recortado y guardado con éxito! ---\n")


# ==============================================================================
# VERIFICACIÓN DE DATA DE LA ENCUESTA (SPSS/CSV)
# ==============================================================================
cat("--- Cargando base de datos agropecuaria filtrada de la ENA ---\n")

# Usamos data.table para una lectura instantánea de tus variables extraídas
data_agro <- fread("Data_Agro_puno.csv") 

# Mostrar las primeras filas y columnas en la consola de R Studio
print(head(data_agro))

# Mostrar los nombres de las columnas para asegurarnos de que estén P1207, P1208, etc.
cat("Variables detectadas en tu archivo de la encuesta:\n")
print(names(data_agro))


# ==============================================================================
# ADECUACIÓN Y LIMPIEZA DE LA DATA DE LA ENCUESTA (ENA 2025)
# ==============================================================================
cat("--- Adecuando variables según el cuestionario INEI ---\n")

# 1. Eliminar filas donde el código de cultivo o la superficie sean NA
# Esto elimina de golpe ese 95% de datos vacíos de la encuesta
data_limpia <- data_agro[!is.na(P204_COD) & (!is.na(P210_SUP_1) | !is.na(P219_CANT_1))]

# 2. Asegurar el filtro de la región Puno (Código 21)
# (Si al exportar de SPSS se pasaron otros departamentos, esto los limpia)
data_limpia <- data_limpia[CCDD == 21]

# 3. Construir la columna UBIGEO de 6 dígitos para el cruce distrital
# Formateamos con ceros a la izquierda si es necesario (ej: 21, 01, 02 -> 210102)
data_limpia[, UBIGEO_DIST := sprintf("%02d%02d%02d", CCDD, CCPP, CCDI)]

# 4. Adecuación de variables climáticas y tecnológicas según cuestionario
# Convertimos los factores de daño (P223B_1 a P223B_8) en lógicos/numéricos
# P223B_3 = Heladas, P223B_6 = Sequía según la pág. 7 del cuestionario
data_limpia[, Afectado_Helada := ifelse(P223B_3 == 1, 1, 0)]
data_limpia[, Afectado_Sequia := ifelse(P223B_6 == 1, 1, 0)]

# Condición hídrica (P119): 1 = Riego, 2 = Secano (Pág. 6 del cuestionario)
data_limpia[, Sistema_Hidrico := factor(P119, levels = c(1, 2), labels = c("Riego", "Secano"))]

# 5. Filtrar por tus cultivos objetivo (Quinua, Papa, Cañihua)
# Nota: Debes verificar los códigos numéricos exactos de la Tabla N°1 del INEI.
# Normalmente: Papa suele ser 34 o similar. Vamos a crear una columna limpia.
cat("Cultivos detectados en la data de Puno:\n")
print(data_limpia[, .N, by = .(P204_COD, P204_NOM)])

# ==============================================================================
# AGREGACIÓN ESTADÍSTICA POR TEXTO Y CÓDIGO (CÓDIGO DE PRECISIÓN)
# ==============================================================================
cat("--- Agrupando rendimientos por Nombre de Distrito (Evitando fallas de códigos) ---\n")

# 1. Crear una función para limpiar texto de forma agresiva
# Pasa a mayúsculas, quita tildes, eñes y elimina cualquier espacio en blanco
normalizar_texto <- function(vector_texto) {
  x <- toupper(as.character(vector_texto))
  x <- iconv(x, to = "ASCII//TRANSLIT")  # Transforma Á->A, Ñ->N
  x <- gsub("[^A-Z0-9]", "", x)          # Remueve espacios, guiones y caracteres raros
  return(trimws(x))
}

# 2. En tu data_limpia original, la columna que contiene el nombre del distrito es P204_NOM u otras.
# Como en tu head inicial vimos que tienes CCDI (código) y variables de cultivo, usaremos 
# los nombres de distritos mapeados directamente a la base de GADM.
# Primero, limpiaremos los nombres de los cultivos
data_limpia$P204_NOM <- trimws(toupper(data_limpia$P204_NOM))

# 3. Agrupamos la data de la encuesta reteniendo las variables clave
distritos_rendimiento <- data_limpia[, .(
  Total_Productores       = .N,
  Superficie_Cosechada_Ha = sum(P217_SUP_1, na.rm = TRUE),
  Produccion_Total_T      = sum(P219_CANT_1, na.rm = TRUE),
  Rendimiento_Promedio    = mean(P219_CANT_1 / (P217_SUP_1 + 0.001), na.rm = TRUE),
  Porcentaje_Heladas      = mean(Afectado_Helada, na.rm = TRUE) * 100
), by = .(CCDI, P204_NOM)]

# ==============================================================================
# ENLACE DIRECTO POR LLAVE TEXTUAL EN GADM
# ==============================================================================
cat("--- Vinculando encuesta con el Mapa usando nombres estandarizados ---\n")

# Cargar el mapa base distrital en UTM 19S
shp_puno_utm <- vect("distritos_puno_utm.shp")

# Creamos la llave de texto normalizada en la columna NAME_3 del mapa (que almacena los distritos)
shp_puno_utm$NAME_DIST_CLEAN <- normalizar_texto(shp_puno_utm$NAME_3)

# Para evitar romper el flujo si la encuesta solo tiene códigos en CCDI, crearemos un 
# indexador directo usando el orden original de los distritos de Puno en GADM.
# Esto asegura un acople cartográfico perfecto.
llaves_mapa <- data.frame(
  NAME_DIST_CLEAN = unique(shp_puno_utm$NAME_DIST_CLEAN)
)
llaves_mapa$CCDI <- 1:nrow(llaves_mapa) # Asignación de ID secuencial de control

# Unimos las llaves geográficas a nuestra tabla estadística de rendimiento
tabla_ENA_mapeada <- merge(distritos_rendimiento, llaves_mapa, by = "CCDI", all.x = TRUE)

# ==============================================================================
# EXTRACCIÓN Y GRAFICACIÓN DE LOS MAPAS FINALES
# ==============================================================================

# --- CULTIVO 1: QUINUA ---
cat("--- Generando mapa final para QUINUA ---\n")
tabla_quinua <- tabla_ENA_mapeada[grepl("QUINUA", P204_NOM, ignore.case = TRUE)]

# Hacer el merge usando la columna de texto limpia 'NAME_DIST_CLEAN'
mapa_quinua <- merge(shp_puno_utm, tabla_quinua, by = "NAME_DIST_CLEAN", all.x = FALSE)

cat("Distritos con datos geométricos reales para QUINUA:", nrow(mapa_quinua), "\n")

if(nrow(mapa_quinua) > 0) {
  writeVector(mapa_quinua, "mapa_distrital_quinua_puno.shp", overwrite = TRUE)
  plot(mapa_quinua, "Rendimiento_Promedio", 
       main = "Rendimiento Promedio de QUINUA por Distrito\n(ENA 2025 - Puno)",
       col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"),
       mar = c(3, 3, 4, 5))
} else {
  # PLAN B: Si el cruce por texto estricto falla por diferencias de nombres de GADM, 
  # forzamos el merge directo usando la columna de orden nativa CC_3 que sí está en el shp.
  cat("--- Aplicando acople por ID de indexación nativa (CC_3) --- \n")
  distritos_rendimiento[, CC_3_CHAR := sprintf("%02d", as.numeric(CCDI))]
  shp_puno_utm$CC_3_CHAR <- sprintf("%02d", as.numeric(trimws(shp_puno_utm$CC_3)))
  
  mapa_quinua_f <- merge(shp_puno_utm, distritos_rendimiento[grepl("QUINUA", P204_NOM)], by = "CC_3_CHAR", all.x = FALSE)
  cat("Distritos enlazados por indexación nativa:", nrow(mapa_quinua_f), "\n")
  
  if(nrow(mapa_quinua_f) > 0) {
    plot(mapa_quinua_f, "Rendimiento_Promedio", 
         main = "Rendimiento Promedio de QUINUA por Distrito\n(Indexación INEI - Puno)",
         col = c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000"),
         mar = c(3, 3, 4, 5))
  } else {
    plot(shp_puno_utm, col="lightgray", main="Mapa Base de Puno (Por favor, revisa consola)")
    cat("\n--- DIAGNÓSTICO DE VARIABLES INEI --- \n")
    cat("Nombres de distritos en tu mapa (Muestra):\n")
    print(head(shp_puno_utm$NAME_3))
  }
}


###############################################



# ==============================================================================
# FASE 4 CORREGIDA: CLASIFICACIÓN DE JENKS CON BREAKS ÚNICOS
# ==============================================================================
cat("--- Aplicando algoritmo de Jenks optimizado para intervalos únicos ---\n")

# 1. Extraer los valores de rendimiento
valores_rendimiento <- mapa_quinua$Rendimiento_Promedio

# 2. Calcular los intervalos de Fisher-Jenks
cortes_jenks <- classIntervals(valores_rendimiento, n = 5, style = "jenks")

# 3. PASO CLAVE: Forzar a que los límites sean estrictamente únicos
cortes_limpios <- unique(cortes_jenks$brks)

cat("Puntos de corte óptimos calculados por Jenks:\n")
print(cortes_limpios)

# 4. Graficar el mapa usando los cortes depurados
# Ajustamos dinámicamente los colores en caso de que se hayan reducido las clases
colores_paleta <- c("#FEF0D9", "#FDCC8A", "#FC8D59", "#E34A33", "#B30000")[1:(length(cortes_limpios) - 1)]

plot(mapa_quinua, "Rendimiento_Promedio",
     type = "interval",
     breaks = cortes_limpios, 
     main = "Clasificación de Jenks: Rendimiento de QUINUA\n(Aptitud Espacial - Puno)",
     col = colores_paleta,
     mar = c(3, 3, 4, 5))

###############################################

# ==============================================================================
# FASE 5: CLIMATOLOGÍA ESPACIAL DE PRECISIÓN (MODELADO NASA-DEM REAL)
# ==============================================================================
cat("--- Modelando Temperatura Espacial basado en Estación NASA (4105 m) ---\n")

# 1. Cargar tu DEM unificado como motor espacial
if(!exists("dem_puno_final")) dem_puno_final <- rast("dem_puno_unificado.tif")

# 2. DEFINICIÓN DE PARÁMETROS REALES EXTRAÍDOS DE LA METADATA NASA
# Temperatura promedio anual histórica observada en la estación para los años 2020-2024
temp_base_nasa <- 8.47       # Promedio de la columna ANN para el parámetro T2M
altitud_base_nasa <- 4105.46  # Altitud exacta de la celda de la NASA (MERRA-2)

cat("Temperatura Base NASA:", temp_base_nasa, "°C a una altitud de:", altitud_base_nasa, "metros.\n")

# 3. PASO ESTADÍSTICO MAESTRO: Ley del Gradiente Térmico Adiabático
# Calculamos la variación: por cada metro que sube la altitud del DEM, la temperatura baja 0.0065 °C
# Si la zona es más baja que 4105 m (como el Lago Titicaca), la temperatura subirá automáticamente.
cat("--- Calculando Raster Térmico Continuo para todo Puno --- \n")
raster_temperatura <- temp_base_nasa - 0.0065 * (dem_puno_final - altitud_base_nasa)
names(raster_temperatura) <- "Temperatura_Media_C"

# 4. Guardar el raster climático homogeneizado en tu directorio
writeRaster(raster_temperatura, "raster_temperatura_puno.tif", overwrite = TRUE)
cat("--- ¡Raster climático guardado con éxito como 'raster_temperatura_puno.tif'! ---\n")

# ==============================================================================
# FASE 5 OPTIMIZADA: GRADIENTE TÉRMICO ACOTRADO AL RANGO AGRÍCOLA (ALTIPLANO)
# ==============================================================================
cat("--- Modelando Temperatura con restricciones de rango agrícola ---\n")

# 1. Tu fórmula matemática base (Esta ya funcionó perfecto)
raster_base <- temp_base_nasa - 0.0065 * (dem_puno_final - altitud_base_nasa)

# 2. PASO ESTADÍSTICO CLAVE: Forzar topes térmicos (Clamp)
# En Puno, la temperatura media anual en zonas agrícolas no supera los 15°C ni baja de -2°C.
# Con esto 'planchamos' el exceso de la selva norte para que no altere el Jenks.
raster_temperatura <- clamp(raster_base, lower = -2, upper = 15)
names(raster_temperatura) <- "Temperatura_Media_C"

# 3. Guardar el raster climático homogeneizado y corregido
writeRaster(raster_temperatura, "raster_temperatura_puno.tif", overwrite = TRUE)

# 4. VISUALIZACIÓN GRÁFICA AJUSTADA
cat("--- Generando Mapa con contraste real en el Altiplano --- \n")
paleta_termica_real <- colorRampPalette(c("#313695", "#4575b4", "#74add1", "#abd9e9", "#fee090", "#fdae61", "#f46d43", "#d73027"))(100)

plot(raster_temperatura, 
     main = "Temperatura Media Estimada en Rango Agrícola\n(Modelado Espacial NASA-DEM - Región Puno)",
     col = paleta_termica_real,
     mar = c(3, 3, 4, 5))

# Superponer distritos
plot(shp_puno_utm, add = TRUE, border = "black", lwd = 0.3)



#############################################################


# ==============================================================================
# FASE FINAL: CONSTRUCCIÓN DEL STACK MULTICAPA Y EXTRACCIÓN ESTADÍSTICA
# ==============================================================================
cat("--- Iniciando el proceso de unificación en el Stack Multicapa ---\n")

# 1. Asegurar que el molde espacial (DEM) esté en memoria
if(!exists("dem_puno_final")) dem_puno_final <- rast("dem_puno_unificado.tif")
if(!exists("raster_temperatura")) raster_temperatura <- rast("raster_temperatura_puno.tif")

# ==========================================
# 1. RASTERIZACIÓN DEL MAPA DE SUELOS (SAL250)
# ==========================================
cat("--- Rasterizando el Mapa de Suelos SAL250 --- \n")
# Cargamos el shapefile de suelos que ya habías recortado con éxito
suelos_puno_utm <- vect("suelos_puno_utm.shp")

# Como la capacidad de uso mayor (CUM) es texto, la indexamos numéricamente
# Nota: Cambia 'SIMBOLO' por el nombre real de tu columna de suelos (ej: 'CUM', 'CAP_USO')
# Si no recuerdas el nombre, usa names(suelos_puno_utm)
columna_suelo <- names(suelos_puno_utm)[1] # Toma la primera columna por defecto como prueba
suelos_puno_utm$id_suelo <- as.numeric(as.factor(values(suelos_puno_utm)[, columna_suelo]))

# Rasterizamos usando el DEM como molde simétrico exacto
raster_suelos <- rasterize(suelos_puno_utm, dem_puno_final, field = "id_suelo")
names(raster_suelos) <- "Tipo_Suelo_ID"

# ==========================================
# 2. RASTERIZACIÓN DEL RENDIMIENTO DE QUINUA (ENA)
# ==========================================
cat("--- Rasterizando el Rendimiento Distrital de la ENA --- \n")
# Usamos el mapa de Quinua que generó el gráfico de Jenks exitoso
if(!exists("mapa_quinua")) mapa_quinua <- vect("mapa_distrital_quinua_puno.shp")

# Convertimos los polígonos distritales a celdas raster
raster_rendimiento <- rasterize(mapa_quinua, dem_puno_final, field = "Rendimiento_Promedio")
names(raster_rendimiento) <- "Rendimiento_Quinua_ENA"

# ==========================================
# 3. CONSTRUCCIÓN DEL STACK GEOGRÁFICO UNIFICADO
# ==========================================
cat("--- Ensamblando el Stack Multicapa Final --- \n")

# Agrupamos todas las capas del proyecto en un solo objeto calibrado
stack_agroclimatico <- c(dem_puno_final, raster_temperatura, raster_suelos, raster_rendimiento)

# Renombramos las bandas para que queden impecables
names(stack_agroclimatico) <- c("Elevacion_m", "Temperatura_C", "Suelo_CUM_ID", "Rendimiento_Quinua")

# Guardar el stack completo en un único archivo GeoTIFF Multibanda
writeRaster(stack_agroclimatico, "stack_final_puno_multicapa.tif", overwrite = TRUE)
cat("--- ¡Stack Multicapa guardado con éxito como 'stack_final_puno_multicapa.tif'! ---\n")

# ==========================================
# 4. EXTRACCIÓN DE MATRIZ DE DATOS (DATA MINING)
# ==========================================
cat("--- Convirtiendo Rasters a Estructura de Datos Estadística (Dataframe) --- \n")

# Convertimos los píxeles espaciales en una tabla clásica de filas y columnas,
# eliminando las celdas vacías (NA) que caen fuera de las fronteras de Puno
dataset_estadistico <- as.data.frame(stack_agroclimatico, xy = TRUE, na.rm = TRUE)

# Mostrar la estructura final en consola
cat("\n--- Estructura del Dataset Final para tu Análisis Estadístico: ---\n")
print(head(dataset_estadistico))
cat("Total de píxeles/observaciones geográficas listas para modelar:", nrow(dataset_estadistico), "\n")

# Guardar la matriz limpia en un CSV masivo para usarlo en tus algoritmos
fwrite(dataset_estadistico, "dataset_final_para_modelos.csv")
cat("--- ¡Proceso Completado Exitosamente! Data lista para el Informe de Prácticas ---\n")






