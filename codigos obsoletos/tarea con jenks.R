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
setwd("C:/Tu_Carpeta_de_Trabajo/Proyecto_Puno")

# ==========================================
# ESPACIO GEOGRÁFICO: LÍMITES DE PUNO
# ==========================================
cat("--- Procesando marco geográfico de Puno ---\n")
shp_peru <- st_read("provincias_peru.shp") # Tu shapefile base de provincias

# Filtramos estrictamente el departamento de Puno (Código 21)
shp_puno <- shp_peru[shp_peru$DEPARTAMEN == "PUNO", ]

# Transformamos a coordenadas métricas (UTM Zona 19S / EPSG:32719) para Puno
shp_puno_utm <- st_transform(shp_puno, 32719)

# ==========================================
# ELEVACIÓN: UNIFICAR LAS 15 PIEZAS DEL DEM
# ==========================================
cat("--- Unificando las 15 piezas del DEM (NASA) ---\n")

# Detectar los 15 archivos .zip en tu carpeta de descargas
zips_dem <- list.files("dem_zips/", pattern = "\\.zip$", full.names = TRUE)

# Descomprimir automáticamente cada uno en una carpeta temporal
carpeta_temp <- "dem_descomprimido"
for(zip in zips_dem) {
  unzip(zip, exdir = carpeta_temp)
}

# Leer todos los archivos .hgt o .tif que se extrajeron
archivos_descomprimidos <- list.files(carpeta_temp, pattern = "\\.(hgt|tif)$", full.names = TRUE)
lista_rasters <- lapply(archivos_descomprimidos, rast)

# Fusionar (Mosaico) de los 15 cuadrantes en un único mapa continuo
dem_completo <- do.call(mosaic, lista_rasters)

# Reproyectar a UTM Zona 19S y recortar con la máscara exacta de Puno
dem_utm <- project(dem_completo, "EPSG:32719")
dem_puno_final <- crop(dem_utm, shp_puno_utm) |> mask(shp_puno_utm)
names(dem_puno_final) <- "Elevacion_Metros"

# Guardar tu DEM unificado para no tener que repetir este proceso pesado
writeRaster(dem_puno_final, "dem_puno_unificado.tif", overwrite = TRUE)

# ==========================================
# EDALFOLOGÍA: RECORTAR SUELOS (SAL250)
# ==========================================
cat("--- Recortando Mapa de Suelos SAL250 ---\n")
suelos_peru <- st_read("suelos_peru.shp") # Tu mapa del MINAM de todo el Perú

# Intersección espacial para extraer solo los suelos que caen dentro de Puno
suelos_puno <- st_intersection(suelos_peru, st_buffer(shp_puno, 0))
suelos_puno_utm <- st_transform(suelos_puno, 32719)

# Guardar el shapefile de suelos exclusivo de Puno
st_write(suelos_puno_utm, "suelos_puno_utm.shp", delete_layer = TRUE)

# ==========================================
# VERIFICACIÓN DE DATA DE LA ENCUESTA (SPSS/CSV)
# ==========================================
cat("--- Cargando base de datos agropecuaria filtrada ---\n")
# Usamos data.table para una lectura instantánea
data_agro <- fread("Data_Agro_puno.csv") 

# Mostrar las primeras filas y estructura en consola para asegurar las variables
print(head(data_agro))