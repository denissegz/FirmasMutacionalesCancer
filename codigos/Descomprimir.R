##CODIGO NECESARIO PARA DESCOMPRIMIR LAS CARPETAS QUE CONTIENEN LOS ARCHIVOS MAF CUANDO SE DESCARGAN LAS DATOS DE MUTACIONES DE TCGA DATA PORTAL
#El codigo descomprimira automaticamente en "Explorador de archivos" las carpetas, por lo que solo se usa este codigo una vez

#SOLO SE NECESITA CAMBIAR LA LINEA SIGUIENTE
# Define el directorio principal donde esta la carpeta con todos los archivos de las mutaciones
directorio_principal <- "../Downloads/Data/Colorectal"




if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils")
}
library(R.utils)

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)
  
# Lista para almacenar la informacion de los archivos .maf
resultados_lista <- list()

# Obtiene la lista de carpetas dentro del directorio principal
carpetas <- list.dirs(directorio_principal, full.names = TRUE, recursive = TRUE)

# Recorre cada carpeta
for (carpeta in carpetas) {
  # Lista los archivos .gz en la carpeta
  archivos_gz <- list.files(carpeta, pattern = "\\.gz$", full.names = TRUE)
  # Recorre cada archivo .gz
  for (archivo_gz in archivos_gz) {
    # Descomprime el archivo .gz
    archivo_maf <- sub("\\.gz$", "", archivo_gz)
    if (!file.exists(archivo_maf)) {
      gunzip(archivo_gz)  # Descomprime el archivo
    } 
    # Lee el archivo .maf descomprimido
    if (file.exists(archivo_maf)) {
      maf_data <- fread(archivo_maf)  # Usa fread para leer el archivo
      # Agrega la información a la lista
      resultados_lista[[length(resultados_lista) + 1]] <- maf_data
    }
  }
}

# Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
mutaciones <- rbindlist(resultados_lista, fill = TRUE)

# Muestra el data.frame final
View(mutaciones)

