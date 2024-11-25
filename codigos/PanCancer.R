##### FIRMAS MUTACIONALES DE CANCER CONTRA OTROS CANCERES
#En cada cancer se tiene que modificar la ruta de los datos clinicos y de mutaciones, este codigo lee archivos de mutaciones de TCGA en formato MAF
#En la parte de "Clasificación binaria" de la seccion "CARGAR DATOS" se cambia al numero del cancer de interes, para que este sea la clase 1 y todos los demas canceres sean de la clase 0
#En la seccion de "VALIDACION" se tiene que cambiar las rutas de los archivos, ese codigo lee mutaciones en archivo CSV (del portal CBioPortal)

##LIBRERIAS
{
  #Si no se tienen, intalar librerias
  libraries <- c("reshape2", "ggplot2", "dplyr", "caret", "tidyr", "rpart", "RColorBrewer",
                 "rpart.plot", "randomForest", "class", "FactoMineR",
                 "caretEnsemble", "car", "pROC", "data.table", 
                 "isotree", "dbscan", "stringr", "factoextra")
  
  # Recorremos la lista de librerías
  for (lib in libraries) {
    if (!requireNamespace(lib, quietly = TRUE)) {
      install.packages(lib)
      library(lib, character.only = TRUE)
      cat("Librería", lib, "instalada y cargada.\n")
    } else {
      cat("Librería", lib, "ya está instalada.\n")
    }
  }
  
  library(reshape2)
  library(ggplot2)
  library(dplyr)
  library(caret)
  library(tidyr)
  library(rpart)
  library(rpart.plot)
  library(randomForest)
  library(caret)       
  library(class)
  library(caretEnsemble)
  library(car)
  library(pROC)
  library(data.table)
  library(isotree)
  library(dbscan)
  library(stringr)
  library(FactoMineR)
  library(factoextra)
  library(RColorBrewer)
}

##CARGAR DATOS
{
  #1.Colorectal
  {
    #Clínico
    {
    clinical1 <- read.csv("../Data/Colorectal/clinical.csv", header = TRUE, sep = ",")
    #View(clinical1)
    # Reemplazar ''--' por NA
    clinical1[clinical1 == "'--"] <- NA
    #Contar estapas clínica y paológica
    count_clinical1 <- sum(!is.na(clinical1$ajcc_clinical_stage))
    print(count_clinical1)
    count_path1 <- sum(!is.na(clinical1$ajcc_pathologic_stage))
    print(count_path1)
    #View(clinical1)
    # Solo filas con stage definida
    clinical_clean1 <- clinical1
      #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
    #Diferencias de etapa clínica y patológica
    diferencias1 <- clinical_clean1 %>%
      filter(ajcc_pathologic_stage != ajcc_clinical_stage) %>%
      select(ajcc_pathologic_stage, ajcc_clinical_stage)
    print(diferencias1)
    #Quedarse con las columnas importantes
    clinical_clean1 <- clinical_clean1 %>%
      select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
    #View(clinical_clean1)
    #Num de pacientes
    num_unicos1 <- clinical_clean1 %>%
      summarise(count = n_distinct(case_id))
    print(num_unicos1)
    # Eliminar duplicados en clinical_clean1
    clinical_clean1_unique <- clinical_clean1 %>%
      distinct(case_id, .keep_all = TRUE)
    #Eliminar metastasis
    clinical_final1 <- clinical_clean1_unique[!grepl("metastasis", clinical_clean1_unique$classification_of_tumor, ignore.case = TRUE), ]
    #View(clinical_final1)
    #num pacinetes
    num_unicos1 <- clinical_final1 %>%
      summarise(count = n_distinct(case_id))
    print(num_unicos1)
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal1 <- "../Data/Colorectal"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista1 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas1 <- list.dirs(directorio_principal1, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta1 in carpetas1) {
        # Lista los archivos .maf en la carpeta
        archivos_maf1 <- list.files(carpeta1, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf1 in archivos_maf1) {
          # Lee el archivo .maf
          maf_data1 <- fread(archivo_maf1)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista1[[length(resultados_lista1) + 1]] <- maf_data1
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones1 <- rbindlist(resultados_lista1, fill = TRUE)
      # Dataframe final
      #View(mutaciones1)
      #Genes no especificados
      diferencias1 <- mutaciones1 %>%
        filter(SYMBOL != Hugo_Symbol) %>%
        select(SYMBOL, Hugo_Symbol)
      print(diferencias1)
      mutaciones_clean1 <- mutaciones1 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somaticas
      mutaciones_unicas1 <- unique(mutaciones1$Mutation_Status)
      print(mutaciones_unicas1)
      mutaciones_clean1 <- mutaciones_clean1 %>%
        filter(!is.na(Mutation_Status))
      #Añadir columna con identificador
      mutaciones_clean1 <- mutaciones_clean1 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      #
      diferencias1 <- mutaciones_clean1 %>%
        filter(Tumor_Seq_Allele1 != Reference_Allele) %>%
        select(Tumor_Seq_Allele1, Reference_Allele)
      print(diferencias1)
      #Quedarse con las columnas importantes
      mutaciones_clean1 <- mutaciones_clean1 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT )
      #View(mutaciones_clean1)
      #Num de pacientes
      num_unicos1 <- mutaciones_clean1 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos1)
    }
    #Unir
    {
      resultado1 <- merge(mutaciones_clean1, clinical_final1, by = "case_id", all=FALSE)
      #View(resultado)
      #Eliminar duplicados
      data1 <- resultado1 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      #Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data1, columns) {
        for (col in columns) {
          data1 <- data1 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data1)
      }
      # Aplicar la función a los datos
      data1 <- clean_columns(data1, columns_to_clean)
      #View(data1)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data1 <- data1 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data1 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data1 <- data1 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data1)
      #Num de pacientes
      num_unicos1 <- data1 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos1)
      
      #Codificación de variables categóricas
      data_final1 <- data1%>%
        mutate(classif_cancer = 1)%>%
        mutate(tipo_cancer = "Colorectal")
      View(data_final1)
    }
  }
  #2.Pancreas
  {
   #Clínicos
    {
      clinical2 <- read.csv("../Data/Pancreas/clinical.csv", header = TRUE, sep = ",")
      #View(clinical2)
      # Reemplazar ''--' por NA
      clinical2[clinical2 == "'--"] <- NA
      # Contar etapas clínica y patológica
      count_clinical2 <- sum(!is.na(clinical2$ajcc_clinical_stage))
      print(count_clinical2)
      count_path2 <- sum(!is.na(clinical2$ajcc_pathologic_stage))
      print(count_path2)
      # View(clinical2)
      # Solo filas con stage definida
      clinical_clean2 <- clinical2
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Diferencias de etapa clínica y patológica
      diferencias2 <- clinical_clean2 %>%
        filter(ajcc_pathologic_stage != ajcc_clinical_stage) %>%
        select(ajcc_pathologic_stage, ajcc_clinical_stage)
      print(diferencias2)
      # Quedarse con las columnas importantes
      clinical_clean2 <- clinical_clean2 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      #View(clinical_clean2)
      # Num de pacientes
      num_unicos2 <- clinical_clean2 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos2)
      # Eliminar duplicados en clinical_clean2
      clinical_clean2_unique <- clinical_clean2 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final2 <- clinical_clean2_unique[!grepl("metastasis", clinical_clean2_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final2)
      # Num pacientes
      num_unicos2 <- clinical_final2 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos2)

    }
    #Mutaciones
    {
      library(R.utils)
      library(data.table)
      # Define el directorio principal donde están las carpetas
      directorio_principal2 <- "../Data/Pancreas"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista2 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas2 <- list.dirs(directorio_principal2, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta2 in carpetas2) {
        # Lista los archivos .maf en la carpeta
        archivos_maf2 <- list.files(carpeta2, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf2 in archivos_maf2) {
          # Lee el archivo .maf
          maf_data2 <- fread(archivo_maf2)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista2[[length(resultados_lista2) + 1]] <- maf_data2
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones2 <- rbindlist(resultados_lista2, fill = TRUE)
      # Dataframe final
      # View(mutaciones2)
      # Genes no especificados
      mutaciones_clean2 <- mutaciones2 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas2 <- unique(mutaciones2$Mutation_Status)
      print(mutaciones_unicas2)
      mutaciones_clean2 <- mutaciones_clean2 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean2 <- mutaciones_clean2 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean2 <- mutaciones_clean2 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean2)
      # Num de pacientes
      num_unicos2 <- mutaciones_clean2 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos2)
      
    }
    #Unir
    {
      resultado2 <- merge(mutaciones_clean2, clinical_final2, by = "case_id", all = FALSE)
      #View(resultado)
      #Eliminar duplicados
      data2 <- resultado2 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      #Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data2, columns) {
        for (col in columns) {
          data2 <- data2 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data2)
      }
      # Aplicar la función a los datos
      data2 <- clean_columns(data2, columns_to_clean)
      #View(data2)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data2 <- data2 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data2 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data2 <- data2 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data2)
      #Num de pacientes
      num_unicos2 <- data2 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos2)
      
      
      #Codificación de variables categóricas
      data_final2 <- data2%>%
        mutate(classif_cancer = 2)%>%
        mutate(tipo_cancer = "Pancreas")
      View(data_final2)
      
    }
  }
  #3.Liver
  {
    #Clinical
    {
      clinical3 <- read.csv("../Data/Liver/clinical.csv", header = TRUE, sep = ",")
      # View(clinical3)
      # Reemplazar ''--' por NA
      clinical3[clinical3 == "'--"] <- NA
      # Contar etapas clínica y patológica
      count_clinical3 <- sum(!is.na(clinical3$ajcc_clinical_stage))
      print(count_clinical3)
      count_path3 <- sum(!is.na(clinical3$ajcc_pathologic_stage))
      print(count_path3)
      # View(clinical3)
      # Solo filas con stage definida
      clinical_clean3 <- clinical3
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean3 <- clinical_clean3 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean3)
      # Num de pacientes
      num_unicos3 <- clinical_clean3 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos3)
      # Eliminar duplicados en clinical_clean3
      clinical_clean3_unique <- clinical_clean3 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final3 <- clinical_clean3_unique[!grepl("metastasis", clinical_clean3_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final3)
      # Num pacientes
      num_unicos3 <- clinical_final3 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos3)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal3 <- "../Data/Liver"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista3 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas3 <- list.dirs(directorio_principal3, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta3 in carpetas3) {
        # Lista los archivos .maf en la carpeta
        archivos_maf3 <- list.files(carpeta3, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf3 in archivos_maf3) {
          # Lee el archivo .maf
          maf_data3 <- fread(archivo_maf3)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista3[[length(resultados_lista3) + 1]] <- maf_data3
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones3 <- rbindlist(resultados_lista3, fill = TRUE)
      # Dataframe final
      # View(mutaciones3)
      # Genes no especificados
      mutaciones_clean3 <- mutaciones3 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas3 <- unique(mutaciones3$Mutation_Status)
      print(mutaciones_unicas3)
      mutaciones_clean3 <- mutaciones_clean3 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean3 <- mutaciones_clean3 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean3 <- mutaciones_clean3 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean3)
      # Num de pacientes
      num_unicos3 <- mutaciones_clean3 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos3)
      
    }
    #Unir
    {
      resultado3 <- merge(mutaciones_clean3, clinical_final3, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data3 <- resultado3 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data3, columns) {
        for (col in columns) {
          data3 <- data3 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data3)
      }
      # Aplicar la función a los datos
      data3 <- clean_columns(data3, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data3 <- data3 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data3 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data3 <- data3 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data3)
      # Num de pacientes
      num_unicos3 <- data3 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos3)
      
      # Codificación de variables categóricas
      data_final3 <- data3%>%
        mutate(classif_cancer = 3)%>%
        mutate(tipo_cancer = "Liver")
      View(data_final3)
      
    }
  }
  #4.Stomach
  {
    #Clinico
    {
      clinical4 <- read.csv("../Data/Stomach/clinical.csv", header = TRUE, sep = ",")
      # View(clinical4)
      # Reemplazar ''--' por NA
      clinical4[clinical4 == "'--"] <- NA
      # View(clinical4)
      # Solo filas con stage definida
      clinical_clean4 <- clinical4 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean4 <- clinical_clean4 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean4)
      # Num de pacientes
      num_unicos4 <- clinical_clean4 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos4)
      # Eliminar duplicados en clinical_clean4
      clinical_clean4_unique <- clinical_clean4 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final4 <- clinical_clean4_unique[!grepl("metastasis", clinical_clean4_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final4)
      # Num pacientes
      num_unicos4 <- clinical_final4 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos4)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal4 <- "../Data/Stomach"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista4 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas4 <- list.dirs(directorio_principal4, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta4 in carpetas4) {
        # Lista los archivos .maf en la carpeta
        archivos_maf4 <- list.files(carpeta4, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf4 in archivos_maf4) {
          # Lee el archivo .maf
          maf_data4 <- fread(archivo_maf4)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista4[[length(resultados_lista4) + 1]] <- maf_data4
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones4 <- rbindlist(resultados_lista4, fill = TRUE)
      # Dataframe final
      # View(mutaciones4)
      # Genes no especificados
      mutaciones_clean4 <- mutaciones4 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas4 <- unique(mutaciones4$Mutation_Status)
      print(mutaciones_unicas4)
      mutaciones_clean4 <- mutaciones_clean4 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean4 <- mutaciones_clean4 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean4 <- mutaciones_clean4 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean4)
      # Num de pacientes
      num_unicos4 <- mutaciones_clean4 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos4)
      
      
    }
    #Unir
    {
      resultado4 <- merge(mutaciones_clean4, clinical_final4, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data4 <- resultado4 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data4, columns) {
        for (col in columns) {
          data4 <- data4 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data4)
      }
      # Aplicar la función a los datos
      data4 <- clean_columns(data4, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data4 <- data4 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data4 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data4 <- data4 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data4)
      # Num de pacientes
      num_unicos4 <- data4 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos4)
      
      # Codificación de variables categóricas
      data_final4 <- data4%>%
        mutate(classif_cancer = 4)%>%
        mutate(tipo_cancer = "Stomach")
      View(data_final4)
      
    }
  }
  #5. Brain
  {
    #Clinico
    {
      clinical5 <- read.csv("../Data/Brain/clinical.csv", header = TRUE, sep = ",")
      # View(clinical5)
      # Reemplazar ''--' por NA
      clinical5[clinical5 == "'--"] <- NA
      # Contar etapas clínica y patológica
      count_clinical5 <- sum(!is.na(clinical5$ajcc_clinical_stage))
      print(count_clinical5)
      count_path5 <- sum(!is.na(clinical5$ajcc_pathologic_stage))
      print(count_path5)
      # View(clinical5)
      # Solo filas con stage definida
      clinical_clean5 <- clinical5
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean5 <- clinical_clean5 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean5)
      # Num de pacientes
      num_unicos5 <- clinical_clean5 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos5)
      # Eliminar duplicados en clinical_clean5
      clinical_clean5_unique <- clinical_clean5 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final5 <- clinical_clean5_unique[!grepl("metastasis", clinical_clean5_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final5)
      # Num pacientes
      num_unicos5 <- clinical_final5 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos5)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal5 <- "../Data/Brain"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista5 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas5 <- list.dirs(directorio_principal5, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta5 in carpetas5) {
        # Lista los archivos .maf en la carpeta
        archivos_maf5 <- list.files(carpeta5, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf5 in archivos_maf5) {
          # Lee el archivo .maf
          maf_data5 <- fread(archivo_maf5)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista5[[length(resultados_lista5) + 1]] <- maf_data5
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones5 <- rbindlist(resultados_lista5, fill = TRUE)
      # Dataframe final
      # View(mutaciones5)
      # Genes no especificados
      mutaciones_clean5 <- mutaciones5 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas5 <- unique(mutaciones5$Mutation_Status)
      print(mutaciones_unicas5)
      mutaciones_clean5 <- mutaciones_clean5 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean5 <- mutaciones_clean5 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean5 <- mutaciones_clean5 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean5)
      # Num de pacientes
      num_unicos5 <- mutaciones_clean5 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos5)
      
      
    }
    #Unir
    {
      resultado5 <- merge(mutaciones_clean5, clinical_final5, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data5 <- resultado5 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data5, columns) {
        for (col in columns) {
          data5 <- data5 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data5)
      }
      # Aplicar la función a los datos
      data5 <- clean_columns(data5, columns_to_clean)
      # View(data)
      #num_unicos5 <- data5 %>%
       # summarise(count = n_distinct(Identificador))
      #print(num_unicos5)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data5 <- data5 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data5 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data5 <- data5 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data5)
      # Num de pacientes
      num_unicos5 <- data5 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos5)
      
      # Codificación de variables categóricas
      data_final5 <- data5%>%
        mutate(classif_cancer = 5)%>%
        mutate(tipo_cancer = "Brain")
      View(data_final5)
      
    }
  }
  #6.Esophagus
  {
    #Clinico
    {
      clinical6 <- read.csv("../Data/Esophagus/clinical.csv", header = TRUE, sep = ",")
      # View(clinical6)
      # Reemplazar ''--' por NA
      clinical6[clinical6 == "'--"] <- NA
      # Contar etapas clínica y patológica
      count_clinical6 <- sum(!is.na(clinical6$ajcc_clinical_stage))
      print(count_clinical6)
      count_path6 <- sum(!is.na(clinical6$ajcc_pathologic_stage))
      print(count_path6)
      # View(clinical6)
      # Solo filas con stage definida
      clinical_clean6 <- clinical6 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean6 <- clinical_clean6 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean6)
      # Num de pacientes
      num_unicos6 <- clinical_clean6 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos6)
      # Eliminar duplicados en clinical_clean6
      clinical_clean6_unique <- clinical_clean6 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final6 <- clinical_clean6_unique[!grepl("metastasis", clinical_clean6_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final6)
      # Num pacientes
      num_unicos6 <- clinical_final6 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos6)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal6 <- "../Data/Esophagus"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista6 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas6 <- list.dirs(directorio_principal6, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta6 in carpetas6) {
        # Lista los archivos .maf en la carpeta
        archivos_maf6 <- list.files(carpeta6, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf6 in archivos_maf6) {
          # Lee el archivo .maf
          maf_data6 <- fread(archivo_maf6)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista6[[length(resultados_lista6) + 1]] <- maf_data6
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones6 <- rbindlist(resultados_lista6, fill = TRUE)
      # Dataframe final
      # View(mutaciones6)
      # Genes no especificados
      mutaciones_clean6 <- mutaciones6 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas6 <- unique(mutaciones6$Mutation_Status)
      print(mutaciones_unicas6)
      mutaciones_clean6 <- mutaciones_clean6 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean6 <- mutaciones_clean6 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean6 <- mutaciones_clean6 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean6)
      # Num de pacientes
      num_unicos6 <- mutaciones_clean6 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos6)
      
      
    }
    #Unir
    {
      resultado6 <- merge(mutaciones_clean6, clinical_final6, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data6 <- resultado6 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data6, columns) {
        for (col in columns) {
          data6 <- data6 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data6)
      }
      # Aplicar la función a los datos
      data6 <- clean_columns(data6, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data6 <- data6 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data6 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data6 <- data6 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data6)
      # Num de pacientes
      num_unicos6 <- data6 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos6)
      
      # Codificación de variables categóricas
      data_final6 <- data6%>%
        mutate(classif_cancer = 6)%>%
        mutate(tipo_cancer = "Esophagus")
      View(data_final6)
      
    }
  }
  #7.Bladder
  {
    #Clinico
    {
      clinical7 <- read.csv("../Data/Bladder/clinical.csv", header = TRUE, sep = ",")
    # View(clinical7)
    # Reemplazar ''--' por NA
    clinical7[clinical7 == "'--"] <- NA
    # Contar etapas clínica y patológica
    count_clinical7 <- sum(!is.na(clinical7$ajcc_clinical_stage))
    print(count_clinical7)
    count_path7 <- sum(!is.na(clinical7$ajcc_pathologic_stage))
    print(count_path7)
    # View(clinical7)
    # Solo filas con stage definida
    clinical_clean7 <- clinical7 
      #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
    # Quedarse con las columnas importantes
    clinical_clean7 <- clinical_clean7 %>%
      select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
    # View(clinical_clean7)
    # Num de pacientes
    num_unicos7 <- clinical_clean7 %>%
      summarise(count = n_distinct(case_id))
    print(num_unicos7)
    # Eliminar duplicados en clinical_clean7
    clinical_clean7_unique <- clinical_clean7 %>%
      distinct(case_id, .keep_all = TRUE)
    # Eliminar metastasis
    clinical_final7 <- clinical_clean7_unique[!grepl("metastasis", clinical_clean7_unique$classification_of_tumor, ignore.case = TRUE), ]
    #View(clinical_final7)
    # Num pacientes
    num_unicos7 <- clinical_final7 %>%
      summarise(count = n_distinct(case_id))
    print(num_unicos7)
    
  }
  #Mutaciones
  {
    # Define el directorio principal donde están las carpetas
    directorio_principal7 <- "../Data/Bladder"
    # Inicializa una lista para almacenar la información de los archivos .maf
    resultados_lista7 <- list()
    # Obtiene la lista de carpetas dentro del directorio principal
    carpetas7 <- list.dirs(directorio_principal7, full.names = TRUE, recursive = TRUE)
    # Recorre cada carpeta
    for (carpeta7 in carpetas7) {
      # Lista los archivos .maf en la carpeta
      archivos_maf7 <- list.files(carpeta7, pattern = "\\.maf$", full.names = TRUE)
      # Recorre cada archivo .maf
      for (archivo_maf7 in archivos_maf7) {
        # Lee el archivo .maf
        maf_data7 <- fread(archivo_maf7)  # fread para leer el archivo
        # Agrega la información a la lista
        resultados_lista7[[length(resultados_lista7) + 1]] <- maf_data7
      }
    }
    # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
    mutaciones7 <- rbindlist(resultados_lista7, fill = TRUE)
    # Dataframe final
    # View(mutaciones7)
    # Genes no especificados
    mutaciones_clean7 <- mutaciones7 %>%
      filter(Hugo_Symbol != "Unknown")
    # Solo filas con mutaciones somáticas
    mutaciones_unicas7 <- unique(mutaciones7$Mutation_Status)
    print(mutaciones_unicas7)
    mutaciones_clean7 <- mutaciones_clean7 %>%
      filter(!is.na(Mutation_Status))
    # Añadir columna con identificador
    mutaciones_clean7 <- mutaciones_clean7 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
    # Quedarse con las columnas importantes
    mutaciones_clean7 <- mutaciones_clean7 %>%
      select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
    #View(mutaciones_clean7)
    # Num de pacientes
    num_unicos7 <- mutaciones_clean7 %>%
      summarise(count = n_distinct(case_id))
    print(num_unicos7)
  }
  #Unir
    {
      resultado7 <- merge(mutaciones_clean7, clinical_final7, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data7 <- resultado7 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data7, columns) {
        for (col in columns) {
          data7 <- data7 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data7)
      }
      # Aplicar la función a los datos
      data7 <- clean_columns(data7, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data7 <- data7 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data7 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data7 <- data7 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data7)
      # Num de pacientes
      num_unicos7 <- data7 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos7)
      
      # Codificación de variables categóricas
      data_final7 <- data7 %>%
        mutate(classif_cancer = 7) %>%
        mutate(tipo_cancer = "Bladder")
      View(data_final7)
    }
  }
  #8.Lung
  {
    #Clinico
    {
      clinical8 <- read.csv("../Data/Lung/clinical.csv", header = TRUE, sep = ",")
      # View(clinical8)
      # Reemplazar ''--' por NA
      clinical8[clinical8 == "'--"] <- NA
      # Contar etapas clínica y patológica
      count_clinical8 <- sum(!is.na(clinical8$ajcc_clinical_stage))
      print(count_clinical8)
      count_path8 <- sum(!is.na(clinical8$ajcc_pathologic_stage))
      print(count_path8)
      # View(clinical8)
      # Solo filas con stage definida
      clinical_clean8 <- clinical8 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean8 <- clinical_clean8 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean8)
      # Num de pacientes
      num_unicos8 <- clinical_clean8 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos8)
      # Eliminar duplicados en clinical_clean8
      clinical_clean8_unique <- clinical_clean8 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final8 <- clinical_clean8_unique[!grepl("metastasis", clinical_clean8_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final8)
      # Num pacientes
      num_unicos8 <- clinical_final8 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos8)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal8 <- "../Data/Lung"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista8 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas8 <- list.dirs(directorio_principal8, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta8 in carpetas8) {
        # Lista los archivos .maf en la carpeta
        archivos_maf8 <- list.files(carpeta8, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf8 in archivos_maf8) {
          # Lee el archivo .maf
          maf_data8 <- fread(archivo_maf8)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista8[[length(resultados_lista8) + 1]] <- maf_data8
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones8 <- rbindlist(resultados_lista8, fill = TRUE)
      # Dataframe final
      # View(mutaciones8)
      # Genes no especificados
      mutaciones_clean8 <- mutaciones8 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas8 <- unique(mutaciones8$Mutation_Status)
      print(mutaciones_unicas8)
      mutaciones_clean8 <- mutaciones_clean8 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean8 <- mutaciones_clean8 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean8 <- mutaciones_clean8 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean8)
      # Num de pacientes
      num_unicos8 <- mutaciones_clean8 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos8)
    }
    #Unir
    {
      resultado8 <- merge(mutaciones_clean8, clinical_final8, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data8 <- resultado8 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data8, columns) {
        for (col in columns) {
          data8 <- data8 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data8)
      }
      # Aplicar la función a los datos
      data8 <- clean_columns(data8, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data8 <- data8 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data8 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data8 <- data8 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data8)
      # Num de pacientes
      num_unicos8 <- data8 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos8)
      
      # Codificación de variables categóricas
      data_final8 <- data8 %>%
        mutate(classif_cancer = 8) %>%
        mutate(tipo_cancer = "Lung")
      View(data_final8)
    }
    
  }
  #9.Breast
  {
    #Clinico
    {
      clinical9 <- read.csv("../Data/Breast/clinical.csv", header = TRUE, sep = ",")
      # View(clinical9)
      # Reemplazar ''--' por NA
      clinical9[clinical9 == "'--"] <- NA
      # Contar etapas clínica y patológica
      count_clinical9 <- sum(!is.na(clinical9$ajcc_clinical_stage))
      print(count_clinical9)
      count_path9 <- sum(!is.na(clinical9$ajcc_pathologic_stage))
      print(count_path9)
      # View(clinical9)
      # Solo filas con stage definida
      clinical_clean9 <- clinical9 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean9 <- clinical_clean9 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean9)
      # Num de pacientes
      num_unicos9 <- clinical_clean9 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos9)
      # Eliminar duplicados en clinical_clean9
      clinical_clean9_unique <- clinical_clean9 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final9 <- clinical_clean9_unique[!grepl("metastasis", clinical_clean9_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final9)
      # Num pacientes
      num_unicos9 <- clinical_final9 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos9)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal9 <- "../Data/Breast"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista9 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas9 <- list.dirs(directorio_principal9, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta9 in carpetas9) {
        # Lista los archivos .maf en la carpeta
        archivos_maf9 <- list.files(carpeta9, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf9 in archivos_maf9) {
          # Lee el archivo .maf
          maf_data9 <- fread(archivo_maf9)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista9[[length(resultados_lista9) + 1]] <- maf_data9
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones9 <- rbindlist(resultados_lista9, fill = TRUE)
      # Dataframe final
      # View(mutaciones9)
      # Genes no especificados
      mutaciones_clean9 <- mutaciones9 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas9 <- unique(mutaciones9$Mutation_Status)
      print(mutaciones_unicas9)
      mutaciones_clean9 <- mutaciones_clean9 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean9 <- mutaciones_clean9 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean9 <- mutaciones_clean9 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean9)
      # Num de pacientes
      num_unicos9 <- mutaciones_clean9 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos9)
    }
    #Unir
    {
      resultado9 <- merge(mutaciones_clean9, clinical_final9, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data9 <- resultado9 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data9, columns) {
        for (col in columns) {
          data9 <- data9 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data9)
      }
      # Aplicar la función a los datos
      data9 <- clean_columns(data9, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data9 <- data9 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data9 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data9 <- data9 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data9)
      # Num de pacientes
      num_unicos9 <- data9 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos9)
  
      # Codificación de variables categóricas
      data_final9 <- data9 %>%
        mutate(classif_cancer = 9) %>%
        mutate(tipo_cancer = "Breast")
      View(data_final9)
    }
    
  }
  #10. Kidney
  {
    #Clinico
    {
      clinical10 <- read.csv("../Data/Kidney/clinical.csv", header = TRUE, sep = ",")
      # View(clinical10)
      # Reemplazar ''--' por NA
      clinical10[clinical10 == "'--"] <- NA
      # View(clinical10)
      # Solo filas con stage definida
      clinical_clean10 <- clinical10 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean10 <- clinical_clean10 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean10)
      # Num de pacientes
      num_unicos10 <- clinical_clean10 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos10)
      # Eliminar duplicados en clinical_clean10
      clinical_clean10_unique <- clinical_clean10 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final10 <- clinical_clean10_unique[!grepl("metastasis", clinical_clean10_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final10)
      # Num pacientes
      num_unicos10 <- clinical_final10 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos10)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal10 <- "../Data/Kidney"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista10 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas10 <- list.dirs(directorio_principal10, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta10 in carpetas10) {
        # Lista los archivos .maf en la carpeta
        archivos_maf10 <- list.files(carpeta10, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf10 in archivos_maf10) {
          # Lee el archivo .maf
          maf_data10 <- fread(archivo_maf10)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista10[[length(resultados_lista10) + 1]] <- maf_data10
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones10 <- rbindlist(resultados_lista10, fill = TRUE)
      # Dataframe final
      # View(mutaciones10)
      # Genes no especificados
      mutaciones_clean10 <- mutaciones10 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas10 <- unique(mutaciones10$Mutation_Status)
      print(mutaciones_unicas10)
      mutaciones_clean10 <- mutaciones_clean10 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean10 <- mutaciones_clean10 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean10 <- mutaciones_clean10 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean10)
      # Num de pacientes
      num_unicos10 <- mutaciones_clean10 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos10)
    }
    #Unir
    {
      resultado10 <- merge(mutaciones_clean10, clinical_final10, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data10 <- resultado10 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data10, columns) {
        for (col in columns) {
          data10 <- data10 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data10)
      }
      # Aplicar la función a los datos
      data10 <- clean_columns(data10, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data10 <- data10 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data10 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data10 <- data10 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data10)
      # Num de pacientes
      num_unicos10 <- data10 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos10)
      
      # Codificación de variables categóricas
      data_final10 <- data10 %>%
        mutate(classif_cancer = 10) %>%
        mutate(tipo_cancer = "Kidney")
      View(data_final10)
    }
    
  }
  #11. Adrenal gland
  {
    #Clinico
    {
      clinical11 <- read.csv("../Data/Adrenal/clinical.csv", header = TRUE, sep = ",")
      # View(clinical11)
      # Reemplazar ''--' por NA
      clinical11[clinical11 == "'--"] <- NA
      # View(clinical11)
      # Solo filas con stage definida
      clinical_clean11 <- clinical11 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean11 <- clinical_clean11 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean11)
      # Num de pacientes
      num_unicos11 <- clinical_clean11 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos11)
      # Eliminar duplicados en clinical_clean11
      clinical_clean11_unique <- clinical_clean11 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final11 <- clinical_clean11_unique[!grepl("metastasis", clinical_clean11_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final11)
      # Num pacientes
      num_unicos11 <- clinical_final11 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos11)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal11 <- "../Data/Adrenal"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista11 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas11 <- list.dirs(directorio_principal11, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta11 in carpetas11) {
        # Lista los archivos .maf en la carpeta
        archivos_maf11 <- list.files(carpeta11, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf11 in archivos_maf11) {
          # Lee el archivo .maf
          maf_data11 <- fread(archivo_maf11)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista11[[length(resultados_lista11) + 1]] <- maf_data11
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones11 <- rbindlist(resultados_lista11, fill = TRUE)
      # Dataframe final
       #View(mutaciones11)
      # Genes no especificados
      mutaciones_clean11 <- mutaciones11 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas11 <- unique(mutaciones11$Mutation_Status)
      print(mutaciones_unicas11)
      mutaciones_clean11 <- mutaciones_clean11 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean11 <- mutaciones_clean11 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean11 <- mutaciones_clean11 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean11)
      # Num de pacientes
      num_unicos11 <- mutaciones_clean11 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos11)
    }
    #Unir
    {
      resultado11 <- merge(mutaciones_clean11, clinical_final11, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data11 <- resultado11 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data11, columns) {
        for (col in columns) {
          data11 <- data11 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data11)
      }
      # Aplicar la función a los datos
      data11 <- clean_columns(data11, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data11 <- data11 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data11 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data11 <- data11 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      #View(data11)
      # Num de pacientes
      num_unicos11 <- data11 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos11)
      
      # Codificación de variables categóricas
      data_final11 <- data11 %>%
        mutate(classif_cancer = 11) %>%
        mutate(tipo_cancer = "Adrenal")
      View(data_final11)
    }
    
  }
  #12. Head and neck
  {
    #Clinico
    {
      clinical12 <- read.csv("../Data/Head/clinical.csv", header = TRUE, sep = ",")
      # View(clinical12)
      # Reemplazar ''--' por NA
      clinical12[clinical12 == "'--"] <- NA
      # View(clinical12)
      # Solo filas con stage definida
      clinical_clean12 <- clinical12 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean12 <- clinical_clean12 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean12)
      # Num de pacientes
      num_unicos12 <- clinical_clean12 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos12)
      # Eliminar duplicados en clinical_clean12
      clinical_clean12_unique <- clinical_clean12 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final12 <- clinical_clean12_unique[!grepl("metastasis", clinical_clean12_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final12)
      # Num pacientes
      num_unicos12 <- clinical_final12 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos12)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal12 <- "../Data/Head"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista12 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas12 <- list.dirs(directorio_principal12, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta12 in carpetas12) {
        # Lista los archivos .maf en la carpeta
        archivos_maf12 <- list.files(carpeta12, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf12 in archivos_maf12) {
          # Lee el archivo .maf
          maf_data12 <- fread(archivo_maf12)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista12[[length(resultados_lista12) + 1]] <- maf_data12
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones12 <- rbindlist(resultados_lista12, fill = TRUE)
      # Dataframe final
      # View(mutaciones12)
      # Genes no especificados
      mutaciones_clean12 <- mutaciones12 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas12 <- unique(mutaciones12$Mutation_Status)
      print(mutaciones_unicas12)
      mutaciones_clean12 <- mutaciones_clean12 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean12 <- mutaciones_clean12 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean12 <- mutaciones_clean12 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean12)
      # Num de pacientes
      num_unicos12 <- mutaciones_clean12 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos12)
    }
    #Unir
    {
      resultado12 <- merge(mutaciones_clean12, clinical_final12, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data12 <- resultado12 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data12, columns) {
        for (col in columns) {
          data12 <- data12 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data12)
      }
      # Aplicar la función a los datos
      data12 <- clean_columns(data12, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data12 <- data12 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data12 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data12 <- data12 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data12)
      # Num de pacientes
      num_unicos12 <- data12 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos12)
      
      # Codificación de variables categóricas
      data_final12 <- data12 %>%
        mutate(classif_cancer = 12) %>%
        mutate(tipo_cancer = "Head and Neck")
      View(data_final12)
    }
    
  }
  #13. Uterus
  {
    #Clinico
    {
      clinical13 <- read.csv("../Data/Uterus/clinical.csv", header = TRUE, sep = ",")
      # View(clinical13)
      # Reemplazar ''--' por NA
      clinical13[clinical13 == "'--"] <- NA
      # View(clinical13)
      # Solo filas con stage definida
      clinical_clean13 <- clinical13 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean13 <- clinical_clean13 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean13)
      # Num de pacientes
      num_unicos13 <- clinical_clean13 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos13)
      # Eliminar duplicados en clinical_clean13
      clinical_clean13_unique <- clinical_clean13 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final13 <- clinical_clean13_unique[!grepl("metastasis", clinical_clean13_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final13)
      # Num pacientes
      num_unicos13 <- clinical_final13 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos13)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal13 <- "../Data/Uterus"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista13 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas13 <- list.dirs(directorio_principal13, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta13 in carpetas13) {
        # Lista los archivos .maf en la carpeta
        archivos_maf13 <- list.files(carpeta13, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf13 in archivos_maf13) {
          # Lee el archivo .maf
          maf_data13 <- fread(archivo_maf13)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista13[[length(resultados_lista13) + 1]] <- maf_data13
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones13 <- rbindlist(resultados_lista13, fill = TRUE)
      # Dataframe final
      # View(mutaciones13)
      # Genes no especificados
      mutaciones_clean13 <- mutaciones13 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas13 <- unique(mutaciones13$Mutation_Status)
      print(mutaciones_unicas13)
      mutaciones_clean13 <- mutaciones_clean13 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean13 <- mutaciones_clean13 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean13 <- mutaciones_clean13 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean13)
      # Num de pacientes
      num_unicos13 <- mutaciones_clean13 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos13)
    }
    #Unir
    {
      resultado13 <- merge(mutaciones_clean13, clinical_final13, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data13 <- resultado13 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data13, columns) {
        for (col in columns) {
          data13 <- data13 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data13)
      }
      # Aplicar la función a los datos
      data13 <- clean_columns(data13, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data13 <- data13 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data13 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data13 <- data13 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data13)
      # Num de pacientes
      num_unicos13 <- data13 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos13)
      
      # Codificación de variables categóricas
      data_final13 <- data13 %>%
        mutate(classif_cancer = 13) %>%
        mutate(tipo_cancer = "Uterus")
      View(data_final13)
    }
    
  }
  #14. Skin
  {
    #Clinico
    {
      clinical14 <- read.csv("../Data/Skin/clinical.csv", header = TRUE, sep = ",")
      # View(clinical14)
      # Reemplazar ''--' por NA
      clinical14[clinical14 == "'--"] <- NA
      # View(clinical14)
      # Solo filas con stage definida
      clinical_clean14 <- clinical14 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean14 <- clinical_clean14 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean14)
      # Num de pacientes
      num_unicos14 <- clinical_clean14 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos14)
      # Eliminar duplicados en clinical_clean14
      clinical_clean14_unique <- clinical_clean14 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final14 <- clinical_clean14_unique[!grepl("metastasis", clinical_clean14_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final14)
      # Num pacientes
      num_unicos14 <- clinical_final14 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos14)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal14 <- "../Data/Skin"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista14 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas14 <- list.dirs(directorio_principal14, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta14 in carpetas14) {
        # Lista los archivos .maf en la carpeta
        archivos_maf14 <- list.files(carpeta14, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf14 in archivos_maf14) {
          # Lee el archivo .maf
          maf_data14 <- fread(archivo_maf14)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista14[[length(resultados_lista14) + 1]] <- maf_data14
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones14 <- rbindlist(resultados_lista14, fill = TRUE)
      # Dataframe final
      # View(mutaciones14)
      # Genes no especificados
      mutaciones_clean14 <- mutaciones14 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas14 <- unique(mutaciones14$Mutation_Status)
      print(mutaciones_unicas14)
      mutaciones_clean14 <- mutaciones_clean14 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean14 <- mutaciones_clean14 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean14 <- mutaciones_clean14 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean14)
      # Num de pacientes
      num_unicos14 <- mutaciones_clean14 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos14)
    }
    #Unir
    {
      resultado14 <- merge(mutaciones_clean14, clinical_final14, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data14 <- resultado14 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data14, columns) {
        for (col in columns) {
          data14 <- data14 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data14)
      }
      # Aplicar la función a los datos
      data14 <- clean_columns(data14, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data14 <- data14 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data14 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data14 <- data14 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data14)
      # Num de pacientes
      num_unicos14 <- data14 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos14)
      
      # Codificación de variables categóricas
      data_final14 <- data14 %>%
        mutate(classif_cancer = 14) %>%
        mutate(tipo_cancer = "Skin")
      View(data_final14)
    }
    
  }
  #15. Thyroid
  {
    #Clinico
    {
      clinical15 <- read.csv("../Data/Thyroid/clinical.csv", header = TRUE, sep = ",")
      # View(clinical15)
      # Reemplazar ''--' por NA
      clinical15[clinical15 == "'--"] <- NA
      # View(clinical15)
      # Solo filas con stage definida
      clinical_clean15 <- clinical15 
        #filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
      # Quedarse con las columnas importantes
      clinical_clean15 <- clinical_clean15 %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      # View(clinical_clean15)
      # Num de pacientes
      num_unicos15 <- clinical_clean15 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos15)
      # Eliminar duplicados en clinical_clean15
      clinical_clean15_unique <- clinical_clean15 %>%
        distinct(case_id, .keep_all = TRUE)
      # Eliminar metastasis
      clinical_final15 <- clinical_clean15_unique[!grepl("metastasis", clinical_clean15_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_final15)
      # Num pacientes
      num_unicos15 <- clinical_final15 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos15)
      
    }
    #Mutaciones
    {
      # Define el directorio principal donde están las carpetas
      directorio_principal15 <- "../Data/Thyroid"
      # Inicializa una lista para almacenar la información de los archivos .maf
      resultados_lista15 <- list()
      # Obtiene la lista de carpetas dentro del directorio principal
      carpetas15 <- list.dirs(directorio_principal15, full.names = TRUE, recursive = TRUE)
      # Recorre cada carpeta
      for (carpeta15 in carpetas15) {
        # Lista los archivos .maf en la carpeta
        archivos_maf15 <- list.files(carpeta15, pattern = "\\.maf$", full.names = TRUE)
        # Recorre cada archivo .maf
        for (archivo_maf15 in archivos_maf15) {
          # Lee el archivo .maf
          maf_data15 <- fread(archivo_maf15)  # fread para leer el archivo
          # Agrega la información a la lista
          resultados_lista15[[length(resultados_lista15) + 1]] <- maf_data15
        }
      }
      # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
      mutaciones15 <- rbindlist(resultados_lista15, fill = TRUE)
      # Dataframe final
      # View(mutaciones15)
      # Genes no especificados
      mutaciones_clean15 <- mutaciones15 %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somáticas
      mutaciones_unicas15 <- unique(mutaciones15$Mutation_Status)
      print(mutaciones_unicas15)
      mutaciones_clean15 <- mutaciones_clean15 %>%
        filter(!is.na(Mutation_Status))
      # Añadir columna con identificador
      mutaciones_clean15 <- mutaciones_clean15 %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      # Quedarse con las columnas importantes
      mutaciones_clean15 <- mutaciones_clean15 %>%
        select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT)
      #View(mutaciones_clean15)
      # Num de pacientes
      num_unicos15 <- mutaciones_clean15 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos15)
    }
    #Unir
    {
      resultado15 <- merge(mutaciones_clean15, clinical_final15, by = "case_id", all = FALSE)
      # View(resultado)
      # Eliminar duplicados
      data15 <- resultado15 %>%
        distinct(case_id, Identificador, .keep_all = TRUE)
      # Limpiar las columnas con vacíos
      columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
      # Función para limpiar las columnas
      clean_columns <- function(data15, columns) {
        for (col in columns) {
          data15 <- data15 %>%
            mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
            mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
        }
        return(data15)
      }
      # Aplicar la función a los datos
      data15 <- clean_columns(data15, columns_to_clean)
      # View(data)
      # Eliminar mutaciones que se repitan menos de 5 veces
      data15 <- data15 %>% filter(!is.na(Identificador))  # Filtrar los NA primero
      mutaciones_repetidas <- data15 %>% group_by(Identificador) %>% summarise(count = n(), .groups = 'drop') %>% filter(count >= 5)
      data15 <- data15 %>% semi_join(mutaciones_repetidas, by = "Identificador")
      # View(data15)
      # Num de pacientes
      num_unicos15 <- data15 %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos15)
      
      # Codificación de variables categóricas
      data_final15 <- data15 %>%
        mutate(classif_cancer = 15) %>%
        mutate(tipo_cancer = "Thyroid")
      View(data_final15)
    }
    
  }
  
  #Unir todas
  {
    # Unir los data_final de los diferentes conjuntos
    {
      data_final1$age_at_index <- as.character(data_final1$age_at_index)
      data_final2$age_at_index <- as.character(data_final2$age_at_index)
      data_final3$age_at_index <- as.character(data_final3$age_at_index)
      data_final4$age_at_index <- as.character(data_final4$age_at_index)
      data_final6$age_at_index <- as.character(data_final6$age_at_index)
      data_final7$age_at_index <- as.character(data_final7$age_at_index)
      data_final8$age_at_index <- as.character(data_final8$age_at_index)
      data_final9$age_at_index <- as.character(data_final9$age_at_index)
      data_final10$age_at_index <- as.character(data_final10$age_at_index)
      data_final11$age_at_index <- as.character(data_final11$age_at_index)
      data_final12$age_at_index <- as.character(data_final12$age_at_index)
      data_final13$age_at_index <- as.character(data_final13$age_at_index)
      data_final14$age_at_index <- as.character(data_final14$age_at_index)
      data_final15$age_at_index <- as.character(data_final15$age_at_index)
    }
    data_final <- bind_rows(data_final1, data_final2, data_final3, data_final4, data_final5, data_final6, 
                            data_final7, data_final8, data_final9, data_final10, data_final11, data_final12, 
                            data_final13, data_final14, data_final15)
    View(data_final)
    # Verificar si hay filas duplicadas
    duplicates <- data_final %>%
      group_by(case_id, Identificador, Hugo_Symbol) %>%
      summarise(count = n(), .groups = 'drop') %>%
      filter(count > 1)
    View(duplicates) 
    # Eliminar las filas duplicadas
    data_final <- data_final %>%
      distinct(case_id, Identificador, Hugo_Symbol, .keep_all = TRUE)
    View(data_final)
  }
  #Clasificacion binaria
  {
    #Se pone el numero del cancer de interes
  data_final <- data_final %>%
    mutate(classif_can = ifelse(classif_cancer %in% c("15"), 1, 0))
    View(data_final)
  }
}

##ANALISIS EXPLORATORIO
{
  # Pacientes por tipo de cáncer, con mutaciones que se repiten al menos 5 veces
  {
  counts <- data_final %>%
    group_by(tipo_cancer) %>%
    summarise(count = n_distinct(case_id))
  # View(counts)
  ggplot(counts, aes(x = tipo_cancer, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Pacientes en cada tipo de cáncer (mutaciones repetidas más de 5 veces)",
         x = "Tipo de Cáncer",
         y = "Conteo de Pacientes") +
    theme_minimal() +
    geom_text(aes(label = count), 
              position = position_stack(vjust = 0.5), 
              color = "white")
  }
  
  #Número de mutaciones por tipo de cáncer
  {
    # Filtrar las mutaciones que se repiten al menos 5 veces
    filtered_data <- data_final %>%
      group_by(Identificador) %>%
      filter(n() >= 5) %>%
      ungroup()
  mutaciones_por_cancer <- filtered_data %>%
    group_by(tipo_cancer) %>%
    summarise(num_mutaciones = n_distinct(Identificador))  # Contamos las mutaciones únicas
  print(mutaciones_por_cancer)
  }
  
  #Upset plot
  {
    library(UpSetR)
    library(tibble)
    #Todas las mutaciones
    {
    mutaciones_stage <- data_final %>%
      select(Identificador, tipo_cancer) %>%
      mutate(value = 1) %>%
      group_by(Identificador, tipo_cancer) %>%
      summarize(count = sum(value), .groups = 'drop') %>%
      pivot_wider(names_from = tipo_cancer, values_from = count, values_fill = list(count = 0))
    # columnas binarias
    mutaciones_binary <- as.data.frame(sapply(mutaciones_stage[-1], function(x) ifelse(x > 0, 1, 0)))
    mutaciones_ready <- cbind(Identificador = mutaciones_stage$Identificador, mutaciones_binary)
    upset(mutaciones_ready, sets = colnames(mutaciones_ready)[-1], order.by = "freq")
    #Num de mutaciones por etapa
    mutaciones_por_etapa <- data_final %>%
      group_by(tipo_cancer) %>%
      summarise(num_mutaciones = n_distinct(Identificador))
    print(mutaciones_por_etapa)
    }
    
    #Upset con mutaciones repetidas al menos 5 veces
    {
    mutaciones_filtradas <- data_final %>%
      group_by(Identificador) %>%
      summarise(num_pacientes = n_distinct(case_id), .groups = 'drop') %>%
      filter(num_pacientes >= 5) %>%
      select(Identificador)
    mutaciones_stage <- data_final %>%
      semi_join(mutaciones_filtradas, by = "Identificador") %>%
      select(Identificador, tipo_cancer) %>%
      mutate(value = 1) %>%
      group_by(Identificador, tipo_cancer) %>%
      summarize(count = sum(value), .groups = 'drop') %>%
      pivot_wider(names_from = tipo_cancer, values_from = count, values_fill = list(count = 0))
    mutaciones_binary <- as.data.frame(sapply(mutaciones_stage[-1], function(x) ifelse(x > 0, 1, 0)))
    mutaciones_ready <- cbind(Identificador = mutaciones_stage$Identificador, mutaciones_binary)
    #Upset plot
    #png("upset_plot.png", width = 800, height = 600)
    upset(mutaciones_ready, sets = colnames(mutaciones_ready)[-1], order.by = "freq")
    #dev.off()
    
    mutaciones_repetidas <- data_final %>%
      group_by(Identificador) %>%
      summarise(num_ocurrencias = n(), .groups = 'drop') %>%
      filter(num_ocurrencias >= 5) %>%
      select(Identificador)
    # Contar número de mutaciones por etapa solo con las filtradas
    mutaciones_por_etapa <- data_final %>%
      semi_join(mutaciones_repetidas, by = "Identificador") %>%
      group_by(STAGE) %>%
      summarise(num_mutaciones = n_distinct(Identificador), .groups = 'drop')
    print(mutaciones_por_etapa)
    }
    }
  
}

##MODELO
{
  #Organizacion apta para el modelo
  {
    #Quedarse con las mutaciones que se repiten al menos 5 veces
    mutation_counts <- table(data_final$Identificador)
    mutations_to_keep <- names(mutation_counts[mutation_counts >= 5])
    filtered_data_final <- data_final[data_final$Identificador %in% mutations_to_keep, ]
    #View(filtered_data_final)
    #DF con el tipo de cancer
    df_modelo <- data_final %>%
      select(classif_can, case_id) %>%
      distinct(case_id, .keep_all = TRUE) 
    #View(df_modelo)
    # DF con mutaciones unicas
    df_gen_mutado <- filtered_data_final %>%
      select(case_id,Identificador) %>%
      mutate(mutado = 1) %>%  # Asignar 1 si la mutación está presente
      distinct() %>%  # no haya duplicados
      pivot_wider(
        names_from = Identificador, 
        values_from = mutado,
        values_fill = list(mutado = 0)  # Rellenar con 0 si la mutación no está presente
      )
    #View(df_gen_mutado)
    # Combinar los dos dataframes
    data_modelo <- merge(df_gen_mutado, df_modelo, by = "case_id")
    View(data_modelo)
  }
  table(data_modelo$classif_can)
  #Entrenamiento y prueba
  {
    set.seed(135)
    trainIndex <- createDataPartition(data_modelo$classif_can, p = .7, 
                                      list = FALSE, 
                                      times = 1)
    df_train <- data_modelo[trainIndex, ]
    df_test  <- data_modelo[-trainIndex, ]
    
    x_train <- df_train %>% select(-case_id, -classif_can)
    y_train <- df_train$classif_can
    x_test <- df_test %>% select(-case_id, -classif_can)
    y_test <- df_test$classif_can
    y_train <- factor(y_train, levels = c("0", "1"))
    y_test <- factor(y_test, levels = c("0", "1"))
    
    
    # Evaluar el balanceo del conjunto de entrenamiento
    train_class_distribution <- table(y_train)
    print(train_class_distribution)
    # Evaluar el balanceo del conjunto de prueba
    test_class_distribution <- table(y_test)
    print(test_class_distribution)
  }
  #Modelo
  {
    ## RANDOM FOREST
    rf_model <- randomForest(x = x_train, y = as.factor(y_train), ntree = 100)
    print(rf_model)
    # Predicciones
    predictions <- predict(rf_model, x_test)
    # Evaluar el rendimiento
    confusionMatrix(predictions, as.factor(df_test$classif_can))
    confusion_mat <- confusionMatrix(predictions, as.factor(y_test))
    print(confusion_mat)
    
    #En caso de tener baja especificidad se puede usar este codigo, ajustando "threshold" (sube la especificidad mientras baja la sensibilidad)
    # Predicciones de probabilidades
    pred_probs <- predict(rf_model, x_test, type = "prob")
    # Umbral personalizado para la clase '1'
    threshold <- 0.2
    #probabilidad de la clase es mayor o igual al umbral
    predicted_classes <- ifelse(pred_probs[, 2] >= threshold, 1, 0)
    head(predicted_classes)
    predicted_classes <- factor(predicted_classes, levels = c(0, 1))
    #Matriz de confusion
    confusionMatrix(predicted_classes, as.factor(y_test))
    
    
    # Grafico de la matriz de confusion
    confusion_df <- as.data.frame(confusion_mat$table)
    ggplot(confusion_df, aes(x = Prediction, y = Reference)) +
      geom_tile(aes(fill = Freq), color = "white") +
      scale_fill_gradient(low = "white", high = "blue") +
      geom_text(aes(label = Freq), color = "black") +
      labs(title = "Matriz de Confusión (datos de prueba)", x = "Predicción", y = "Referencia")
    
    #Curva ROC y AUC
    y_test <- factor(y_test, levels = c(0, 1))  # "0" para etapa tardía, "1" para etapa temprana
    #Probabilidades de predicción para la clase positiva
    prob_predictions <- predict(rf_model, x_test, type = "prob")
    #Curva ROC
    png("curva_roc.png", width = 800, height = 800)
    roc_curve <- roc(y_test, prob_predictions[, 2], plot = TRUE, col = "blue", main = "Curva ROC (para cáncer de tiroides)", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
    dev.off()
    # Calcular el AUC (area bajo la curva)
    auc_value <- auc(roc_curve)
    print(paste("AUC:", auc_value))
    
    #Importancia de las variables
    importance_df <- as.data.frame(importance(rf_model))
    importance_df$Variable <- rownames(importance_df)
    top_importance_df <- importance_df %>%
      arrange(desc(MeanDecreaseGini)) %>%
      head(10)
    ggplot(top_importance_df, aes(x = reorder(Variable, MeanDecreaseGini), y = MeanDecreaseGini)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Importancia de las Variables", x = "Variable", y = "Importancia")
    #     ggsave("importancia_variables.png", width = 8, height = 6)
    
    
  }
}
#PCA
{
  #Entrenamiento
  {
    {
    #Quedarse con las mutaciones que se repiten al menos 5 veces
    mutation_counts <- table(data_final$Identificador)
    mutations_to_keep <- names(mutation_counts[mutation_counts >= 5])
    filtered_data_final <- data_final[data_final$Identificador %in% mutations_to_keep, ]
    #View(filtered_data_final)
    #DF con el tipo de cancer
    df_modelo <- data_final %>%
      select(tipo_cancer, case_id) %>%
      distinct(case_id, .keep_all = TRUE) 
    #View(df_modelo)
    # DF con mutaciones unicas
    df_gen_mutado <- filtered_data_final %>%
      select(case_id,Identificador) %>%
      mutate(mutado = 1) %>%  # Asignar 1 si la mutación está presente
      distinct() %>%  # no haya duplicados
      pivot_wider(
        names_from = Identificador, 
        values_from = mutado,
        values_fill = list(mutado = 0)  # Rellenar con 0 si la mutación no está presente
      )
    #View(df_gen_mutado)
    # Combinar los dos dataframes
    data_modelo <- merge(df_gen_mutado, df_modelo, by = "case_id")
    View(data_modelo)
  }
    
    # Seleccionar las mutaciones para realizar el PCA
    df_pca <- data_modelo %>%
      select(-case_id)  # Excluimos 'case_id' pero no 'tipo_cancer'
    # Estandarizar los datos de mutaciones (sin la columna tipo_cancer)
    df_pca_scaled <- scale(df_pca[, -which(names(df_pca) == "tipo_cancer")])
    # Realizar el PCA
    pca_result <- prcomp(df_pca_scaled, center = TRUE, scale. = TRUE)
    # Ver el resumen del PCA para entender la varianza explicada
    summary(pca_result)
    # Extraer las primeras dos componentes principales
    pca_scores <- pca_result$x[, 1:2]
    # Crear un dataframe con las componentes principales y la clasificación del cáncer (tipo_cancer)
    pca_df <- data.frame(pca_scores, tipo_cancer = data_modelo$tipo_cancer)
    # Ahora, 'pca_df' contiene las componentes principales junto con la clasificación de cáncer
    # Graficar los resultados del PCA, coloreados por tipo de cáncer
    colors <- c(brewer.pal(12, "Set3"), brewer.pal(3, "Set1"))
    ggplot(pca_df, aes(x = PC1, y = PC2, color = tipo_cancer)) +
      geom_point() +
      scale_color_manual(values = colors) +
      labs(title = "PCA: Conjunto de entrenamiento",
           x = "Componente Principal 1",
           y = "Componente Principal 2") +
      theme_minimal()
    
    
    {
    # Paso 1: Preprocesar los datos
    # Seleccionamos las mutaciones (columnas de mutaciones) y normalizamos los datos
    pca_data <- data_modelo %>%
      select(-case_id, -tipo_cancer)  # Excluir 'case_id' y 'tipo_cancer' para el PCA
    # Paso 2: Realizar el PCA
    # Normalizamos y realizamos el PCA sobre las mutaciones
    pca_result <- prcomp(pca_data, center = TRUE, scale. = TRUE)
    # Paso 3: Resumen de los resultados del PCA
    summary(pca_result)
    # Paso 4: Proporción de varianza explicada por cada componente
    explained_variance <- summary(pca_result)$importance[2, ]
    print(explained_variance)
    # Paso 5: Visualización de la varianza explicada
    fviz_eig(pca_result,  addlabels = TRUE)
    # Paso 6: Visualización de los individuos (muestras) en las dos primeras componentes principales
    fviz_pca_ind(pca_result, 
                 geom = "point", 
                 col.ind = data_modelo$tipo_cancer,  # Colorear por 'tipo_cancer'
                 palette = "Set3",  # Cambiar la paleta de colores para distinguir mejor
                 addEllipses = TRUE,  # Añadir elipses de concentración para cada tipo de cáncer
                 legend.title = "Cancer Type")
    # Graficar los resultados del PCA ajustando márgenes y tamaño del gráfico
    fviz_pca_ind(pca_result, 
                 geom = "point", 
                 col.ind = data_modelo$tipo_cancer,  # Colorear por 'tipo_cancer'
                 palette = "Set3",  # Cambiar la paleta de colores
                 addEllipses = TRUE,  # Añadir elipses de concentración
                 legend.title = "Cancer Type") +
      theme_minimal() +
      theme(
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5),  # Ajustar márgenes: top, right, bottom, left
        aspect.ratio = 1  # Ajustar la proporción del gráfico (aspect ratio)
      )
    
    # Paso 7: Visualización de las variables (mutaciones) en el espacio de componentes principales
    fviz_pca_var(pca_result, 
                 col.var = "contrib",  # Colorear las variables según su contribución a las componentes
                 gradient.cols = c("blue", "red"),  # Colores para las contribuciones (de azul a rojo)
                 repel = TRUE)  # Evitar superposición de etiquetas
      }
    
  }
}
##VALIDACION
{
  #Cargar datos de validacion
{
  {
    #OrigiMed, Nature 2022
    {
    #Datos tipo de cancer
    clinical_val1 <- read.csv("../Data/Validacion/OrigiMed, Nature 2022/data_clinical_sample.csv", header = TRUE, sep = ",")
    clinical_val1 <- clinical_val1 %>%
      select(SAMPLE_ID, CANCER_TYPE)
    #View(clinical_val1)
    #Datos de mutaciones
    mutaciones_val1<-read.csv("../Data/Validacion/OrigiMed, Nature 2022/data_mutations.csv", header = TRUE, sep = ",")
    #View(mutaciones_val1)
    mutaciones_val1 <- mutaciones_val1 %>%
      mutate(Identificador = paste0("chr", Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Tumor_Seq_Allele1, Tumor_Seq_Allele2)%>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)%>%
    mutate(Chromosome = paste0("chr", Chromosome))
    #View(mutaciones_val1)
    #Merge
    data_final_val1 <- merge(clinical_val1, mutaciones_val1, by = "SAMPLE_ID", all = FALSE)
    data_final_val1 <- data_final_val1 %>%
      filter(!is.na(CANCER_TYPE) & !is.na(Identificador))
    #View(data_final_val1)
    }
    
    #TCGA, Firehose Legacy (tiroides)
    {
    #Datos tipo de cancer
      clinical_val2 <- read.csv("../Data/Validacion/TCGA, Firehose Legacy/data_clinical_sample.csv", header = TRUE, sep = ",")
      clinical_val2 <- clinical_val2 %>%
        select(SAMPLE_ID, CANCER_TYPE)
      #View(clinical_val2)
    #Datos de mutaciones
      mutaciones_val2 <- read.csv("../Data/Validacion/TCGA, Firehose Legacy/data_mutations.csv", header = TRUE, sep = ",")
       mutaciones_val2 <- mutaciones_val2 %>%
        mutate(Identificador = paste0("chr", Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
        select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Tumor_Seq_Allele1, Tumor_Seq_Allele2) %>%
        rename(SAMPLE_ID = Tumor_Sample_Barcode)%>%
         mutate(Chromosome = paste0("chr", Chromosome))
      #View(mutaciones_val2)
      #Merge
      data_final_val2 <- merge(clinical_val2, mutaciones_val2, by = "SAMPLE_ID", all = FALSE)
      data_final_val2 <- data_final_val2 %>%
        filter(!is.na(CANCER_TYPE) & !is.na(Identificador))
      #View(data_final_val2)
    }
    
    #Unir archivos
    {
      data_final_val <- rbind(data_final_val1, data_final_val2)
      #View(data_final_val)
      # Verificar duplicados en las columnas SAMPLE_ID e Identificador
      duplicados <- data_final_val[duplicated(data_final_val[, c("SAMPLE_ID", "Identificador", "CANCER_TYPE")]), ]
      #View(duplicados)
      # Eliminar filas duplicadas
      data_final_val <- data_final_val %>%
        distinct(SAMPLE_ID, Identificador, CANCER_TYPE, .keep_all = TRUE)
      View(data_final_val)
    }
    #Archivo para pasar de GRCh37 a GRCh38 en https://genome.ucsc.edu/cgi-bin/hgLiftOver
    {
      #BED format
      data_genoma <- data_final_val %>%
        mutate(
          Chromosome = str_extract(Identificador, "chr\\d+"),  # Extrae el cromosoma (ej. chr12)
          Start_Position = as.numeric(str_extract(Identificador, "(?<=:)(\\d+)(?=_)")),  # Extrae la primera posición
          End_Position = as.numeric(str_extract(Identificador, "(?<=_)(\\d+)(?=_)"))  # Extrae la segunda posición
        ) %>%
        select(Chromosome, Start_Position, End_Position)  # Seleccionamos solo las nuevas columnas
      View(data_genoma)
      #Guardar archivo
      write.table(data_genoma, 
                  file = "data_genoma.bed",  
                  sep = "\t",               # Usar tabulaciones como delimitador
                  row.names = FALSE,        # No incluir números de filas
                  col.names = FALSE,        # No incluir nombres de columnas
                  quote = FALSE)            # No poner comillas alrededor de los valore
      #Archivo con posiciones anteriores (data_genoma)
      {
        data_genoma37 <- data_final_val %>%
          mutate(
            Start_Position_1 = Start_Position + 1,
            Identificador2 = paste(Chromosome, Start_Position_1, sep = ":"),
            Identificador2 = paste(Identificador2, End_Position, sep = "-") 
          ) %>%
          select(Identificador, Identificador2, Start_Position,Start_Position_1, End_Position, Chromosome, Tumor_Seq_Allele1, Tumor_Seq_Allele2)
        View(data_genoma37)
      }
    }
      #Archivo procesado GRCh38
      {
        bed_data <- read.table("../Data/Validacion/BED.bed", header = FALSE, sep = "\t")
          bed_data <- bed_data[, -ncol(bed_data)]
          colnames(bed_data) <- c("Chromosome2", "Start_Position2", "End_Position2", "Identificador2")
          View(bed_data)
          
        
        data_val <- merge(bed_data, data_genoma37, by = "Identificador2", all = TRUE)
        data_val <- distinct(data_val)
        data_val <- data_val %>%
          mutate(
            Identificador_val = paste(Chromosome2, 
                                      Start_Position2, 
                                      sep = ":") %>%
              paste(., Start_Position2, sep = "_") %>%
              paste(., Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "_")
          )%>%
          select(Identificador, Identificador_val)
        View(data_val)
        
        #Coincidencias de mutaciones entre el conjunto de validacion anterior y nuevo
        coincidencias <- inner_join(data_val, data_final_val, by = "Identificador")
        num_coincidencias <- nrow(coincidencias)
        print(num_coincidencias)
        
        #Merge con el dataset con Identificador anterior para que coincida con los pacientes
        data_validacion <- merge(data_final_val, data_val, by = "Identificador", all = FALSE)
        data_validacion <- data_validacion %>%
          select(SAMPLE_ID, Identificador_val, CANCER_TYPE, Hugo_Symbol)
        data_validacion <- data_validacion %>%
          rename(Identificador = Identificador_val)
        View(data_validacion)
      }
    }
}
  # Identificadores comunes en ambos conjuntos (validacion y entrenamiento)
  identificadores_comunes <- intersect(data_validacion$Identificador, data_final$Identificador)
  length(identificadores_comunes)
  
  # Ver las clases únicas en 'CANCER_TYPE'
  unique(data_validacion$CANCER_TYPE)
  
  #Seleccionar el cancer a validar
  {
  # Cancer de tiroides
  {
    data_validacion <- data_validacion %>%
      mutate(classif_can = ifelse(CANCER_TYPE %in% c("Thyroid Carcinoma", "Thyroid Cancer"), 1, 0))
  }
  
  # Cancer de riñón
  {
    data_validacion <- data_validacion %>%
      mutate(classif_can = ifelse(CANCER_TYPE %in% c("Kidney Renal Cell Carcinoma"), 1, 0))  
  }
  
  # Cancer de hígado
  {
    data_validacion <- data_validacion %>%
      mutate(classif_can = ifelse(CANCER_TYPE %in% c("Liver Hepatocellular Carcinoma"), 1, 0))   
  }
  }
  
  # Mutaciones totales del conjunto de validación
  num_mutaciones_unicas <- data_validacion %>%
    summarise(count = n_distinct(Identificador))
  print(num_mutaciones_unicas)
  
  # Organizar los datos de validación para el modelo (binario), agregando las mutaciones del conjunto de entrenamiento aunque no se tengan en el conjunto de validacion para poder hacer el proceso de predicciones
  {
    mutation_counts_val <- table(data_validacion$Identificador)
    mutations_to_keep_val <- names(mutation_counts_val[mutation_counts_val >= 5])
    filtered_data_validacion <- data_validacion[data_validacion$Identificador %in% mutations_to_keep_val, ]
    
    df_modelo_val <- data_validacion %>%
      select(classif_can, SAMPLE_ID) %>%
      distinct(SAMPLE_ID, .keep_all = TRUE) %>%
      rename(case_id = SAMPLE_ID)  
    
    df_gen_mutado_val <- filtered_data_validacion %>%
      select(SAMPLE_ID, Identificador) %>%
      mutate(mutado = 1) %>%
      distinct() %>%
      pivot_wider(
        names_from = Identificador, 
        values_from = mutado,
        values_fill = list(mutado = 0)  
      )
    
    mutation_columns <- colnames(data_modelo)[!colnames(data_modelo) %in% c("classif_can", "case_id")]
    
    missing_columns <- setdiff(mutation_columns, colnames(df_gen_mutado_val))
    
    if(length(missing_columns) > 0) {
      missing_columns_df <- data.frame(SAMPLE_ID = df_gen_mutado_val$SAMPLE_ID)
      for(col in missing_columns) {
        missing_columns_df[[col]] <- 0  
      }
      df_gen_mutado_val <- left_join(df_gen_mutado_val, missing_columns_df, by = "SAMPLE_ID")
    }
    
    df_gen_mutado_val <- df_gen_mutado_val %>%
      rename(case_id = SAMPLE_ID)
    
    data_modelo_val <- merge(df_gen_mutado_val, df_modelo_val, by = "case_id")
    
    View(data_modelo_val)
  }
  # Revisar la distribución de las clases en el conjunto de validación
  table(data_modelo_val$classif_can)
  # Predicciones de la validación
  {
    predictions <- predict(rf_model, newdata = data_modelo_val %>% select(-classif_can, -case_id))
    confusion_mat_val <- confusionMatrix(predictions, as.factor(data_modelo_val$classif_can))
    print(confusion_mat_val)
    
    #ROC
    prob_predictions <- predict(rf_model, newdata = data_modelo_val %>% select(-classif_can, -case_id), type = "prob")
    y_true <- as.factor(data_modelo_val$classif_can)  
    roc_curve_val <- roc(y_true, prob_predictions[, 2], plot = TRUE, col = "blue", main = "Curva ROC - Validación", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
    png("curva_roc_validacion.png", width = 800, height = 800)
    roc_curve_val <- roc(y_true, prob_predictions[, 2], plot = TRUE, col = "blue", main = "Curva ROC - Validación", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
    plot(roc_curve_val, col = "blue", main = "Curva ROC - Validación", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
    dev.off()
     auc_value_val <- auc(roc_curve_val)
    print(paste("AUC - Validación:", auc_value_val))
    
    #Ajustar el umbral si se tiene baja especificidad
    {
      pred_probs <- predict(rf_model, newdata = data_modelo_val %>% select(-classif_can, -case_id), type = "prob")
      threshold <- 0.1
      predicted_classes <- ifelse(pred_probs[, 2] >= threshold, 1, 0)
      predicted_classes <- factor(predicted_classes, levels = c(0, 1))  
      data_modelo_val$classif_can <- factor(data_modelo_val$classif_can, levels = c(0, 1))  
      confusion_mat_val <- confusionMatrix(predicted_classes, data_modelo_val$classif_can)
      print(confusion_mat_val)
      y_true <- as.factor(data_modelo_val$classif_can)
      
      # Generar la curva ROC
      png("curva_roc_validacion.png", width = 800, height = 800)
      roc_curve <- roc(y_true, pred_probs[, 2], plot = TRUE, col = "blue", main = "Curva ROC - Validación", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad" )
      plot(roc_curve, col = "blue", main = "Curva ROC - Validación", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
      dev.off()
      # Calcular AUC
      auc_value <- auc(roc_curve)
      print(paste("AUC - Validación:", auc_value))
                   
    }
    
    
  }
  

}
