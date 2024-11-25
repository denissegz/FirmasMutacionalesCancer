###FIRMAS MUTACIONALES DE ETAPAS TEMPRANAS DE CÁNCER
#Este codigo esta hecho para leer los archivos descomprimidos MAF del portal TCGA
#Se puede hacer para cualquier cancer, solo se cambia la direccion de los archivos clínicos y de mutaciones (clinical y directorio_principal), en la sección "CARGAR ARCHIVOS"
#En la parte de "Codificacion de etapas" en la sección "CARGAR ARCHIVOS" se puede escoger que etapas tomar en cuenta segun el estudio que se quiera hacer

#En caso de usar archivos de mutaciones de CBioPortal los archivos se cargaran simplemente con:
# mutaciones<-read.csv("../cBioPortal/MSK, JNCI 2021/data_mutations.csv", header = TRUE, sep = ",")


### LIBRERÍAS
{
  # Lista de paquetes que queremos cargar
  packages <- c("reshape2", "ggplot2", "dplyr", "caret", "tidyr", "rpart", "rpart.plot", 
                "randomForest", "e1071", "gbm", "class", "nnet", "glmnet", 
                "caretEnsemble", "car", "pROC", "data.table", "UpSetR", "tibble", 
                "isotree", "dbscan")
  
  # Comprobar si los paquetes están instalados y, si no lo están, instalarlos
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)  # Instalar el paquete si no está instalado
    }
    library(pkg, character.only = TRUE)  # Cargar el paquete
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
  library(e1071)
  library(gbm)
  library(class)
  library(nnet)
  library(glmnet)
  library(caretEnsemble)
  library(car)
  library(pROC)
  library(data.table)
  library(UpSetR)
  library(tibble)
  library(isotree)
  library(dbscan)
  
  
}

### CARGAR ARCHIVOS
{
#CLINICOS
  clinical <- read.csv("../Data/Colorectal/clinical.csv", header = TRUE, sep = ",")
{
  #View(clinical)
  # Reemplazar ''--' por NA
  clinical[clinical == "'--"] <- NA
  #Contar estapas clínica y paológica
  count_clinical <- sum(!is.na(clinical$ajcc_clinical_stage))
  print(count_clinical)
  count_path <- sum(!is.na(clinical$ajcc_pathologic_stage))
  print(count_path)
  #View(clinical)
  #Solo pacientes con stage definida
  clinical_clean <- clinical %>%
    filter(!is.na(ajcc_pathologic_stage) & ajcc_pathologic_stage != "Not Reported" & ajcc_pathologic_stage != "Stage X" & ajcc_pathologic_stage != "Unknown")
  #Diferencias de etapa clínica y patológica
  diferencias <- clinical_clean %>%
    filter(ajcc_pathologic_stage != ajcc_clinical_stage) %>%
    select(ajcc_pathologic_stage, ajcc_clinical_stage)
  print(diferencias)
  #Quedarse con las columnas importantes
  clinical_clean2 <- clinical_clean %>%
    select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
  #View(clinical_clean2)
  #Num de pacientes
  num_unicos <- clinical_clean2 %>%
    summarise(count = n_distinct(case_id))
  print(num_unicos)
  # Eliminar duplicados en clinical_clean2
  clinical_clean2_unique <- clinical_clean2 %>%
    distinct(case_id, .keep_all = TRUE)
  #Eliminar metastasis
  clinical_final <- clinical_clean2_unique[!grepl("metastasis", clinical_clean2_unique$classification_of_tumor, ignore.case = TRUE), ]
  #View(clinical_final)
  #num pacinetes
  num_unicos <- clinical_final %>%
    summarise(count = n_distinct(case_id))
  print(num_unicos)
}

#MUTACIONES SOMÁTICAS
  directorio_principal <- "../Data/Colorectal"
{
  # Inicializa una lista para almacenar la información de los archivos .maf
  resultados_lista <- list()
  # Obtiene la lista de carpetas dentro del directorio principal
  carpetas <- list.dirs(directorio_principal, full.names = TRUE, recursive = TRUE)
  # Recorre cada carpeta
  for (carpeta in carpetas) {
    # Lista los archivos .maf en la carpeta
    archivos_maf <- list.files(carpeta, pattern = "\\.maf$", full.names = TRUE)
    # Recorre cada archivo .maf
    for (archivo_maf in archivos_maf) {
      # Lee el archivo .maf
      maf_data <- fread(archivo_maf)  # fread para leer el archivo
      # Agrega la información a la lista
      resultados_lista[[length(resultados_lista) + 1]] <- maf_data
    }
  }
  # Combina todos los data.frames de la lista en uno solo, llenando columnas faltantes
  mutaciones <- rbindlist(resultados_lista, fill = TRUE)
  # Dataframe final
  #View(mutaciones)
  #Genes no especificados
  mutaciones_clean <- mutaciones %>%
    filter(Hugo_Symbol != "Unknown")
  # Solo filas con mutaciones somaticas
  mutaciones_unicas <- unique(mutaciones$Mutation_Status)
  print(mutaciones_unicas)
  mutaciones_clean <- mutaciones_clean %>%
    filter(!is.na(Mutation_Status))
  #Añadir columna con identificador
  mutaciones_clean <- mutaciones_clean %>%
    mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
  #Quedarse con las columnas importantes
  mutaciones_clean <- mutaciones_clean %>%
    select(case_id, Identificador, Hugo_Symbol, Gene, Entrez_Gene_Id, Consequence, Variant_Classification, Variant_Type, BIOTYPE, Mutation_Status, SIFT, PolyPhen, CLIN_SIG, IMPACT )
  #View(mutaciones_clean)
  #Num de pacientes
  num_unicos2 <- mutaciones_clean %>%
    summarise(count = n_distinct(case_id))
  print(num_unicos2)
}

#Unir mutaciones con clinico
{
resultado <- merge(mutaciones_clean, clinical_final, by = "case_id", all=FALSE)
#View(resultado)
#Eliminar duplicados
data <- resultado %>%
  distinct(case_id, Identificador, .keep_all = TRUE)
#Limpiar las columnas con vacíos
columns_to_clean <- c("PolyPhen", "SIFT", "CLIN_SIG")
# Función para limpiar las columnas
clean_columns <- function(data, columns) {
  for (col in columns) {
    data <- data %>%
      mutate(!!sym(col) := na_if(trimws(!!sym(col)), "")) %>%  # Reemplaza cadenas vacías con NA
      mutate(!!sym(col) := gsub("\\s*\\(.*?\\)", "", !!sym(col)))  # Elimina texto entre paréntesis
  }
  return(data)
}
# Aplicar la función a los datos
data_final <- clean_columns(data, columns_to_clean)
#View(data_final)
#Num de pacientes
num_unicos3 <- data_final %>%
  summarise(count = n_distinct(case_id))
print(num_unicos3)
}

#Codificacion de etapas
  {
#Etapa 1 y 2, 1, las demas 0
data_final <- data_final %>%
  mutate(classif_stage = ifelse(ajcc_pathologic_stage %in% c("Stage 0", "Stage I", "Stage IA", "Stage IA1", "Stage IA2","Stage IA3", "Stage IB", "Stage IIA", "Stage IIB", "Stage II", "Stage IIC"), 1, 0))
#Etapa 1, 1, las demas 0
data_final <- data_final %>%
  mutate(classif_stage = ifelse(ajcc_pathologic_stage %in% c("Stage 0", "Stage I", "Stage IA", "Stage IB", "Stage IA1", "Stage IA2","Stage IA3"), 1, 0))
#Etapa 2, 1, las demas 0
data_final <- data_final %>%
  mutate(classif_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC"), 1, 0))
#Etapa 3, 1, las demas 0
data_final <- data_final %>%
  mutate(classif_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IIIA", "Stage IIIB", "Stage III", "Stage IIIC" , "Stage IIIC1", "Stage IIIC2"), 1, 0))
#Etapa 4, 1, las demas 0
data_final <- data_final %>%
  mutate(classif_stage = ifelse(ajcc_pathologic_stage %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB"), 1, 0))
View(data_final)
}
}

### ANALISIS EXPLORATORIO
{
  # Contar mutaciones que estén en al menos 5 pacientes
  {
  mutaciones_repetidas <- data_final %>%
    group_by(Identificador) %>%  # Agrupar por mutación (Identificador)
    summarise(pacientes_unicos = n_distinct(case_id)) %>%  # Contar pacientes únicos
    filter(pacientes_unicos >= 5)  # Filtrar mutaciones presentes en al menos 5 pacientes
  numero_mutaciones <- nrow(mutaciones_repetidas)
  print(numero_mutaciones)
  }
  #Contar mutaciones únicas
  {
  mutaciones_unicas <- n_distinct(data_final$Identificador)
  print(mutaciones_unicas)
  }
  #Pacientes por etapa
  {
  counts <- data_final %>%
    select(ajcc_pathologic_stage, case_id) %>%
    distinct(case_id, .keep_all = TRUE) %>%
    group_by(ajcc_pathologic_stage, case_id) %>%
    summarise(count = n())
  #View(counts)
  ggplot(counts, aes(x = ajcc_pathologic_stage)) +
    geom_bar(fill = "steelblue") +
    labs(title = "Pacientes en cada etapa patológica (Pulmón)",
         x = "Etapa Patológica AJCC",
         y = "Conteo") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    geom_text(stat = "count", aes(label = ..count..), 
              position = position_stack(vjust = 0.5), 
              color = "white")
  #Agrupados por las 4 etapas
  counts <- counts %>%
    mutate(stage_group = case_when(
      ajcc_pathologic_stage == "Stage 0" ~ "Stage 0",
      ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IA1", "Stage IA2","Stage IA3", "Stage IB") ~ "Stage I",
      ajcc_pathologic_stage %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC") ~ "Stage II",
      ajcc_pathologic_stage %in% c("Stage IIIC", "Stage IIIB", "Stage III", "Stage IIIA", "Stage IIIC1", "Stage IIIC2") ~ "Stage III",
      ajcc_pathologic_stage %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB") ~ "Stage IV",
      TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
    ))
  counts <- counts %>%
    group_by(stage_group) %>%
    summarise(count = n())
  #View(data_graf)
  ggplot(counts, aes(x = stage_group, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Pacientes en cada etapa patológica (Glándula adrenal)",
         x = "Etapa Patológica AJCC",
         y = "Conteo") +
    theme_minimal() +
    geom_text(aes(label = count), 
              position = position_stack(vjust = 0.5), 
              color = "white")
  }
  #Upset plot entre etapas
  {
    data_graf <- data_final %>%
      mutate(stage_group = case_when(
        ajcc_pathologic_stage == "Stage 0" ~ "Stage 0",
        ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage IA1", "Stage IA2", "Stage IA3") ~ "Stage I",
        ajcc_pathologic_stage %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC") ~ "Stage II",
        ajcc_pathologic_stage %in% c("Stage IIIC", "Stage IIIB", "Stage III", "Stage IIIA", "Stage IIIC1", "Stage IIIC2") ~ "Stage III",
        ajcc_pathologic_stage %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB") ~ "Stage IV",
        TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
      ))
    #Todas las mutaciones
    {
      mutaciones_stage <- data_graf %>%
        select(Identificador, stage_group) %>%
        mutate(value = 1) %>%
        group_by(Identificador, stage_group) %>%
        summarize(count = sum(value), .groups = 'drop') %>%
        pivot_wider(names_from = stage_group, values_from = count, values_fill = list(count = 0))
      # columnas binarias
      mutaciones_binary <- as.data.frame(sapply(mutaciones_stage[-1], function(x) ifelse(x > 0, 1, 0)))
      mutaciones_ready <- cbind(Identificador = mutaciones_stage$Identificador, mutaciones_binary)
      upset(mutaciones_ready, sets = colnames(mutaciones_ready)[-1], order.by = "freq")
      #Num de mutaciones por etapa
      mutaciones_por_etapa <- data_graf %>%
        group_by(stage_group) %>%
        summarise(num_mutaciones = n_distinct(Identificador))
      print(mutaciones_por_etapa)
      
    }
    
    #Upset con mutaciones repetidas al menos 5 veces
    {
      mutaciones_filtradas <- data_graf %>%
        group_by(Identificador) %>%
        summarise(num_pacientes = n_distinct(case_id), .groups = 'drop') %>%
        filter(num_pacientes >= 5) %>%
        select(Identificador)
      mutaciones_stage <- data_graf %>%
        semi_join(mutaciones_filtradas, by = "Identificador") %>%
        select(Identificador, stage_group) %>%
        mutate(value = 1) %>%
        group_by(Identificador, stage_group) %>%
        summarize(count = sum(value), .groups = 'drop') %>%
        pivot_wider(names_from = stage_group, values_from = count, values_fill = list(count = 0))
      mutaciones_binary <- as.data.frame(sapply(mutaciones_stage[-1], function(x) ifelse(x > 0, 1, 0)))
      mutaciones_ready <- cbind(Identificador = mutaciones_stage$Identificador, mutaciones_binary)
      #Upset plot
      #png("upset_plot.png", width = 800, height = 600)
      upset(mutaciones_ready, sets = colnames(mutaciones_ready)[-1], order.by = "freq")
      #dev.off()
      
      
      mutaciones_repetidas <- data_graf %>%
        group_by(Identificador) %>%
        summarise(num_ocurrencias = n(), .groups = 'drop') %>%
        filter(num_ocurrencias >= 5) %>%
        select(Identificador)
      # Contar número de mutaciones por etapa solo con las filtradas
      mutaciones_por_etapa <- data_graf %>%
        semi_join(mutaciones_repetidas, by = "Identificador") %>%
        group_by(stage_group) %>%
        summarise(num_mutaciones = n_distinct(Identificador), .groups = 'drop')
      
      print(mutaciones_por_etapa)
    }
    }
  #Graficas de mutaciones agrupadas segun diferentes caracteristicas
  {
    # Según PolyPhen
    df_freq2 <- data_final %>%
      filter(!is.na(PolyPhen) & PolyPhen != "unknown")  # Sin NA ni unknown
    df_freq2 <- df_freq2 %>%
      group_by(Hugo_Symbol, PolyPhen, Identificador) %>%
      summarise(Frecuencia = n(), .groups = 'drop') 
    #View(df_freq2)
    # Filtrar los genes más frecuentes
    top_genes2 <- df_freq2 %>%
      group_by(Identificador) %>%
      summarise(Total_Frecuencia = sum(Frecuencia)) %>%
      top_n(50, Total_Frecuencia) %>%
      pull(Identificador)
    df_top_freq2 <- df_freq2 %>%
      filter(Identificador %in% top_genes2)
    # Histograma
    colors <- c("probably_damaging" = "#FF8C00",  
                "possibly_damaging" = "#FFD700",
                "benign" = "#40E0D0") 
    ggplot(df_top_freq2, aes(x = Hugo_Symbol, y = Frecuencia, fill = PolyPhen)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Genes mutados agrupados por PolyPhen",
           x = "Genes",
           y = "Frecuencia",
           fill = "Polymorphism Phenotyping (PolyPhen)") +
      scale_fill_manual(values = colors) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
    
    
    # Según SIFT
    df_freq3 <- data_final %>%
      filter(!is.na(SIFT) & SIFT != "unknown")  # Sin NA ni unknown
    df_freq3 <- df_freq3 %>%
      group_by(Hugo_Symbol, SIFT, Identificador) %>%
      summarise(Frecuencia = n(), .groups = 'drop') 
    #View(df_freq3)
    # Filtrar los genes más frecuentes
    top_genes3 <- df_freq3 %>%
      group_by(Identificador) %>%
      summarise(Total_Frecuencia = sum(Frecuencia)) %>%
      top_n(50, Total_Frecuencia) %>%
      pull(Identificador)
    df_top_freq3 <- df_freq3 %>%
      filter(Identificador %in% top_genes3)
    # Histograma
    ggplot(df_top_freq3, aes(x = Hugo_Symbol, y = Frecuencia, fill = SIFT)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Genes mutados agrupados por SIFT",
           x = "Genes",
           y = "Frecuencia",
           fill = "Sorting Intolerant From Tolerant (SIFT)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    # Según CLIN_SIG
    df_freq4 <- data_final %>%
      filter(!is.na(CLIN_SIG))  # Sin NA 
    df_freq4 <- df_freq4 %>%
      group_by(Hugo_Symbol, CLIN_SIG, Identificador) %>%
      summarise(Frecuencia = n(), .groups = 'drop')
    #View(df_freq4)
    # Filtrar los genes más frecuentes
    top_genes4 <- df_freq4 %>%
      group_by(Identificador) %>%
      summarise(Total_Frecuencia = sum(Frecuencia), .groups = 'drop') %>%
      top_n(60, Total_Frecuencia) %>%
      pull(Identificador)
    df_top_freq4 <- df_freq4 %>%
      filter(Identificador %in% top_genes4)
    # Histograma
    ggplot(df_top_freq4, aes(x = Hugo_Symbol, y = Frecuencia, fill = CLIN_SIG)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Genes mutados agrupados por CLIN_SIG",
           x = "Genes",
           y = "Frecuencia",
           fill = "Clinical significance (CLIN_SIG)") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    # Según Variant_Classification
    df_freq5 <- data_final %>%
      group_by(Hugo_Symbol, Variant_Classification, Identificador) %>%
      summarise(Frecuencia = n(), .groups = 'drop')
    #View(df_freq5)
    # Filtrar los genes más frecuentes
    top_genes5 <- df_freq5 %>%
      group_by(Identificador) %>%
      summarise(Total_Frecuencia = sum(Frecuencia), .groups = 'drop') %>%
      top_n(50, Total_Frecuencia) %>%
      pull(Identificador)
    df_top_freq5 <- df_freq5 %>%
      filter(Identificador %in% top_genes5)
    # Histograma
    ggplot(df_top_freq5, aes(x = Hugo_Symbol, y = Frecuencia, fill = Variant_Classification)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(title = "Genes mutados agrupados por la clasificación de la variante",
           x = "Genes",
           y = "Frecuencia",
           fill = "Variant Classification") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  }
}

### CLASIFICADOR
{
  #Organización apta para el modelo
  {
    #Quedarse con las mutaciones que se repiten al menos 5 veces
    mutation_counts <- table(data_final$Identificador)
    mutations_to_keep <- names(mutation_counts[mutation_counts >= 5])
    filtered_data_final <- data_final[data_final$Identificador %in% mutations_to_keep, ]
    #View(filtered_data_final)
    #DF con la etapa
    df_modelo <- filtered_data_final %>%
      select(classif_stage, case_id) %>%
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
  table(data_modelo$classif_stage)
  #Entrenamiento y prueba
  {
    set.seed(135)
    trainIndex <- createDataPartition(data_modelo$classif_stage, p = .7, 
                                      list = FALSE, 
                                      times = 1)
    df_train <- data_modelo[trainIndex, ]
    df_test  <- data_modelo[-trainIndex, ]
    
    x_train <- df_train %>% select(-case_id, -classif_stage)
    y_train <- df_train$classif_stage
    x_test <- df_test %>% select(-case_id, -classif_stage)
    y_test <- df_test$classif_stage
    #y_test <- factor(y_test, levels = c("0", "1"))
    
    
    
    # Evaluar el balanceo del conjunto de entrenamiento
    train_class_distribution <- table(y_train)
    print(train_class_distribution)
    # Evaluar el balanceo del conjunto de prueba
    test_class_distribution <- table(y_test)
    print(test_class_distribution)
  }

  #Modelo
  {
    # RANDOM FOREST
    rf_model <- randomForest(x = x_train, y = as.factor(y_train), ntree = 100)
    print(rf_model)
    # Predicciones
    predictions <- predict(rf_model, x_test)
    # Evaluar el rendimiento
    confusionMatrix(predictions, as.factor(df_test$classif_stage))
    confusion_mat <- confusionMatrix(predictions, as.factor(y_test))
    
    #Si se obtiene una especificdad baja se pueda ajustar el umbral para que aumente, a la vez que disminuye la sensibilidad
    # Predicciones de probabilidad
    prob_predictions <- predict(rf_model, x_test, type = "prob")
    # Establecer el umbral
    umbral <- 0.3
    y_pred_class <- ifelse(prob_predictions[,2] >= umbral, "1", "0")
    # Convertir las predicciones a factores
    y_pred_class <- factor(y_pred_class, levels = c("0", "1"))
    # Evaluar el rendimiento
    confusionMatrix(y_pred_class, as.factor(y_test))
    
    
    # Gráfico de la matriz de confusión
    confusion_df <- as.data.frame(confusion_mat$table)
    ggplot(confusion_df, aes(x = Prediction, y = Reference)) +
      geom_tile(aes(fill = Freq), color = "white") +
      scale_fill_gradient(low = "white", high = "blue") +
      geom_text(aes(label = Freq), color = "black") +
      labs(title = "Matriz de Confusión (datos de prueba)", x = "Predicción", y = "Referencia")
    
    #Curva ROC y AUC
    y_test <- factor(y_test, levels = c(0, 1))  # "0" para etapa tardía, "1" para etapa temprana
    #Probabilidades de predicción para la clase positiva (por ejemplo, "1" es la clase positiva)
    prob_predictions <- predict(rf_model, x_test, type = "prob")
    #Curva ROC (se guarda en archivos)
    png("curva_roc.png", width = 800, height = 800)
    roc_curve <- roc(y_test, prob_predictions[, 2], plot = TRUE, col = "blue", main = "Curva ROC del Modelo", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
    dev.off()
    # Calcular el AUC (Área bajo la curva)
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
    # Este ultimo renglon es para guardarl en archivos la grafica en caso que RStudio diga que es muy grande
    #     ggsave("importancia_variables.png", width = 8, height = 6)
  
  }
}
