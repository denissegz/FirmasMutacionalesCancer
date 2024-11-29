##LIBRERÍAS
{
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
  library(stringr)
  library(isotree)
  library(dbscan)
  
  
}

##CARGAR ARCHIVOS
{
  
#MSK, JNCI 2021 
  {
    #Datos clínicos (etapa)
  clinical1 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, JNCI 2021/data_clinical_patient.csv", header = TRUE, sep = ",")
  #View(clinical1)
  #Stage definida
  clinical1 <- clinical1 %>%
    select(PATIENT_ID, STAGE_AT_DX)%>%
    filter(!is.na(STAGE_AT_DX))
  #View(clinical1)
  
  #Datos de muestra (para unir etapa y mutación)
  sample1<-read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, JNCI 2021/data_clinical_sample.csv", header = TRUE, sep = ",")
  #View(sample1)
  #Sin metastasis
  sample1 <- sample1[!grepl("Metastasis", sample1$SAMPLE_TYPE, ignore.case = TRUE), ]
  sample1 <- sample1[grepl("Colon Adenocarcinoma", sample1$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
  sample1 <- sample1%>%
    select(PATIENT_ID, SAMPLE_ID)
  #View(sample1)
  
  #Datos de mutaciones
   mutaciones1<-read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, JNCI 2021/data_mutations.csv", header = TRUE, sep = ",")
   #View(mutaciones1) 
   mutaciones1 <- mutaciones1 %>%
     mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))%>%
     select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode)%>%
   rename(SAMPLE_ID = Tumor_Sample_Barcode)
   #View(mutaciones1)
   
#MERGE
   data1 <- merge(sample1, mutaciones1, by = "SAMPLE_ID")
   data1 <- merge(data1, clinical1, by = "PATIENT_ID")%>%
   rename(STAGE = STAGE_AT_DX)
   View(data1)
   
   }
  #MSK, JCO Precis Oncol 2022
  {
    #Datos clínicos (etapa)
    clinical2 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, JCO Precis Oncol 2022/data_clinical_patient.csv", header = TRUE, sep = ",")
    #View(clinical2)
    #Sin metastasis
    clinical2 <- clinical2[!grepl("Yes", clinical2$METASTATIC, ignore.case = TRUE), ]
    clinical2 <- clinical2 %>%
      mutate(STAGE_AT_DIAGNOSIS = case_when(
        STAGE_AT_DIAGNOSIS == 1 ~ "I",
        STAGE_AT_DIAGNOSIS == 2 ~ "II",
        STAGE_AT_DIAGNOSIS == 3 ~ "III",
        STAGE_AT_DIAGNOSIS == 4 ~ "IV",
        TRUE ~ as.character(STAGE_AT_DIAGNOSIS)  # Mantener otros valores sin cambios
      ))
    clinical2 <- clinical2 %>%
      select(PATIENT_ID, STAGE_AT_DIAGNOSIS)
    #View(clinical2)
    
    #Datos de muestra (para unir etapa y mutación)
    sample2<-read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, JCO Precis Oncol 2022/data_clinical_sample.csv", header = TRUE, sep = ",")
   # View(sample2)
    #Adenocarcinoma
    sample2 <- sample2[grepl("Colon Adenocarcinoma", sample2$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    #Sin metastasis
    sample2 <- sample2%>%
      select(PATIENT_ID, SAMPLE_ID)
   # View(sample2)
    
    #Datos de mutaciones
    mutaciones2<-read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, JCO Precis Oncol 2022/data_mutations.csv", header = TRUE, sep = ",")
   # View(mutaciones2) 
    mutaciones2 <- mutaciones2 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))%>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode)%>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
    #View(mutaciones2)
    
    #MERGE
    data2 <- merge(sample2, mutaciones2, by = "SAMPLE_ID")
    data2 <- merge(data2, clinical2, by = "PATIENT_ID")%>%
      rename(STAGE = STAGE_AT_DIAGNOSIS)
    View(data2)
}
  #MSK, Gastroenterology 2020
  {
    # Datos clínicos (etapa)
    clinical3 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, Gastroenterology 2020/data_clinical_patient.csv", header = TRUE, sep = ",")
    #View(clinical3)
    #Sin metastasis
    clinical3 <- clinical3 %>%
      filter(
        METASTASIS_LIVER == "No" &
          METASTASIS_LUNG == "No" &
          METASTASIS_BONE == "No" &
          METASTASIS_OTHER == "No"
      )
    clinical3 <- clinical3 %>%
      select(PATIENT_ID)
    #View(clinical3)
    
    # Datos de muestra (para unir etapa y mutación)
    sample3 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, Gastroenterology 2020/data_clinical_sample.csv", header = TRUE, sep = ",")
    #View(sample3)
    #Adenoca
    sample3 <- sample3[grepl("Colon Adenocarcinoma", sample3$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    sample3 <- sample3 %>%
      select(PATIENT_ID, SAMPLE_ID, STAGE_AT_DIAGNOSIS)%>%
      filter(!is.na(STAGE_AT_DIAGNOSIS))
   # View(sample3)
    
    # Datos de mutaciones
    mutaciones3 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/MSK, Gastroenterology 2020/data_mutations.csv", header = TRUE, sep = ",")
    #View(mutaciones3)
    mutaciones3 <- mutaciones3 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
   # View(mutaciones3)
    
    # MERGE
    data3 <- merge(sample3, mutaciones3, by = "SAMPLE_ID")
    data3 <- merge(data3, clinical3, by = "PATIENT_ID")%>%
      rename(STAGE = STAGE_AT_DIAGNOSIS)
    View(data3)
  }
  #CAS Shanghai, Cancer Cell 2020
  {
    # Datos clínicos (etapa)
    clinical4 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/CAS Shanghai, Cancer Cell 2020/data_clinical_patient.csv", header = TRUE, sep = ",")
    #View(clinical4)
    # Stage definida
    clinical4 <- clinical4 %>%
      select(PATIENT_ID, STAGE) %>%
      filter(!is.na(STAGE))
    #View(clinical4)

    # Datos de muestra (para unir etapa y mutación)
    sample4 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/CAS Shanghai, Cancer Cell 2020/data_clinical_sample.csv", header = TRUE, sep = ",")
    #View(sample4)
    #Sin metastasis
    sample4 <- sample4[grepl("Non-Metastatic", sample4$METASTATIC_STATUS, ignore.case = TRUE), ]
    # adenoca
    sample4 <- sample4[grepl("Colorectal Adenocarcinoma", sample4$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    sample4 <- sample4 %>%
      select(PATIENT_ID, SAMPLE_ID)
   # View(sample4)
    
    # Datos de mutaciones
    mutaciones4 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/CAS Shanghai, Cancer Cell 2020/data_mutations.csv", header = TRUE, sep = ",")
   # View(mutaciones4)
    mutaciones4 <- mutaciones4 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
    #View(mutaciones4)
    
    # MERGE
    data4 <- merge(sample4, mutaciones4, by = "SAMPLE_ID")
    data4 <- merge(data4, clinical4, by = "PATIENT_ID")
    View(data4)
    
  }
  #TCGA, PanCancer Atlas
  {
    # Datos clínicos (etapa)
    clinical5 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, PanCancer Atlas/data_clinical_patient.csv", header = TRUE, sep = ",")
    #View(clinical5)
    # Stage definida
    clinical5 <- clinical5 %>%
    mutate(STAGE = case_when(
      AJCC_PATHOLOGIC_TUMOR_STAGE == "STAGE 0" ~ "0",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE I", "STAGE IA") ~ "I",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IIA", "STAGE IIB", "STAGE II", "STAGE IIC") ~ "II",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IIIC", "STAGE IIIB", "STAGE III", "STAGE IIIA") ~ "III",
      AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("STAGE IVA", "STAGE IV", "STAGE IVC", "STAGE IVB") ~ "IV",
      TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
    ))
    clinical5 <- clinical5[!grepl("Otro", clinical5$STAGE, ignore.case = TRUE), ]
    clinical5 <- clinical5 %>%
      select(PATIENT_ID, STAGE) %>%
      filter(!is.na(STAGE))
   # View(clinical5)
    
    # Datos de muestra (para unir etapa y mutación)
    sample5 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, PanCancer Atlas/data_clinical_sample.csv", header = TRUE, sep = ",")
   # View(sample5)
    # Sin metastasis
    sample5 <- sample5[!grepl("Metastasis", sample5$SAMPLE_TYPE, ignore.case = TRUE), ]
    #Adeno
    sample5 <- sample5[grepl("Colon Adenocarcinoma", sample5$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    sample5 <- sample5 %>%
      select(PATIENT_ID, SAMPLE_ID)
   # View(sample5)
    
    # Datos de mutaciones
    mutaciones5 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, PanCancer Atlas/data_mutations.csv", header = TRUE, sep = ",")
   # View(mutaciones5)
    mutaciones5 <- mutaciones5 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
   # View(mutaciones5)
    
    # MERGE
    data5 <- merge(sample5, mutaciones5, by = "SAMPLE_ID")
    data5 <- merge(data5, clinical5, by = "PATIENT_ID")
    View(data5)
    
  }
  #TCGA, Firehose Legacy
  {
    # Datos clínicos (etapa)
    clinical6 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, Firehose Legacy/data_clinical_patient.csv", header = TRUE, sep = ",")
   # View(clinical6)
    # Stage definida
    clinical6 <- clinical6 %>%
      mutate(STAGE = case_when(
        AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage 0" ~ "0",
        AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage I", "Stage IA") ~ "I",
        AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC") ~ "II",
        AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage IIIC", "Stage IIIB", "Stage III", "Stage IIIA") ~ "III",
        AJCC_PATHOLOGIC_TUMOR_STAGE %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB") ~ "IV",
        TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
      ))
    clinical6 <- clinical6[!grepl("Otro", clinical6$STAGE, ignore.case = TRUE), ]
    clinical6 <- clinical6 %>%
      select(PATIENT_ID, STAGE) %>%
      filter(!is.na(STAGE))
   # View(clinical6)
    
    # Datos de muestra (para unir etapa y mutación)
    sample6 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, Firehose Legacy/data_clinical_sample.csv", header = TRUE, sep = ",")
   # View(sample6)
    # Sin metastasis
    sample6 <- sample6[!grepl("Metastasis", sample6$SAMPLE_TYPE, ignore.case = TRUE), ]
    # Adeno
    sample6 <- sample6[grepl("Colon Adenocarcinoma", sample6$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    sample6 <- sample6 %>%
      select(PATIENT_ID, SAMPLE_ID)
    #View(sample6)
    
    # Datos de mutaciones
    mutaciones6 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, Firehose Legacy/data_mutations.csv", header = TRUE, sep = ",")
    #View(mutaciones6)
    mutaciones6 <- mutaciones6 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
   # View(mutaciones6)
    
    # MERGE
    data6 <- merge(sample6, mutaciones6, by = "SAMPLE_ID")
    data6 <- merge(data6, clinical6, by = "PATIENT_ID")
    View(data6)
    
        
  }
  #TCGA, Nature 2012
  {
    # Datos clínicos (etapa)
   # clinical7 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, Nature 2012/data_clinical_patient.csv", header = TRUE, sep = ",")
    # View(clinical7)
    
    # Datos de muestra (para unir etapa y mutación)
    sample7 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, Nature 2012/data_clinical_sample.csv", header = TRUE, sep = ",")
    # View(sample7)
     sample7 <- sample7 %>%
       mutate(STAGE = case_when(
         TUMOR_STAGE_2009 == "Stage 0" ~ "0",
         TUMOR_STAGE_2009 %in% c("Stage I", "Stage IA") ~ "I",
         TUMOR_STAGE_2009 %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC") ~ "II",
         TUMOR_STAGE_2009 %in% c("Stage IIIC", "Stage IIIB", "Stage III", "Stage IIIA") ~ "III",
         TUMOR_STAGE_2009 %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB") ~ "IV",
         TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
       ))
     sample7 <- sample7[!grepl("Otro", sample7$STAGE, ignore.case = TRUE), ]
    # Sin metastasis
    sample7 <- sample7[!grepl("Metastasis", sample7$SAMPLE_TYPE, ignore.case = TRUE), ]
    # Adeno
    sample7 <- sample7[grepl("Colon Adenocarcinoma", sample7$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    sample7 <- sample7 %>%
      select(PATIENT_ID, SAMPLE_ID, STAGE)
    # View(sample7)
    
    # Datos de mutaciones
    mutaciones7 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/TCGA, Nature 2012/data_mutations.csv", header = TRUE, sep = ",")
    # View(mutaciones7)
    mutaciones7 <- mutaciones7 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
    # View(mutaciones7)
    
    # MERGE
    data7 <- merge(sample7, mutaciones7, by = "SAMPLE_ID")
    View(data7)
    
  }
  #DFCI, Cell Reports 2016
  {
    # Datos clínicos (etapa)
    #clinical8 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/DFCI, Cell Reports 2016/data_clinical_patient.csv", header = TRUE, sep = ",")
    #View(clinical8)
    
    # Datos de muestra (para unir etapa y mutación)
    sample8 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/DFCI, Cell Reports 2016/data_clinical_sample.csv", header = TRUE, sep = ",")
    # View(sample8)
    #NO ESPECIFICA METASTASIS
    # Adeno
    sample8 <- sample8[grepl("Colorectal Adenocarcinoma", sample8$CANCER_TYPE_DETAILED, ignore.case = TRUE), ]
    sample8 <- sample8 %>%
      select(PATIENT_ID, SAMPLE_ID, TUMOR_STAGE)%>%
      filter(!is.na(TUMOR_STAGE))%>%
      rename(STAGE = TUMOR_STAGE)
    # View(sample8)
    
    # Datos de mutaciones
    mutaciones8 <- read.csv("C:/Users/Denisse González/Downloads/cBioPortal/DFCI, Cell Reports 2016/data_mutations.csv", header = TRUE, sep = ",")
    # View(mutaciones8)
    mutaciones8 <- mutaciones8 %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2)) %>%
      select(Identificador, Hugo_Symbol, Tumor_Sample_Barcode) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)
    # View(mutaciones8)
    
    # MERGE
    data8 <- merge(sample8, mutaciones8, by = "SAMPLE_ID")
    View(data8)
    
  }
  #TCGA Genomic Data Commons Data Portal
  {
    #Clinico
    {
    clinical <- read.csv("C:/Users/Denisse González/Downloads/Data/Colorectal/clinical.csv", header = TRUE, sep = ",")
    #View(clinical)
    # Reemplazar ''--' por NA
    clinical[clinical == "'--"] <- NA
    # Solo filas con stage definida
    clinical_clean <- clinical %>%
      filter(!is.na(ajcc_pathologic_stage))
    #Quedarse con las columnas importantes
    clinical_clean2 <- clinical_clean %>%
      select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
    #View(clinical_clean2)
    # Eliminar duplicados en clinical_clean2
    clinical_clean2_unique <- clinical_clean2 %>%
      distinct(case_id, .keep_all = TRUE)
    #Eliminar metastasis
    clinical_clean2_unique <- clinical_clean2_unique[!grepl("metastasis", clinical_clean2_unique$classification_of_tumor, ignore.case = TRUE), ]
    #View(clinical_clean2_unique)
    #Qudarse con Adenocarcinoma, NOS
    clinical_final <- clinical_clean2_unique %>%
      filter(primary_diagnosis %in% c("Adenocarcinoma, NOS"))
    clinical_final <- clinical_final %>%
      mutate(STAGE = case_when(
        ajcc_pathologic_stage == "Stage 0" ~ "0",
        ajcc_pathologic_stage %in% c("Stage I", "Stage IA") ~ "I",
        ajcc_pathologic_stage %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC") ~ "II",
        ajcc_pathologic_stage %in% c("Stage IIIC", "Stage IIIB", "Stage III", "Stage IIIA") ~ "III",
        ajcc_pathologic_stage %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB") ~ "IV",
        TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
      ))
    #Quedarse con las columnas importantes
    clinical_final <- clinical_final %>%
      select(case_id, STAGE)
    #View(clinical_final)
  }
  
  #MUTACIONES SOMÁTICAS
  {
    library(data.table)
    # Define el directorio principal donde están las carpetas
    directorio_principal <- "C:/Users/Denisse González/Downloads/Data/Colorectal"
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
    mutaciones_clean <- mutaciones %>%
      filter(Hugo_Symbol != "Unknown")
    # Solo filas con mutaciones somaticas
 unique(mutaciones$Mutation_Status)
    mutaciones_clean <- mutaciones_clean %>%
      filter(!is.na(Mutation_Status))
    #Añadir columna con identificador
    mutaciones_clean <- mutaciones_clean %>%
      mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
    #Quedarse con las columnas importantes
    mutaciones_clean <- mutaciones_clean %>%
      select(case_id, Tumor_Sample_Barcode, Identificador, Hugo_Symbol, Chromosome, Start_Position, End_Position, Tumor_Seq_Allele1, Tumor_Seq_Allele2 )
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
  data9 <- resultado %>%
    distinct(case_id, Identificador, .keep_all = TRUE)%>%
  rename(SAMPLE_ID = Tumor_Sample_Barcode)%>%
    rename(PATIENT_ID = case_id)
  data9$Identificador <- gsub("chr", "", data9$Identificador)
  View(data9)
    }
#Cambiar a genoma 37
    {
        #BED format
      {
      data_genoma <- data9 %>%
        select (Chromosome, Start_Position, End_Position)
        View(data_genoma)
        #Guardar archivo
        write.table(data_genoma, 
                    file = "data_genoma.bed",  
                    sep = "\t",               # Usar tabulaciones como delimitador
                    row.names = FALSE,        # No incluir números de filas
                    col.names = FALSE,        # No incluir nombres de columnas
                    quote = FALSE)            # No poner comillas alrededor de los valore
      }
        #Archivo con posiciones anteriores (data_genoma)
        {
          data_genoma38 <- data9 %>%
            mutate(
              Start_Position_1 = Start_Position + 1,
              Identificador2 = paste(Chromosome, Start_Position_1, sep = ":"),
              Identificador2 = paste(Identificador2, End_Position, sep = "-") 
            ) %>%
            select(Identificador, Identificador2, Start_Position,Start_Position_1, End_Position, Chromosome, Tumor_Seq_Allele1, Tumor_Seq_Allele2)
          View(data_genoma38)
        }
      #Archivo procesado a GRCh37
      {
        bed_data <- read.table("C:/Users/Denisse González/Downloads/Data/Colorectal/BED.bed", header = FALSE, sep = "\t")
        bed_data <- bed_data[, -ncol(bed_data)]
        colnames(bed_data) <- c("Chromosome2", "Start_Position2", "End_Position2", "Identificador2")
        View(bed_data)
        data_val <- merge(bed_data, data_genoma38, by = "Identificador2", all = TRUE)
        data_val <- distinct(data_val)
        data_val <- data_val %>%
          mutate(
            Identificador_ = paste(
              str_remove(Chromosome2, "^chr"),  # Elimina "chr" del inicio de Chromosome2
              Start_Position2, 
              sep = ":"
            ) %>%
              paste(., Start_Position2, sep = "_") %>%
              paste(., Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "_")
          ) %>%
          select(Identificador, Identificador_)
        View(data_val)
        
      }
    }
    #Merge con el dataset con Identificador anterior para que coincida con los pacientes
    {
    data_9 <- merge(data9, data_val, by = "Identificador", all = FALSE)
      data_9 <- data_9 %>%
      select(SAMPLE_ID, PATIENT_ID, STAGE, Identificador_, Hugo_Symbol)
      data_9 <- data_9 %>%
      rename(Identificador = Identificador_)
    View(data_9)
    }
 }
  #UNIR TODAS
  {
    data_final <- bind_rows(data1, data2, data3, data4, data5, data6, data7, data8, data_9)
    # Verificar si hay filas duplicadas
    duplicates <- data_final %>%
      group_by(SAMPLE_ID, PATIENT_ID, STAGE, Identificador, Hugo_Symbol) %>%
      summarise(count = n(), .groups = 'drop') %>%
      filter(count > 1)
    View(duplicates) 
    # Eliminar las filas duplicadas
    data_final$STAGE[data_final$STAGE == 0] <- "I"
    data_final <- data_final %>%
      distinct(SAMPLE_ID, Identificador, .keep_all = TRUE)
    View(data_final)
  }
  #Clasificación binaria
  {
    #Etapa 1 y 2, 1, las demas 0
    data_final <- data_final %>%
      mutate(classif_stage = ifelse(STAGE %in% c("I", "II"), 1, 0))
    #Etapa 1, 1, las demas 0
    data_final <- data_final %>%
      mutate(classif_stage = ifelse(STAGE %in% c("I"), 1, 0))
    #Etapa2, 1, las demas 0
    #data_final <- data_final %>%
    #  mutate(classif_stage = ifelse(STAGE %in% c("II"), 1, 0))
    #Etapa 3, 1, las demas 0
    #data_final <- data_final %>%
    #  mutate(classif_stage = ifelse(STAGE %in% c("III"), 1, 0))
    #Etapa 4, 1, las demas 0
    #data_final <- data_final %>%
    #  mutate(classif_stage = ifelse(STAGE %in% c("IV"), 1, 0))
    unique(data_final$STAGE)
  }
}

##GRÁFICAS
{
  #Pacientes por etapa
  counts <- data_final %>%
    group_by(STAGE) %>%
    summarise(count = n())
  #View(counts)
  ggplot(counts, aes(x = STAGE, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Pacientes en cada etapa patológica (Adenocarcinoma, NOS)",
         x = "Etapa Patológica AJCC",
         y = "Conteo") +
    theme_minimal() +
    geom_text(aes(label = count), 
              position = position_stack(vjust = 0.5), 
              color = "white")
  
  #Intersecciones entre las mutaciones y las etapas
  library(UpSetR)
  library(tibble)
  mutaciones_stage <- data_final %>%
    select(Identificador, STAGE) %>%
    mutate(value = 1) %>%
    group_by(Identificador, STAGE) %>%
    summarize(count = sum(value), .groups = 'drop') %>%
    pivot_wider(names_from = STAGE, values_from = count, values_fill = list(count = 0))
  # columnas binarias
  mutaciones_binary <- as.data.frame(sapply(mutaciones_stage[-1], function(x) ifelse(x > 0, 1, 0)))
  mutaciones_ready <- cbind(Identificador = mutaciones_stage$Identificador, mutaciones_binary)
  upset(mutaciones_ready, sets = colnames(mutaciones_ready)[-1], order.by = "freq")
  #Num de mutaciones por etapa
  mutaciones_por_etapa <- data_final %>%
    group_by(STAGE) %>%
    summarise(num_mutaciones = n_distinct(Identificador))
  print(mutaciones_por_etapa)
  
  
  
  
  mutaciones_filtradas <- data_final %>%
    group_by(Identificador) %>%
    summarise(num_pacientes = n_distinct(PATIENT_ID), .groups = 'drop') %>%
    filter(num_pacientes >= 5) %>%
    select(Identificador)
  mutaciones_stage <- data_final %>%
    semi_join(mutaciones_filtradas, by = "Identificador") %>%
    select(Identificador, STAGE) %>%
    mutate(value = 1) %>%
    group_by(Identificador, STAGE) %>%
    summarize(count = sum(value), .groups = 'drop') %>%
    pivot_wider(names_from = STAGE, values_from = count, values_fill = list(count = 0))
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

##CLASIFICADOR
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
      select(classif_stage, PATIENT_ID) %>%
      distinct(PATIENT_ID, .keep_all = TRUE) 
    #View(df_modelo)
    # DF con mutaciones unicas
    df_gen_mutado <- filtered_data_final %>%
      select(PATIENT_ID,Identificador) %>%
      mutate(mutado = 1) %>%  # Asignar 1 si la mutación está presente
      distinct() %>%  # no haya duplicados
      pivot_wider(
        names_from = Identificador, 
        values_from = mutado,
        values_fill = list(mutado = 0)  # Rellenar con 0 si la mutación no está presente
      )
    #View(df_gen_mutado)
    # Combinar los dos dataframes
    data_modelo <- merge(df_gen_mutado, df_modelo, by = "PATIENT_ID")
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
    
    x_train <- df_train %>% select(-PATIENT_ID, -classif_stage)
    y_train <- df_train$classif_stage
    x_test <- df_test %>% select(-PATIENT_ID, -classif_stage)
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
    ## RANDOM FOREST
    rf_model <- randomForest(x = x_train, y = as.factor(y_train), ntree = 500)
    print(rf_model)
    # Predicciones
    predictions <- predict(rf_model, x_test)
    # Evaluar el rendimiento
    confusionMatrix(predictions, as.factor(df_test$classif_stage))
    confusion_mat <- confusionMatrix(predictions, as.factor(y_test))
    print(confusion_mat)
    
    
    # Predicciones de probabilidad
    prob_predictions <- predict(rf_model, x_test, type = "prob")
    # Establecer el umbral
    umbral <- 0.4
    y_pred_class <- ifelse(prob_predictions[,2] >= umbral, "1", "0")
    # Convertir las predicciones a factores
    y_pred_class <- factor(y_pred_class, levels = c("0", "1"))
    # Evaluar el rendimiento
    confusionMatrix(y_pred_class, as.factor(y_test))
    
    
    #Gráficas
    confusion_df <- as.data.frame(confusion_mat$table)
    # Gráfico de la matriz de confusión
    ggplot(confusion_df, aes(x = Prediction, y = Reference)) +
      geom_tile(aes(fill = Freq), color = "white") +
      scale_fill_gradient(low = "white", high = "blue") +
      geom_text(aes(label = Freq), color = "black") +
      labs(title = "Matriz de Confusión (datos de prueba)", x = "Predicción", y = "Referencia")

    #Curva ROC y AUC
    
    
    
    # Asegúrate de que y_test esté correctamente formateado como factor y con las clases en el orden correcto
    y_test <- factor(y_test, levels = c(0, 1))  # "0" para etapa tardía, "1" para etapa temprana
    # Obtener las probabilidades de predicción para la clase positiva (por ejemplo, "1" es la clase positiva)
    prob_predictions <- predict(rf_model, x_test, type = "prob")
    # Ahora que las clases están bien definidas, calculamos la curva ROC
    png("curva_roc.png", width = 800, height = 800)
    roc_curve <- roc(y_test, prob_predictions[, 2], plot = TRUE, col = "blue", main = "Curva ROC del Modelo", legacy.axes = TRUE, xlab = "Especificidad", ylab = "Sensibilidad")
    dev.off()
    # Calcular el AUC (Área bajo la curva)
    auc_value <- auc(roc_curve)
    print(paste("AUC:", auc_value))
  
    
    
    
    # Predicciones de probabilidad
    predictions_prob <- predict(rf_model, x_test, type = "prob")[,2]  # Probabilidad de la clase 1
    # Crear la curva ROC
    roc_curve <- roc(y_test, predictions_prob)
    # Calcular el AUC
    auc_value <- auc(roc_curve)
    print(paste("AUC:", auc_value))
    # Dibujar la curva ROC
    png("curva_roc.png", width = 800, height = 800)
    plot(roc_curve, main = "Curva ROC", col = "black", lwd = 2, xlim = c(1, 0), ylim = c(0, 1))  # Ajuste de ejes
    lines(x = c(1, 0), y = c(0, 1), col = "red", lty = 2)  # Línea diagonal de (1,0) a (0,1)
    dev.off()  # Cierra el dispositivo gráfico
  
    
    
    
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
#Mejora de datos para el modelo
{
  
  # Instalar y cargar la librería isotree si no está instalada
  if (!require(isotree)) install.packages("isotree")
  library(isotree)
  
  # Función para detectar anomalías por etapa
  detect_anomalies_by_stage <- function(data, stage_column = "STAGE", num_trees = 100, threshold = 0.5) {
    
    # Lista para guardar los resultados por cada etapa
    anomaly_results <- list()
    
    # Iterar sobre las etapas únicas en la columna STAGE
    for (stage in unique(data[[stage_column]])) {
      
      # Filtrar los datos para la etapa actual (usando la columna STAGE)
      stage_data <- subset(data, data[[stage_column]] == stage)
      
      # Eliminar columnas que no son útiles para el modelo (como case_id, classif_stage, PATIENT_ID, STAGE)
      stage_data_numeric <- stage_data[, !names(stage_data) %in% c("PATIENT_ID", "STAGE", "case_id", "classif_stage")]
      
      # Aplicar Isolation Forest usando isotree
      iso_forest <- isolation.forest(stage_data_numeric, ntree = num_trees)
      
      # Obtener los scores de anomalía para cada paciente
      anomaly_scores <- predict(iso_forest, stage_data_numeric)
      
      # Agregar los puntajes de anomalía al dataframe
      stage_data$anomaly_score <- anomaly_scores  # Los puntajes de anomalía
      stage_data$anomaly <- ifelse(anomaly_scores > threshold, 1, 0)  # 1 si es anómalo, 0 si no
      
      # Guardar los resultados para esa etapa
      anomaly_results[[stage]] <- stage_data
    }
    
    # Unir los resultados de todas las etapas en un único dataframe
    data_with_anomalies <- do.call(rbind, anomaly_results)
    
    # Retornar el dataframe con los resultados de anomalías
    return(data_with_anomalies)
  }
  
  # Aplicar la función a tus datos
  data_with_anomalies <- detect_anomalies_by_stage(data_modelo, "STAGE", num_trees = 100, threshold = 0.5)
  
  # Ver los resultados
  View(data_with_anomalies)
  
  # Ver la distribución de anomalías por etapa
  table(data_with_anomalies$STAGE, data_with_anomalies$anomaly)
  
  
  
  
  
  # Filtrar los datos para que solo queden las instancias normales (anomaly = 0)
  data_cleaned <- data_with_anomalies[data_with_anomalies$anomaly == 0, ]
  View(data_cleaned)
  data_modelo <- data_cleaned %>%
    select(-STAGE, -anomaly_score, -anomaly)
}

##VALIDACION
{
  #TCGA Genomic Data Commons Data Portal
  {
    #Clinico
    {
      clinical <- read.csv("C:/Users/Denisse González/Downloads/Data/Colorectal/clinical.csv", header = TRUE, sep = ",")
      View(clinical)
      # Reemplazar ''--' por NA
      clinical[clinical == "'--"] <- NA
      # Solo filas con stage definida
      clinical_clean <- clinical %>%
        filter(!is.na(ajcc_pathologic_stage))
      #Quedarse con las columnas importantes
      clinical_clean2 <- clinical_clean %>%
        select(case_id, age_at_index, ethnicity, gender, race, ajcc_pathologic_stage, primary_diagnosis, tissue_or_organ_of_origin, classification_of_tumor)
      #View(clinical_clean2)
      # Eliminar duplicados en clinical_clean2
      clinical_clean2_unique <- clinical_clean2 %>%
        distinct(case_id, .keep_all = TRUE)
      #Eliminar metastasis
      clinical_clean2_unique <- clinical_clean2_unique[!grepl("metastasis", clinical_clean2_unique$classification_of_tumor, ignore.case = TRUE), ]
      #View(clinical_clean2_unique)
      #Qudarse con Adenocarcinoma, NOS
      clinical_final <- clinical_clean2_unique %>%
        filter(primary_diagnosis %in% c("Adenocarcinoma, NOS"))
      clinical_final <- clinical_final %>%
        mutate(STAGE = case_when(
          ajcc_pathologic_stage == "Stage 0" ~ "0",
          ajcc_pathologic_stage %in% c("Stage I", "Stage IA") ~ "I",
          ajcc_pathologic_stage %in% c("Stage IIA", "Stage IIB", "Stage II", "Stage IIC") ~ "II",
          ajcc_pathologic_stage %in% c("Stage IIIC", "Stage IIIB", "Stage III", "Stage IIIA") ~ "III",
          ajcc_pathologic_stage %in% c("Stage IVA", "Stage IV", "Stage IVC", "Stage IVB") ~ "IV",
          TRUE ~ "Otro"  # Para manejar cualquier caso no contemplado
        ))
      #Quedarse con las columnas importantes
      clinical_final <- clinical_final %>%
        select(case_id, STAGE)
      View(clinical_final)
    }
    
    #MUTACIONES SOMÁTICAS
    {
      library(data.table)
      # Define el directorio principal donde están las carpetas
      directorio_principal <- "C:/Users/Denisse González/Downloads/Data/Colorectal"
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
      mutaciones_clean <- mutaciones %>%
        filter(Hugo_Symbol != "Unknown")
      # Solo filas con mutaciones somaticas
      unique(mutaciones$Mutation_Status)
      mutaciones_clean <- mutaciones_clean %>%
        filter(!is.na(Mutation_Status))
      #Añadir columna con identificador
      mutaciones_clean <- mutaciones_clean %>%
        mutate(Identificador = paste0(Chromosome, ":", Start_Position, "_", End_Position, "_", Tumor_Seq_Allele1, "_", Tumor_Seq_Allele2))
      #Quedarse con las columnas importantes
      mutaciones_clean <- mutaciones_clean %>%
        select(case_id, Tumor_Sample_Barcode, Identificador, Hugo_Symbol )
      #View(mutaciones_clean)
      #Num de pacientes
      num_unicos2 <- mutaciones_clean %>%
        summarise(count = n_distinct(case_id))
      print(num_unicos2)
    }
    
    #Unir mutaciones con clinico
    resultado <- merge(mutaciones_clean, clinical_final, by = "case_id", all=FALSE)
    #View(resultado)
    #Eliminar duplicados
    data9 <- resultado %>%
      distinct(case_id, Identificador, .keep_all = TRUE)%>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode)%>%
      rename(PATIENT_ID = case_id)
    data9$Identificador <- gsub("chr", "", data9$Identificador)
    View(data9)
    
  }
#Preparar datos de validacion
  {
    # Contar cuántas veces se repite cada mutación
    mutation_counts_val <- table(data9$Identificador)
    # Filtrar las mutaciones que se repiten al menos 5 veces
    mutations_to_keep_val <- names(mutation_counts_val[mutation_counts_val >= 5])
    # Filtrar el conjunto de validación para mantener solo las mutaciones relevantes
    filtered_data9 <- data9[data9$Identificador %in% mutations_to_keep_val, ]
    # DF con mutaciones únicas
    df_gen_mutado_val <- filtered_data9 %>%
      select(PATIENT_ID, Identificador) %>%
      mutate(mutado = 1) %>%  # Asignar 1 si la mutación está presente
      distinct() %>%  # Eliminar duplicados
      pivot_wider(
        names_from = Identificador, 
        values_from = mutado,
        values_fill = list(mutado = 0)  # Rellenar con 0 si la mutación no está presente
      )
    
    # Ver el resultado (opcional)
    # View(df_gen_mutado_val)
    # Extraer las columnas necesarias de data9 para obtener la variable objetivo
    df_modelo_val <- filtered_data9 %>%
      select(PATIENT_ID, STAGE) %>%
      distinct(PATIENT_ID, .keep_all = TRUE)  # Asegurarnos de que solo haya una fila por paciente
    
    # Unir las mutaciones con la variable objetivo
    data_modelo_val <- merge(df_gen_mutado_val, df_modelo_val, by = "PATIENT_ID")
    
    # Ver el resultado final del conjunto de datos
    View(data_modelo_val)
    
  }
  #Predicciones
  {
    # Separar las características (mutaciones) y la variable objetivo (STAGE)
    x_val <- data_modelo_val %>% select(-PATIENT_ID, -STAGE)
    y_val <- data_modelo_val$STAGE
    
    # Realizar las predicciones con el modelo Random Forest entrenado
    predictions_val <- predict(rf_model, newdata = x_val)
    
    # Evaluar el rendimiento con la matriz de confusión
    library(caret)
    confusion_matrix_val <- confusionMatrix(predictions_val, as.factor(y_val))
    print(confusion_matrix_val)
    
    # Calcular la precisión
    accuracy_val <- mean(predictions_val == y_val)
    print(paste("Precisión en el conjunto de validación:", accuracy_val))
    
  }
  
  
  
}
