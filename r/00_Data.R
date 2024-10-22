rm(list = ls())

library(GenomicFeatures)
library(caret)

tx_db <- makeTxDbFromGFF("data/gencode.v46.chr_patch_hapl_scaff.annotation.gtf", format = "gtf")
exons_list_per_gene <- exonsBy(tx_db, by = "gene")

gene_lengths <- data.frame(
  ID = sub("\\..*", "", names(exons_list_per_gene)), 
  Length = sapply(exons_list_per_gene, function(exons) {
    sum(width(reduce(exons)))
  })
)


meta_data <- read.csv(file = "data/meta_purity.csv", row.names = 1)
dim(meta_data) # 564  39
meta_data <- meta_data[
  meta_data$LibraryType == "Stranded_Total" & 
    meta_data$Local_Recurrence.. != "FALSE" & 
    as.logical(meta_data$Meets.Purity.Threshold) & 
    (meta_data$Non.Surgical.Treatment == "Radiotherapy and TMZ" | meta_data$Non.Surgical.Treatment == "Radiotherapy and TMZ +") & 
    meta_data$IDH == 0, 
]
dim(meta_data) # 208  39
meta_data$Group <- ifelse(meta_data$Cohort == "EORTC", "Validation", "Discovery")
sum(meta_data$Group == "Discovery") # 92

meta_validation <- meta_data[
  meta_data$Group == "Validation", 
]
nrow(meta_validation)  # 116
colData_validation <- data.frame(ID = row.names(meta_validation), 
                                 Stage = meta_validation[, "Stage"],
                                 Source = meta_validation[, "Sample.Source"])

meta_discovery <- meta_data[
  meta_data$Group == "Discovery", 
]
nrow(meta_discovery)  # 100
colData_discovery <- data.frame(ID = row.names(meta_discovery), 
                                Stage = meta_discovery[, "Stage"], 
                                Source = meta_discovery[, "Sample.Source"])


raw_data <- read.csv(file = "data/raw_data.csv", row.names = 1)
dim(raw_data) # 608 67718
raw_data <- raw_data[, colnames(raw_data) %in% gene_lengths$ID]
dim(raw_data) # 608 62903
raw_data <- raw_data[rownames(meta_data), , drop = FALSE]
dim(raw_data) # 208 62903

countData_validation <- t(raw_data[meta_data$Group == "Validation", ])
dim(countData_validation) # 62903   116
countData_validation <- na.omit(countData_validation) # "ENSG00000130283"
dim(countData_validation) # 62902   116

countData_discovery <- t(raw_data[meta_data$Group == "Discovery", ])
dim(countData_discovery) # 62903   92
countData_discovery <- na.omit(countData_discovery)  # "ENSG00000130283"
dim(countData_discovery) # 62902   92

label_validation <- data.frame(sample = row.names(meta_validation), 
                               label = meta_validation[, "Stage"],
                               batch = meta_validation[, "Sample.Source"])
label_validation$label <- ifelse(label_validation$label == "Primary", 0, 1)
label_validation$subject = factor(gsub("_P|_R", "", colData_validation$ID))

label_discovery <- data.frame(sample = row.names(meta_discovery), 
                              label = meta_discovery[, "Stage"],
                              batch = meta_discovery[, "Sample.Source"])
label_discovery$label <- ifelse(label_discovery$label == "Primary", 0, 1)
label_discovery$subject = factor(gsub("_P|_R", "", colData_discovery$ID))

set.seed(17)

unique_sources <- unique(meta_data$Sample.Source)
train_patients <- c()

for (source in unique_sources) {
  source_data <- meta_data[meta_data$Sample.Source == source, ]
  unique_patients <- unique(source_data$Patient)
  train_patients <- c(train_patients, sample(unique_patients, size = floor(0.78 * length(unique_patients))))
}

meta_train <- meta_data[meta_data$Patient %in% train_patients, ]
meta_test <- meta_data[!meta_data$Patient %in% train_patients, ]

colData_train <- data.frame(ID = row.names(meta_train), 
                            Stage = meta_train$Stage, 
                            Source = meta_train$Sample.Source)
colData_test <- data.frame(ID = row.names(meta_test), 
                           Stage = meta_test$Stage, 
                           Source = meta_test$Sample.Source)

countData_train <- t(raw_data[rownames(meta_train), ])
countData_train <- na.omit(countData_train)
countData_test <- t(raw_data[rownames(meta_test), ])
countData_test <- na.omit(countData_test)

label_train <- data.frame(sample = row.names(meta_train), 
                          label = ifelse(meta_train$Stage == "Primary", 0, 1), 
                          batch = meta_train$Sample.Source,
                          subject = factor(gsub("_P|_R", "", row.names(meta_train))))
label_test <- data.frame(sample = row.names(meta_test), 
                         label = ifelse(meta_test$Stage == "Primary", 0, 1), 
                         batch = meta_test$Sample.Source,
                         subject = factor(gsub("_P|_R", "", row.names(meta_test))))

save(colData_test, colData_train, 
     countData_test, countData_train, 
     meta_test, meta_train, gene_lengths, 
     file = "step0output.Rdata")

write.csv(label_test, 'data/label_test.csv', row.names = FALSE)
write.csv(label_train, 'data/label_train.csv', row.names = FALSE)
