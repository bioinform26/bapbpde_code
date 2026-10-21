rm(list=ls())
setwd("/raw_data")

library(GEOquery)
library(tidyverse)
library(affy)
library(sva)
library(hgu133plus2.db)
library(data.table)
library(WGCNA)

#### GSE30219####
# read CEL files
FileName <- list.files(path = "GSE30219_RAW/")
sample_id <- str_sub(FileName, 1, 10)# Adjust the GSM ID length if needed.
Ids1 <- data.frame(FileName = FileName,
                   sample_id = sample_id,
                   stringsAsFactors = F)
Ids1<-Ids1[!duplicated(Ids1$sample_id),]

Ids=read.csv('GSE30219_group_clinic.csv',sep=',',header = T)
rownames(Ids)<-Ids[,1]
Ids<-Ids[Ids1$sample_id,]

Ids$FileName<-Ids1$FileName

GSE30219_cel <- ReadAffy(filenames = Ids$FileName,
                         celfile.path = "GSE30219_RAW",
                         phenoData = Ids) 
GSE30219_rma <- affy::rma(GSE30219_cel)

#correct batches and annotate probes.
Ids$batch <- str_sub(GSE30219_cel@protocolData@data$ScanDate, 1, 8)
mod <- model.matrix(~as.factor(tissue_type), data = Ids)
GSE30219_combat1 <- ComBat(dat = exprs(GSE30219_rma), batch = Ids$batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
anno_gpl570 <- toTable(hgu133plus2SYMBOL) %>%
  set_names("probe_id", "symbol")
source('Function/AnnoTate.R')
expr_anno <- AnnoTate(GSE30219_combat1, anno_gpl570, fun = mean)
#Adenocarcinoma
expr_anno <- as.data.frame(t(expr_anno))
stopifnot(all(Ids$sample_id %in% rownames(expr_anno)))
expr_anno <- expr_anno[Ids$sample_id, , drop = FALSE]
expr_anno <- expr_anno[
  which(Ids$tissue_type == "Adenocarcinoma"),
  , drop = FALSE
]
#zscore
expr_anno <- as.data.frame(scale(expr_anno))
clinic <- read.csv('GSE30219_clinical.csv')
stopifnot(!anyDuplicated(rownames(expr_anno)), !anyDuplicated(clinic$ID),
          all(clinic$ID %in% rownames(expr_anno)))
expr_anno <- expr_anno[clinic$ID, , drop = FALSE]
data <- cbind(clinic[, 1:3], expr_anno)
# Preserve the original month-to-day conversion. Verify the input unit.
data$OS.time <- data$OS.time * 30
write.csv(data, 'GSE30219_OS_exp.csv', row.names = FALSE)

#### GSE31210####
#read CEL files
FileName <- list.files(path = "GSE31210_RAW/")
sample_id <- str_sub(FileName, 1, 9)# Adjust the GSM ID length if needed.
Ids1 <- data.frame(FileName = FileName,
                   sample_id = sample_id,
                   stringsAsFactors = F)
Ids1<-Ids1[!duplicated(Ids1$sample_id),]

Ids=read.csv('GSE31210_group_clinic.csv',sep=',',header = T)
rownames(Ids)<-Ids[,1]
Ids<-Ids[Ids1$sample_id,]
Ids$FileName<-Ids1$FileName

GSE31210_cel <- ReadAffy(filenames = Ids$FileName,
                         celfile.path = "GSE31210_RAW",
                         phenoData = Ids) 
GSE31210_rma <- affy::rma(GSE31210_cel)
#correct batches and annotate probes.
Ids$batch <- str_sub(GSE31210_cel@protocolData@data$ScanDate, 1, 8)
mod <- model.matrix(~as.factor(tissue_type), data = Ids)
GSE31210_combat1 <- ComBat(dat = exprs(GSE31210_rma), batch = Ids$batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
anno_gpl570 <- toTable(hgu133plus2SYMBOL) %>%
  set_names("probe_id", "symbol")
expr_anno <- AnnoTate(GSE31210_combat1, anno_gpl570, fun = mean)
#Adenocarcinoma
expr_anno <- as.data.frame(t(expr_anno))
stopifnot(all(Ids$sample_id %in% rownames(expr_anno)))
expr_anno <- expr_anno[Ids$sample_id, , drop = FALSE]
expr_anno <- expr_anno[
  which(Ids$tissue_type == "Adenocarcinoma"),
  , drop = FALSE
]
#zscore
expr_anno <- as.data.frame(scale(expr_anno))

clinic <- read.csv('GSE31210_clinical.csv')
stopifnot(!anyDuplicated(rownames(expr_anno)), !anyDuplicated(clinic$ID),
          all(clinic$ID %in% rownames(expr_anno)))
expr_anno <- expr_anno[clinic$ID, , drop = FALSE]
data <- cbind(clinic[, 1:3], expr_anno)
write.csv(data, 'GSE31210_OS_exp.csv', row.names = FALSE)

#### GSE72094####
#: read CEL files
FileName <- list.files(path = "GSE72094_RAW/")
sample_id <- str_sub(FileName, 1, 10)# Adjust the GSM ID length if needed.
Ids1 <- data.frame(FileName = FileName,
                   sample_id = sample_id,
                   stringsAsFactors = F)
Ids1<-Ids1[!duplicated(Ids1$sample_id),]

Ids=read.csv('GSE72094_group_clinic.csv',sep=',',header = T)
rownames(Ids)<-Ids[,1]
Ids<-Ids[Ids1$sample_id,]
Ids$FileName<-Ids1$FileName

# One-time CDF setup; uncomment only if the package is not installed.
# library(makecdfenv)
# make.cdf.package("GPL15048_HuRSTA_2a520709.CDF.gz", "hursta2a520709cdf",
#                  species = "Homo sapiens", compress = TRUE, unlink = TRUE)
# install.packages("hursta2a520709cdf", repos = NULL, type = "source")
library(hursta2a520709cdf)
GSE72094_cel <- ReadAffy(filenames = Ids$FileName,
                         celfile.path = "GSE72094_RAW",
                         phenoData = Ids,
                         cdfname = "hursta2a520709cdf") 
GSE72094_rma <- affy::rma(GSE72094_cel)
#correct batches and annotate probes.
Ids$batch <- str_sub(GSE72094_cel@protocolData@data$ScanDate, 1, 8)
mod <- model.matrix(~1, data = Ids)
GSE72094_combat1 <- ComBat(dat = exprs(GSE72094_rma), batch = Ids$batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
gpl <- getGEO("GPL15048", destdir = ".")
anno_GPL15048 <- Table(gpl)[, c(1, 4)]
anno_GPL15048$probe_id <- anno_GPL15048$ID
anno_GPL15048$symbol <- anno_GPL15048$GeneSymbol
clean <- anno_GPL15048 %>% filter(symbol != "")
expr_anno <- AnnoTate(GSE72094_combat1, clean, fun = mean)
#Adenocarcinoma
expr_anno <- as.data.frame(t(expr_anno))
stopifnot(all(Ids$sample_id %in% rownames(expr_anno)))
expr_anno <- expr_anno[Ids$sample_id, , drop = FALSE]
expr_anno <- expr_anno[
  which(Ids$tissue_type == "Adenocarcinoma"),
  , drop = FALSE
]
#zscore
expr_anno <- as.data.frame(scale(expr_anno))
clinic <- read.csv('GSE72094_clinical.csv')
stopifnot(!anyDuplicated(rownames(expr_anno)), !anyDuplicated(clinic$ID),
          all(clinic$ID %in% rownames(expr_anno)))
expr_anno <- expr_anno[clinic$ID, , drop = FALSE]
data <- cbind(clinic[, 1:3], expr_anno)
write.csv(data, 'GSE72094_OS_exp.csv', row.names = FALSE)

#### TCGA####
load('TCGA.Rdata')
stopifnot(nrow(tcga_df) > 0, ncol(tcga_df) > 1)

# Check gene identifiers.
gene_id <- trimws(as.character(tcga_df[[1]]))

stopifnot(
  !anyNA(gene_id),
  all(nzchar(gene_id)),
  !anyDuplicated(gene_id)
)

# Check expression values before conversion.
expr_df <- tcga_df[, -1, drop = FALSE]

stopifnot(all(vapply(expr_df, is.numeric, logical(1))))

tcga <- as.matrix(expr_df)
storage.mode(tcga) <- "double"
rownames(tcga) <- gene_id

stopifnot(
  all(is.finite(tcga)),
  all(tcga >= 0)
)

# Standardize and validate TCGA sample barcodes.
sample_id <- gsub(".", "-", colnames(tcga), fixed = TRUE)

stopifnot(
  !anyDuplicated(sample_id),
  all(grepl(
    "^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}-[0-9]{2}[A-Za-z]($|-)",
    sample_id
  ))
)

# Positions 14–15 encode sample type; position 16 encodes vial.
sample_type <- substr(sample_id, 14, 15)
patient_id <- substr(sample_id, 1, 12)

# Retain primary tumors.
keep_tumor <- sample_type == "01"

stopifnot(sum(keep_tumor) > 1)

sample_info <- data.frame(
  sample_id = sample_id[keep_tumor],
  patient_id = patient_id[keep_tumor],
  sample_type = sample_type[keep_tumor],
  stringsAsFactors = FALSE
)

tcga <- tcga[, keep_tumor, drop = FALSE]

# Do not silently choose among multiple samples from one patient.
duplicated_patient <- unique(
  sample_info$patient_id[
    duplicated(sample_info$patient_id) |
      duplicated(sample_info$patient_id, fromLast = TRUE)
  ]
)

if (length(duplicated_patient) > 0) {
  stop(
    "Multiple primary-tumor samples found for these patients. ",
    "Resolve them using a predefined sample-quality rule: ",
    paste(duplicated_patient, collapse = ", ")
  )
}

colnames(tcga) <- sample_info$patient_id

# Log transformation; no automatic cluster-based sample removal.
tcga <- log2(tcga + 1)

# Load survival information.
clinical <- read_table(file.path(
  clinical_dir,
  "TCGA_clinical.csv"
))

stopifnot(all(c("ID", "OS.time", "OS") %in% names(clinical)))

clinical <- as.data.frame(
  clinical[, c("ID", "OS.time", "OS")]
)
clinical$ID <- trimws(as.character(clinical$ID))

# Require patient-level TCGA IDs and valid clinical values.
stopifnot(
  !anyNA(clinical),
  all(grepl(
    "^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}$",
    clinical$ID
  )),
  !anyDuplicated(clinical$ID),
  is.numeric(clinical$OS.time),
  all(is.finite(clinical$OS.time)),
  all(clinical$OS.time >= 0),
  all(clinical$OS %in% c(0, 1))
)

# OS.time must have been verified to use the intended time unit.

missing_id <- setdiff(clinical$ID, colnames(tcga))

if (length(missing_id) > 0) {
  stop(
    "Clinical patients without a primary-tumor expression sample: ",
    paste(missing_id, collapse = ", ")
  )
}

# Match expression columns to clinical rows explicitly.
tcga <- tcga[, match(clinical$ID, colnames(tcga)), drop = FALSE]

stopifnot(
  ncol(tcga) > 1,
  identical(colnames(tcga), clinical$ID),
  all(is.finite(tcga))
)
sample_info <- sample_info[
  match(clinical$ID, sample_info$patient_id),
  ,
  drop = FALSE
]

# Avoid duplicate column names after binding.
stopifnot(!any(rownames(tcga) %in% names(clinical)))


# Gene-wise Z-score across matched tumor patients.
gene_mean <- rowMeans(tcga)
gene_sd <- apply(tcga, 1, sd)

# Remove genes with zero or invalid SD to avoid NaN/Inf.
keep_gene <- is.finite(gene_sd) & gene_sd > 0

tcga <- tcga[keep_gene, , drop = FALSE]
gene_mean <- gene_mean[keep_gene]
gene_sd <- gene_sd[keep_gene]

stopifnot(nrow(tcga) > 0)

tcga <- sweep(tcga, 1, gene_mean, "-")
tcga <- sweep(tcga, 1, gene_sd, "/")

stopifnot(all(is.finite(tcga)))

# Bind clinical information and standardized expression.
result <- cbind(
  clinical,
  as.data.frame(t(tcga), check.names = FALSE)
)

write.csv(result,"TCGA_OS_exp.csv",row.names = FALSE)



