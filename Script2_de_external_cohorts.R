######################################################
# Script2_de_external_cohorts.R
# created on Feb 23, 2024
# lucile.neyton@ucsf.edu
######################################################
rm(list = ls())

# load libraries
library(ggplot2)
library(edgeR)
library(ggrepel)
library(stringr)
library(scales)
library(reshape2)
library(ggpubr)
library(biomaRt)

# set ggplot2 theme
theme_update(text = element_text(size=12),
             axis.text.x = element_text(size=12),
             axis.text.y = element_text(size=12),
             axis.title.x = element_text(size=12),
             axis.title.y = element_text(size=12),
             legend.text = element_text(size=12),
             legend.title = element_text(size=12),
             strip.text.x = element_text(size=12),
             strip.text.y = element_text(size=12),
             plot.title = element_text(size=12))

# set paths
gains_data_path <- "Inputs/GAINS/"

results_path <- "Outputs/"

######################### GAINS #########################
#########################
# DATA LOADING
#########################
count_data_filt <- read.csv(paste0(gains_data_path, "geneLevel_data.csv"), row.names = 1, check.names = F)
meta_data <- read.csv(paste0(gains_data_path, "class_labels.csv"), row.names = 1)
gene_data <- read.csv(paste0(gains_data_path, "../gene_attr.csv"), row.names = 1)

#########################
# DATA PRE-PROCESSING
#########################
# get gene symbols from biomart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
entrez_map <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
      filters = "ensembl_gene_id",
      values = colnames(count_data_filt), 
      mart = ensembl)
entrez_map <- entrez_map[!duplicated(entrez_map$ensembl_gene_id), ]
rownames(entrez_map) <- entrez_map$ensembl_gene_id

# filter genes
count_data_filt <- t(count_data_filt)
count_data_filt <- count_data_filt[rownames(count_data_filt) %in% gene_data$gene_id, ]

# filter meta data
meta_data <- meta_data[colnames(count_data_filt), ]

# add death as binary variable
meta_data$deceased <- (meta_data$X28d_mortality=="non-survivor")

#########################
# PREPARE META DATA
#########################
# change to numeric
meta_data$age_scaled <- scale(as.numeric(meta_data$age))

# change to factor
meta_data$deceased <- as.factor(meta_data$deceased)

# relevel factor levels
meta_data$deceased <- relevel(meta_data$deceased, ref = "FALSE")

#########################
# DIFFERENCES BETWEEN SURVIVORS AND NON-SURVIVORS
#########################
dge_counts <- DGEList(counts = count_data_filt)
design <- model.matrix(~0 + deceased + age_scaled + sex, data = meta_data)

fit <- lmFit(dge_counts$counts, design)

cont_mat <- makeContrasts("deceasedTRUE - deceasedFALSE", levels=design)
fit.de <- contrasts.fit(fit, cont_mat)
fit.de <- eBayes(fit.de, robust=TRUE, trend=T)

top_table <- topTable(fit.de, sort.by = "P", n = Inf)

# add gene symbols
top_table$gene_name <- entrez_map[rownames(top_table), "hgnc_symbol"]

# volcano plot
write.csv(top_table, paste0(results_path, "Supplementary_Data6_GAINS_DEGenes.csv"))








