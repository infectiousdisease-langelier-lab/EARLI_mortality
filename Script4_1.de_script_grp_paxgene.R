######################################################
#
######################################################

rm(list = ls())

# load libraries
library(limma)
library(edgeR) # 3.30.3
library(DESeq2) # 1.28.1
library(stringr)
library(tidyverse)
library(reshape2)
library(data.table)

# set comp
comp_ <- c("G1", "G2")

group_dict <- c('G1'='1_Sepsis+BldCx+',
                'G2'='2_Sepsis+OtherCx+',
                'G4'='4_NO_Sepsis')

suffix_ <- paste(comp_, collapse='') 

vars_to_adjust <- c()
#vars_to_adjust <- c("apacheiii_scaled", "microbial_mass")
if (length(vars_to_adjust)>0){
  vars_to_adjust[1] <- paste0("+ ", vars_to_adjust[1])
}

#########################
# SET PARAMETERS
#########################
# target variable
target_cat <- "Hospital_Death"

suffix_ <- paste(suffix_, target_cat, sep="_")

# set paths
data_path <- "Inputs/"
results_path <-
  paste0("Outputs/")

#########################
# DATA LOADING & PRE-PROCESSING
#########################
# counts data
paxgene_cnt_data <- read.csv(paste0(data_path, "earli_counts_kallisto_mortality.csv"), row.names=1)

# meta data
paxgene_meta_data <- read.csv(paste0(data_path, "Cleaned_Metadata_070324.csv"))
rownames(paxgene_meta_data) <- paste0("EARLI_", paxgene_meta_data$Barcode)

# samples in common
paxgene_samples <- intersect(rownames(paxgene_meta_data),
                             colnames(paxgene_cnt_data))

# count data
paxgene_cnt_data <- paxgene_cnt_data[, paxgene_samples]

# meta data
paxgene_meta_data <- paxgene_meta_data[paxgene_samples, ]

# format count data
# drop duplicated ENSG identifiers
paxgene_cnt_data <- paxgene_cnt_data[!grepl(pattern = ".*(\\.).*(\\.).*", rownames(paxgene_cnt_data)), ]
rownames(paxgene_cnt_data) <- sapply(rownames(paxgene_cnt_data), function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])

# format meta data
paxgene_meta_data$EARLI_Barcode <- rownames(paxgene_meta_data)
paxgene_meta_data$Gender[paxgene_meta_data$Barcode == "50413"] <- "Male"

paxgene_meta_data$Gender[paxgene_meta_data$Gender=="Female"] <- 2
paxgene_meta_data$Gender[paxgene_meta_data$Gender!=2] <- 1

write.csv(paxgene_cnt_data, 
          paste0(results_path, paste0("paxgene_cnt_data.csv")))
write.csv(paxgene_meta_data, 
          paste0(results_path, paste0("paxgene_meta_data.csv")))

# filter for groups of interest
paxgene_meta_data <- paxgene_meta_data[paxgene_meta_data$Group %in% group_dict[comp_], ]
paxgene_meta_data <- paxgene_meta_data[sort(rownames(paxgene_meta_data)), ]

# rm patient without plasma data
paxgene_meta_data <- paxgene_meta_data[rownames(paxgene_meta_data) != "EARLI_11942" , ]

paxgene_cnt_data <- paxgene_cnt_data[, rownames(paxgene_meta_data)]

#########################
# DATA ANALYSIS
#########################
for (i in c(1:10)){
  # read test labels
  test_probs_file <- read.csv(paste0(results_path, paste0("test_labels_", paste0(paste(i, suffix_, sep="_"), ".csv"))), 
                              row.names = 1)
  test_samples_to_exclude <- test_probs_file$X0
  
  # filter for train only
  paxgene_cnt_data_filt <- paxgene_cnt_data[, !(colnames(paxgene_cnt_data) %in% test_samples_to_exclude)]
  
  # filter the meta data to match the selected samples
  # make sure our variable of interest is a factor
  paxgene_meta_data_filt <- paxgene_meta_data[!(rownames(paxgene_meta_data) %in% test_samples_to_exclude), ]

  paxgene_meta_data_filt$Group <- as.factor(paxgene_meta_data_filt$Group)
  
  paxgene_meta_data_filt[, target_cat] <- as.factor(paxgene_meta_data_filt[, target_cat])
  
  # standard scaling on Age
  paxgene_meta_data_filt$age_scaled <- scale(paxgene_meta_data_filt$Age)
  
  # make sure gender is a factor
  paxgene_meta_data_filt$gender <- as.factor(paxgene_meta_data_filt$Gender)
  
  # save filtered counts and metadata
  # make sure test samples included
  write.csv(paxgene_cnt_data[rownames(paxgene_cnt_data_filt), !(colnames(paxgene_cnt_data) %in% test_samples_to_exclude)], 
            paste0(data_path, paste0("paxgene_cnt_data_", paste(i, suffix_, sep="_")), ".csv"))
  write.csv(paxgene_meta_data[!(rownames(paxgene_meta_data) %in% test_samples_to_exclude), ], 
            paste0(data_path, paste0("paxgene_meta_data_", paste(i, suffix_, sep="_")), ".csv"))
  
  # DE analysis
  paxgene_meta_data_filt[, target_cat] <- relevel(paxgene_meta_data_filt[, target_cat], ref = "0")

  dds_paxgene_tmp <- DGEList(counts = paxgene_cnt_data_filt)
  design <- model.matrix(as.formula(paste(paste("~", paste0(" 0 + ", target_cat)), paste0("+ age_scaled + gender ", paste(vars_to_adjust, collapse = " + ")))), 
                         data = paxgene_meta_data_filt)
  
  dds_paxgene_tmp <- calcNormFactors(dds_paxgene_tmp)
  dds_paxgene_tmp <- voom(dds_paxgene_tmp,design, plot=T)

  dds_paxgene_tmp <- lmFit(dds_paxgene_tmp, design)
  
  cont_mat <- makeContrasts(paste(paste0(target_cat, 1), paste0(target_cat, 0), sep=" - "), 
                            levels=design)

  dds_paxgene_tmp <- contrasts.fit(dds_paxgene_tmp, cont_mat)
  dds_paxgene_tmp <- eBayes(dds_paxgene_tmp, robust=TRUE)
  res_paxgene <- topTable(dds_paxgene_tmp, sort.by = "P", n = Inf)
  
  # sort the genes from lowest to highest given adjusted p-values
  res_paxgene <- res_paxgene[order(res_paxgene$adj.P.Val, decreasing = F), ]
  
  # replace NA values with 1s and keep only significant genes
  res_paxgene$adj.P.Val[is.na(res_paxgene$adj.P.Val)] <- 1
  
  # filter
  sig_res_paxgene <- data.frame(res_paxgene)
  
  # save the output as a CSV file
  if (length(vars_to_adjust)>0){
    write.csv(sig_res_paxgene, 
              paste0(results_path, paste0("sig_res_paxgene_", paste0(paste(i, paste(suffix_, paste(str_replace(vars_to_adjust, "\\+ ", ""), collapse = "_AND_"), sep="_"), sep="_"), ".csv"))))
    
  }else{
    write.csv(sig_res_paxgene, 
              paste0(results_path, paste0("sig_res_paxgene_", paste0(paste(i, suffix_, sep="_"), ".csv"))))
  }
}

