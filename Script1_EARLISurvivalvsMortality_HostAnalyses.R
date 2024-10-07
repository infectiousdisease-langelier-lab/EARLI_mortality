library(dplyr)
library(tibble)
library(magrittr)
library(limma)
library(hablar)
library(edgeR) # 3.30.3
library(DESeq2) # 1.28.1
library(stringr)
library(tidyverse)
library(reshape2)
library(data.table)
library(ggplot2) # 3.3.4
library(ggpubr) # 0.4.0
library(ggrepel) # 0.9.1
library(ggfortify) # 0.4.14
library(pheatmap) # 1.0.12
library(reshape2) # 1.4.4
library(VennDiagram) # 1.6.20
library(scales) # 1.1.1
library(singscore)
library(SepstratifieR)

####Code Overview - Host analyses for sepsis patients and Non-infectious critical illness
###Steps:
  ###1. Make DE gene tables for both groups (becomes 2A, 2B for Sepsis, 4B and Supp Data for No Sepsis) 
  ###2. Generate Venn diagram of DE genes between four groups (4A)
  ###3. Cytokine analyses of patients (2D, 4C) 
####4. SRSq predictions for patients  based on host gene expression (5A)

#########################
# DEFINE FUNCTIONS
#########################
plot_volcano <- function(res_table, results_folder, prefix_output){
  res_table$sig <- res_table$adj.P.Val < 0.1
  res_table$sig[res_table$adj.P.Val < 0.1 & res_table$logFC>0] <- "deceased"
  res_table$sig[res_table$adj.P.Val < 0.1 & res_table$logFC<0] <- "survivors"
  res_table <- as.data.frame(res_table)
  
  if (sum(res_table$sig!=FALSE)==0){
    label_data <- res_table[res_table$sig, ]
  }else{
    label_data <- res_table[1:min(c(25, nrow(res_table[res_table$sig, ]))), ]
    
    label_data_logfc <- res_table[order(abs(res_table$logFC), decreasing=T), ]
    label_data_logfc <- label_data_logfc[1:25, ][label_data_logfc[1:25, "sig"], ]
    
    label_data <- rbind(label_data, label_data_logfc)
    label_data <- unique(label_data)
  }
  
  pdf(paste0(results_folder, paste0(prefix_output, "_volcano_plot.pdf")), width = 7, height = 7)
  p <- ggplot(res_table, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col = sig)) +
    scale_color_manual(breaks=c("deceased", "survivors", "FALSE"), values = c("goldenrod", "steelblue", "black")) + 
    #geom_text_repel(data = label_data,
    #                aes(label = gene_name),
    #                max.overlaps=20) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 P-value")  +
    theme_bw() +
    theme(legend.position = "none")
  print(p)
  dev.off()
}

###1. Make DE gene tables for patients with sepsis who survived vs those who did not
# Read data
###Gene count data
genecountspc<-read.csv("Inputs/earli_counts_kallisto_mortality.csv", header=TRUE)
###Metadata
metadata<-read.csv("Inputs/Cleaned_Metadata_070324.csv")
##Table of ENSEMBL gene IDs and gene names
gene_attr<-read.csv("Inputs/gene_attr.csv", header=T)

####Process the gene count data
genecountspc<-as.matrix(genecountspc)
rownames(genecountspc)<-genecountspc[,1]
genecountspc<-genecountspc[,-1]
genecountspc <- apply(genecountspc, c(1,2), function(x) { (as.integer(x))})
head(genecountspc)
rank_col <- "logFC"

####remove the 'earli' from the colnames 
colnames(genecountspc)<-substr(colnames(genecountspc), 7,11)

###Select patients with sepsis for primary analysis (this becomes Figure 2A, 2B, 2C) #3438 genes
metadata_analysis<-metadata %>% filter(Group == "1_Sepsis+BldCx+"| Group =="2_Sepsis+OtherCx+")
###scale Age
metadata_analysis$scaleAge <- scale(metadata_analysis$Age)
#Add gender as a factor
metadata_analysis$Gender<-as.factor(metadata_analysis$Gender) 
rownames(metadata_analysis)<-metadata_analysis$Barcode
###Order the gene count data in the same order as the metadata
counts<-genecountspc
x<- c(rownames(metadata_analysis))
y<-c(colnames(counts))
counts<-counts[,colnames(counts) %in% x]
metadata_analysis<-metadata_analysis[rownames(metadata_analysis) %in% y,]
counts<-counts[,order(noquote(colnames(counts)))]
metadata_analysis<-metadata_analysis[order(noquote(rownames(metadata_analysis))),]
identical(rownames(metadata_analysis), colnames(counts))
#Hospital death should be a factor
metadata_analysis$Hospital_Death<-as.factor(metadata_analysis$Hospital_Death)
##Set up target category
target_cat <- "Hospital_Death" 
###Create design model for differential expression
design <- model.matrix(as.formula(paste(paste("~", paste0(" 0 + ", target_cat)), paste0("+ scaleAge + Gender "))), data = metadata_analysis)
dds_paxgene_tmp <- DGEList(counts = counts)
#normalization factor, including library death and number of genes expressed. Linear model fit  
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
rownames(gene_attr) <- gene_attr$gene_id
res_paxgene$gene_name <- gene_attr[rownames(res_paxgene), "gene_name"]
# replace NA values with 1s and keep only significant genes
res_paxgene$adj.P.Val[is.na(res_paxgene$adj.P.Val)] <- 1
# filter
sig_res_paxgene <- data.frame(res_paxgene)
sig<-sig_res_paxgene %>% filter (adj.P.Val <0.1)
sigup<-sig %>% filter (logFC>0)
sigdown<-sig %>% filter (logFC<0)
dim(sig)
write.csv(sig_res_paxgene, "Outputs/Supplementary_Data1_Sepsis_DEGenes.csv")

####Figure 2A: Create volcano plot for sepsis DE Genes 
plot_volcano(res_paxgene, "Outputs", "/Sepsis")

###Figure 2B. Selection of specific genes of interest to plot - becomes Figure 2B 
specific<-counts 
###Pick out genes of interest including CD3D(ENSG00000167286), BPI (ENSG00000101425), MMP8 (ENSG00000118113)
specific<-counts[c("ENSG00000167286", "ENSG00000101425", "ENSG00000118113"),]
flip<-as.data.frame(t(specific))
flip$Barcode<-NA
for (i in (1:nrow(flip))){flip$Barcode[i]<-rownames(flip)[i]}
outcome<-metadata_analysis[,c("Barcode", "Hospital_Death")]
genespecific<-merge(outcome, flip, by="Barcode")
head(genespecific) 
write.csv(genespecific, "Outputs/Source_Data1_SpecificGenesSepsis.csv")

#Figure 2C: made by WebGestalt of Supplementary Table 1 DE Genes



####Figure 4A, B - Repeat DE analysis for ICU Controls with age, sex as covariates 
metadata_analysis<-metadata %>% filter(Group == "4_NO_Sepsis") 
###scale Age
metadata_analysis$scaleAge <- scale(metadata_analysis$Age)
#Add gender as a factor
metadata_analysis$Gender<-as.factor(metadata_analysis$Gender) 

rownames(metadata_analysis)<-metadata_analysis$Barcode
###Order the gene count data in the same order as the metadata
counts<-genecountspc
x<- c(rownames(metadata_analysis))
y<-c(colnames(counts))
counts<-counts[,colnames(counts) %in% x]
metadata_analysis<-metadata_analysis[rownames(metadata_analysis) %in% y,]
counts<-counts[,order(noquote(colnames(counts)))]
metadata_analysis<-metadata_analysis[order(noquote(rownames(metadata_analysis))),]
identical(rownames(metadata_analysis), colnames(counts))
#Hospital death should be a factor
metadata_analysis$Hospital_Death<-as.factor(metadata_analysis$Hospital_Death)
##Set up target category
target_cat <- "Hospital_Death" 
###Create design model for differential expression
design <- model.matrix(as.formula(paste(paste("~", paste0(" 0 + ", target_cat)), paste0("+ scaleAge + Gender "))), data = metadata_analysis)
dds_paxgene_tmp <- DGEList(counts = counts)
#normalization factor, including library death and number of genes expressed. Linear model fit  
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
rownames(gene_attr) <- gene_attr$gene_id
res_paxgene$gene_name <- gene_attr[rownames(res_paxgene), "gene_name"]
# replace NA values with 1s and keep only significant genes
res_paxgene$adj.P.Val[is.na(res_paxgene$adj.P.Val)] <- 1
# filter
sig_res_paxgene <- data.frame(res_paxgene)
sig<-sig_res_paxgene %>% filter (adj.P.Val <0.1)
sigup<-sig %>% filter (logFC>0)
sigdown<-sig %>% filter (logFC<0)
dim(sig)
write.csv(sig_res_paxgene, "Outputs/Supplementary_Data8_ICUControls_DEGenes.csv")



####Supplemental Figure 2 - Repeat DE analysis for ICU Controls with age and sex as covariates, adding cardiac arrest as a factor 
metadata_analysis<-metadata %>% filter(Group == "4_NO_Sepsis") 
#Create a metadata column for cardiac arrest
metadata_analysis$cardiacarrest<-0

for (i in (1:nrow(metadata_analysis))) {
  if (metadata_analysis$Primary.Diagnosis[i] == "Cardiovascular:Cardiac arrest")
  {metadata_analysis$cardiacarrest[i] = 1}
}
###scale Age
metadata_analysis$scaleAge <- scale(metadata_analysis$Age)
#Add gender as a factor
metadata_analysis$Gender<-as.factor(metadata_analysis$Gender) 
#add cardiac arrest as a factor 
metadata_analysis$cardiacarrest<-as.factor(metadata_analysis$cardiacarrest) 
rownames(metadata_analysis)<-metadata_analysis$Barcode
###Order the gene count data in the same order as the metadata
counts<-genecountspc
x<- c(rownames(metadata_analysis))
y<-c(colnames(counts))
counts<-counts[,colnames(counts) %in% x]
metadata_analysis<-metadata_analysis[rownames(metadata_analysis) %in% y,]
counts<-counts[,order(noquote(colnames(counts)))]
metadata_analysis<-metadata_analysis[order(noquote(rownames(metadata_analysis))),]
identical(rownames(metadata_analysis), colnames(counts))
#Hospital death should be a factor
metadata_analysis$Hospital_Death<-as.factor(metadata_analysis$Hospital_Death)
##Set up target category
target_cat <- "Hospital_Death" 
###Create design model for differential expression
design <- model.matrix(as.formula(paste(paste("~", paste0(" 0 + ", target_cat)), paste0("+ scaleAge + Gender + cardiacarrest"))), data = metadata_analysis)
dds_paxgene_tmp <- DGEList(counts = counts)
#normalization factor, including library death and number of genes expressed. Linear model fit  
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
rownames(gene_attr) <- gene_attr$gene_id
res_paxgene$gene_name <- gene_attr[rownames(res_paxgene), "gene_name"]
# replace NA values with 1s and keep only significant genes
res_paxgene$adj.P.Val[is.na(res_paxgene$adj.P.Val)] <- 1
# filter
sig_res_paxgene <- data.frame(res_paxgene)
sig<-sig_res_paxgene %>% filter (adj.P.Val <0.1)
sigup<-sig %>% filter (logFC>0)
sigdown<-sig %>% filter (logFC<0)
dim(sig)
write.csv(sig_res_paxgene, "Outputs/Supplementary_Data10_ICUControls_DEGenes_CardiacArrestasCovariate.csv")



########4A. Venn diagram creation - this becomes Figure 4B
sepsispax<-read.csv("Outputs/Supplementary_Data1_Sepsis_DEGenes.csv")
icucontrolspax<-read.csv("Outputs/Supplementary_Data8_ICUControls_DEGenes.csv")

#Pull out the significant genes at FDR <0.1
sepsispax_sig<-sepsispax %>% filter(adj.P.Val<0.1)
icucontrolspax_sig<-icucontrolspax %>% filter(adj.P.Val<0.1)

##create columns that say if genes are significant in either group or in both
sepsispax_sig$overlap_with_icucontrols<-sepsispax_sig$gene_name %in% icucontrolspax_sig$gene_name
sepsispax_sig$significant_in_sepsis<-TRUE
sepsispax_sig_short<-sepsispax_sig[,c(8,10,9)]
icucontrolspax_sig$significant_in_icucontrols<-TRUE
icucontrolspax_sig<-icucontrolspax_sig[,c(8,9)]

merge<-merge(sepsispax_sig_short, icucontrolspax_sig, by="gene_name", all.x=TRUE, all.y=TRUE)
for (i in (1:nrow(merge))) {
  if (is.na(merge$significant_in_sepsis[i] == TRUE))
  {merge$significant_in_sepsis[i] = FALSE}
}

for (i in (1:nrow(merge))) {
  if (is.na(merge$overlap_with_icucontrols[i] == TRUE))
  {merge$overlap_with_icucontrols[i] = FALSE}
}

for (i in (1:nrow(merge))) {
  if (is.na(merge$significant_in_icucontrols[i] == TRUE))
  {merge$significant_in_icucontrols[i] = FALSE}
}

table(merge$significant_in_icucontrols)
table(merge$significant_in_sepsis)
table(merge$overlap_with_icucontrols)

write.csv(merge, "Outputs/Source_Data7_DEGene_Overlap.csv")




######Figures 2D,Cytokine analyses
cytomort<-metadata[,c("Barcode", "Hospital_Death", "Group", "PAI1", "IL6", "sTNFr1", "RAGE", "IL8", "ANG2", "Proc_C", "ICAM")]
#cconvert all to numeric
cytomort$PAI1<-as.numeric(cytomort$PAI1)
cytomort$IL6<-as.numeric(cytomort$IL6)
cytomort$sTNFr1<-as.numeric(cytomort$sTNFr1)
cytomort$RAGE<-as.numeric(cytomort$RAGE)
cytomort$IL8<-as.numeric(cytomort$IL8)
cytomort$ANG2<-as.numeric(cytomort$ANG2)
cytomort$Proc_C<-as.numeric(cytomort$Proc_C)
cytomort$ICAM<-as.numeric(cytomort$ICAM)
cytomort_sepsis <-cytomort %>% filter(Group == "1_Sepsis+BldCx+"|Group == "2_Sepsis+OtherCx+")
cytomort_icucontrols <-cytomort %>% filter(Group == "4_NO_Sepsis")
dim(cytomort)
write.csv(cytomort, "Outputs/Source_Data2_CytokineandOutcomeData.csv")

#####Pvalue calculations by group - SEPSIS
###Input data for Table 4
sepsis_deceased<-cytomort_sepsis %>% filter(Hospital_Death == 1)
sepsis_survived<-cytomort_sepsis %>% filter(Hospital_Death == 0)
quantile(sepsis_survived$PAI1, na.rm=TRUE)
quantile(sepsis_survived$IL6, na.rm=TRUE)
quantile(sepsis_survived$sTNFr1, na.rm=TRUE)
quantile(sepsis_survived$RAGE, na.rm=TRUE)
quantile(sepsis_survived$IL8, na.rm=TRUE)
quantile(sepsis_survived$ANG2, na.rm=TRUE)
quantile(sepsis_survived$Proc_C, na.rm=TRUE)
quantile(sepsis_survived$ICAM, na.rm=TRUE)

quantile(sepsis_deceased$PAI1, na.rm=TRUE)
quantile(sepsis_deceased$IL6, na.rm=TRUE)
quantile(sepsis_deceased$sTNFr1, na.rm=TRUE)
quantile(sepsis_deceased$RAGE, na.rm=TRUE)
quantile(sepsis_deceased$IL8, na.rm=TRUE)
quantile(sepsis_deceased$ANG2, na.rm=TRUE)
quantile(sepsis_deceased$Proc_C, na.rm=TRUE)
quantile(sepsis_deceased$ICAM, na.rm=TRUE)

##Fold changes 
fold_PAI1<-(mean(sepsis_deceased$PAI1,na.rm=TRUE))/(mean(sepsis_survived$PAI1, na.rm=TRUE))
fold_IL6<-(mean(sepsis_deceased$IL6,na.rm=TRUE))/(mean(sepsis_survived$IL6, na.rm=TRUE))
fold_sTNFr1<-(mean(sepsis_deceased$sTNFr1,na.rm=TRUE))/(mean(sepsis_survived$sTNFr1, na.rm=TRUE))
fold_RAGE<-(mean(sepsis_deceased$RAGE,na.rm=TRUE))/(mean(sepsis_survived$RAGE, na.rm=TRUE))
fold_IL8<-(mean(sepsis_deceased$IL8,na.rm=TRUE))/(mean(sepsis_survived$IL8, na.rm=TRUE))
fold_ANG2<-(mean(sepsis_deceased$ANG2,na.rm=TRUE))/(mean(sepsis_survived$ANG2, na.rm=TRUE))
fold_Proc_C<-(mean(sepsis_deceased$Proc_C,na.rm=TRUE))/(mean(sepsis_survived$Proc_C, na.rm=TRUE))
fold_ICAM<-(mean(sepsis_deceased$ICAM,na.rm=TRUE))/(mean(sepsis_survived$ICAM, na.rm=TRUE))
fold_PAI1
fold_IL6
fold_sTNFr1
fold_RAGE
fold_IL8
fold_ANG2
fold_Proc_C
fold_ICAM

##P values corrected for multiple testing 
p_vals <- c()
for (col_ in colnames(cytomort_sepsis)){
  if (!(col_ %in% c("Barcode", "HospitalDeath", "Group"))){
    cytomort_sepsis_tmp <- cytomort_sepsis
    cytomort_sepsis_tmp$cell_t <-   cytomort_sepsis_tmp[, col_]
    p_vals[col_] <- wilcox.test(cell_t~Hospital_Death,   cytomort_sepsis_tmp)$p.value
  }
}
p_adj <- p.adjust(p_vals,method="BH")
p_adj


###ICU Controls
#####Pvalue calculations by group 
###Input data for Table 4
icucontrol_deceased<-cytomort_icucontrols %>% filter(Hospital_Death == 1)
icucontrol_survived<-cytomort_icucontrols %>% filter(Hospital_Death == 0)
quantile(icucontrol_survived$PAI1, na.rm=TRUE)
quantile(icucontrol_survived$IL6, na.rm=TRUE)
quantile(icucontrol_survived$sTNFr1, na.rm=TRUE)
quantile(icucontrol_survived$RAGE, na.rm=TRUE)
quantile(icucontrol_survived$IL8, na.rm=TRUE)
quantile(icucontrol_survived$ANG2, na.rm=TRUE)
quantile(icucontrol_survived$Proc_C, na.rm=TRUE)
quantile(icucontrol_survived$ICAM, na.rm=TRUE)

quantile(icucontrol_deceased$PAI1, na.rm=TRUE)
quantile(icucontrol_deceased$IL6, na.rm=TRUE)
quantile(icucontrol_deceased$sTNFr1, na.rm=TRUE)
quantile(icucontrol_deceased$RAGE, na.rm=TRUE)
quantile(icucontrol_deceased$IL8, na.rm=TRUE)
quantile(icucontrol_deceased$ANG2, na.rm=TRUE)
quantile(icucontrol_deceased$Proc_C, na.rm=TRUE)
quantile(icucontrol_deceased$ICAM, na.rm=TRUE)

##Fold changes 
fold_PAI1<-(mean(icucontrol_deceased$PAI1,na.rm=TRUE))/(mean(icucontrol_survived$PAI1, na.rm=TRUE))
fold_IL6<-(mean(icucontrol_deceased$IL6,na.rm=TRUE))/(mean(icucontrol_survived$IL6, na.rm=TRUE))
fold_sTNFr1<-(mean(icucontrol_deceased$sTNFr1,na.rm=TRUE))/(mean(icucontrol_survived$sTNFr1, na.rm=TRUE))
fold_RAGE<-(mean(icucontrol_deceased$RAGE,na.rm=TRUE))/(mean(icucontrol_survived$RAGE, na.rm=TRUE))
fold_IL8<-(mean(icucontrol_deceased$IL8,na.rm=TRUE))/(mean(icucontrol_survived$IL8, na.rm=TRUE))
fold_ANG2<-(mean(icucontrol_deceased$ANG2,na.rm=TRUE))/(mean(icucontrol_survived$ANG2, na.rm=TRUE))
fold_Proc_C<-(mean(icucontrol_deceased$Proc_C,na.rm=TRUE))/(mean(icucontrol_survived$Proc_C, na.rm=TRUE))
fold_ICAM<-(mean(icucontrol_deceased$ICAM,na.rm=TRUE))/(mean(icucontrol_survived$ICAM, na.rm=TRUE))
fold_PAI1
fold_IL6
fold_sTNFr1
fold_RAGE
fold_IL8
fold_ANG2
fold_Proc_C
fold_ICAM

##P values corrected for multiple testing 
p_vals <- c()
for (col_ in colnames(cytomort_icucontrols)){
  if (!(col_ %in% c("Barcode", "HospitalDeath", "Group"))){
    cytomort_icucontrols_tmp <- cytomort_icucontrols
    cytomort_icucontrols_tmp$cell_t <-   cytomort_icucontrols_tmp[, col_]
    p_vals[col_] <- wilcox.test(cell_t~Hospital_Death,   cytomort_icucontrols_tmp)$p.value
  }
}
p_adj <- p.adjust(p_vals,method="BH")
p_adj




#########################
# SRS PREDICTIONS - becomes part of Figure 5
#########################
library("SepstratifieR")

counts_unfmt <- read.csv("Inputs/earli_counts_kallisto_mortality.csv", header=TRUE, row.names = 1)

paxgene_cnt_logcpm <- edgeR::cpm(t(counts_unfmt), log=T)

#genecountspc<-read.csv("Inputs/earli_counts_kallisto_mortality.csv", header=TRUE)
srs_preds <- stratifyPatients(paxgene_cnt_logcpm, gene_set = "extended")

write.csv(srs_preds@SRSq, "Outputs/Sepsis_SRSq.csv")

#########################
# SWEENEY PAPER GENES
#########################
gene_attr<-read.csv("Inputs/gene_attr.csv", header=T, row.names = 1)

dds_score <- DGEList(counts = counts_unfmt)
voom_score <- voom(dds_score)

# external bact signature
logcpm_score_ranks <- rankGenes(voom_score$E)
death1_genes_up <- rownames(gene_attr[gene_attr$gene_name %in% c("TRIB1", "CKS2", "MKI67", "POLD3", "PLK1"), ])
death1_genes_dn <- rownames(gene_attr[gene_attr$gene_name %in% c("TGFBI", "LY86", "CST3", "CBFA2T3", "RCBTB2", 
                                                                          "TST", "CX3CR1", "CD5", "MTMR11", "CLEC10A", 
                                                                          "EMR3", "DHRS7B", "CEACAM8"), ])

death1_scores <- simpleScore(logcpm_score_ranks, upSet = death1_genes_up, downSet = death1_genes_dn)$TotalScore

death2_genes_up <- rownames(gene_attr[gene_attr$gene_name %in% c("CFD", "DDIT4", "DEFA4", "IFI27", "IL1R2", 
                                                                          "IL8", "MAFF", "OCLN", "RGS1"), ])
death2_genes_dn <- rownames(gene_attr[gene_attr$gene_name %in% c("AIM2", "APH1A", "CCR2", "EIF5A", "GSTM1", 
                                                                          "HIST1H3H", "NT5E", "RAB40B", "VNN3"), ])

death2_scores <- simpleScore(logcpm_score_ranks, upSet = death2_genes_up, downSet = death2_genes_dn)$TotalScore

death3_genes_up <- rownames(gene_attr[gene_attr$gene_name %in% c("B4GALT4", "BPI", "CD24", "CEP55", "CTSG",
                                                                          "DDIT4", "G0S2", "MPO", "MT1G", "NDUFV2", 
                                                                          "PAM", "PSMA6", "SEPP1"), ])
death3_genes_dn <- rownames(gene_attr[gene_attr$gene_name %in% c("ABCB4", "CTSS", "IKZF2", "NT5E"), ])

death3_scores <- simpleScore(logcpm_score_ranks, upSet = death3_genes_up, downSet = death3_genes_dn)$TotalScore

death4_genes_up <- rownames(gene_attr[gene_attr$gene_name %in% c("DEFA4", "CD163", "PER1", "RGS1", "HIF1A",
                                                                          "SEPP1", "C11orf74", "CIT"), ])
death4_genes_dn <- rownames(gene_attr[gene_attr$gene_name %in% c("LY86", "TST", "OR52R1", "KCNJ2"), ])

death4_scores <- simpleScore(logcpm_score_ranks, upSet = death4_genes_up, downSet = death4_genes_dn)$TotalScore

death_scores_df <- as.data.frame(list(death1_scores, death2_scores, death3_scores, death4_scores))
rownames(death_scores_df) <- colnames(logcpm_score_ranks)
colnames(death_scores_df) <- c("death1_scores", "death2_scores", "death3_scores", "death4_scores")

write.csv(death_scores_df, "Outputs/Sepsis_death_scores.csv")

