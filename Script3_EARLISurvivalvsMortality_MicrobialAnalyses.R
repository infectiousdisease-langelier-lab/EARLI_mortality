####Script 3 - Microbial analyses for sepsis patient and Non-infectious critical illness
###Steps:
###For  Sepsis 
###1. Compare microbial mass between survival and mortality 
###2. Compare microbial dominance between survival and mortality 
###3. Compare beta diversity between survival and mortality

###For  ICU Controls 
###4. Compare microbial mass between survival and mortality 
###5. Compare microbial dominance between survival and mortality
###6. Compare beta diversity between survival and mortality

##For Sepsis Only:
###7. Generate per-pathogen presence/absence and outcome data for all RBM-proved bacterial
#pathogens and also DNA viruses and calculate a fisher exact test to examine survival vs mortality for each pathogen 
###8 Compare gene expression data between sepsis patients with high or low microbial dominance metrics 

library(phyloseq)
library(mia)
library(dplyr)
library(tidyr)
library(data.table)
library(ggpubr)
library(vegan)
library(rstatix)
library(ggplot2)
library(edgeR)
library(ggrepel)
library(reshape2)

results_path <- "Outputs/"

######Read in files
metadata<-read.csv("Inputs/Cleaned_Metadata_070324.csv")
# background-corrected data 
bgc_data <- read.csv("Inputs/bgccounts_mortality.csv")
# taxonomic data from NCBI
tax_data <- readRDS("Inputs/tax_data.RDS")

#########################
# DATA PRE-PROCESSING
#########################
metadata$EARLI_Barcode <- paste0("EARLI_", metadata$Barcode)
rownames(metadata) <- metadata$EARLI_Barcode
metadata_full <- metadata

# format for phyloseq
ranklist_ <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon")
tax_data <- rbindlist(tax_data, idcol = T) %>%
  group_by(.id) %>%
  distinct(rank, .keep_all = T) %>%
  dplyr::select(-id) %>%
  pivot_wider(values_from = "name", names_from = "rank") %>%
  dplyr::rename(taxon = .id) %>%
  dplyr::select(all_of(ranklist_)) %>%
  ungroup

# format and melt
bgc_data$EARLI_Barcode <- sapply(bgc_data$sample_name, function(x) paste0("EARLI_", strsplit(x, "_", fixed=T)[[1]][2]))
bgc_data_full <- bgc_data

# taxonomic data
tax_data <- as.data.frame(tax_data)
rownames(tax_data) <- tax_data$taxon

tax_data_full <- tax_data

for (suffix_ in c("Sepsis", "NonInfectiousCriticalIllness")){
  
  if (suffix_ == "Sepsis"){
    metadata<-metadata_full %>% filter(Group == "1_Sepsis+BldCx+"| Group =="2_Sepsis+OtherCx+") 
    name_out <- "Supplementary_Data4_sepsis_microbialmass_and_pathdominance.csv"
  }else{
    if (suffix_ == "NonInfectiousCriticalIllness"){
      metadata<-metadata_full %>% filter(Group == "4_NO_Sepsis")
      name_out <- "Supplementary_Data12_ICUControls_microbialmass_and_pathdominance.csv"
    }
  }
  
  # included samples
  bgc_data <- bgc_data_full[bgc_data_full$EARLI_Barcode %in% metadata$EARLI_Barcode, ]
  
  # remove duplicate sample
  #bgc_data <- bgc_data[bgc_data$sample_name!="EARLI_11255_PLaD1_L6_OPS046_S92", ]
  bgc_data <- bgc_data[bgc_data$sample_name!="EARLI_11255_PLaD1_A10_OPS046_S145", ]
  
  bgc_data_cast <- dcast(bgc_data, tax_id~EARLI_Barcode, value.var="nt_rpm")
  bgc_data_cast[is.na(bgc_data_cast)] <- 0
  rownames(bgc_data_cast) <- bgc_data_cast$tax_id
  bgc_data_cast <- bgc_data_cast[, colnames(bgc_data_cast)!="tax_id"]
  
  tax_data <- tax_data_full[rownames(bgc_data_cast), ]
  
  # filter meta data
  bgc_meta_data_filt <- metadata[colnames(bgc_data_cast), ]
  bgc_meta_data_filt$Hospital_Death <- factor(bgc_meta_data_filt$Hospital_Death, levels=c("0", "1"))
  
  # create phyloseq object
  # need OTU, taxonomy, and metadata
  bgc_data_cast <- otu_table(as.matrix(bgc_data_cast), taxa_are_rows = TRUE)
  tax_data <- tax_table(as.matrix(tax_data))
  meta_data_phyloseq <- sample_data(bgc_meta_data_filt)
  microb_data <- phyloseq(bgc_data_cast, tax_data, meta_data_phyloseq)
  
  # aggregate at genus level
  microb_data_genus <- tax_glom(microb_data, taxrank = "genus")
  microb_data_family <- tax_glom(microb_data, taxrank = "family")
  microb_data_glom <- list(genus=microb_data_genus, family=microb_data_family)
  
  #########################
  # DATA ANALYSIS
  ###1. Compare microbial mass between survival and mortality - Figures 3A, 4C, Supp Figure 3
  ###2. Compare microbial dominance between survival and mortality   for sepsis and ICU controls - Figures 3B, 4D,  Supp Figure 3
  ###3. Compare beta diversity between survival and mortality - Figure 3C
  #########################
  # rank_names(microb_data)
  tax_data %>%
    as.data.frame %>%
    filter(superkingdom == "Bacteria") %>% 
    rownames() -> bacterialist
  
  microb_data <- list(bacteria=bacterialist)
  
  # for bacteria at the genus level
  microb_type <- "bacteria"
  microb_list <- microb_data[[microb_type]]
  
  rank_ <- "genus"
  print(rank_)
  
  microb_data_rank <- microb_data_glom[[rank_]]
  
  # data viz
  bact_data_rank <- prune_taxa(microb_list, microb_data_rank)
  
  # dominance
  dom_df <- estimateDominance(makeTreeSummarizedExperimentFromPhyloseq(bact_data_rank), 
                              abund_values = "counts", 
                              index="relative")
  
  print(wilcox.test(colData(dom_df)$relative ~ sample_data(bact_data_rank)$Hospital_Death))
  median(colData(dom_df)$relative[sample_data(bact_data_rank)$Hospital_Death=="1"])
  median(colData(dom_df)$relative[sample_data(bact_data_rank)$Hospital_Death=="0"])
  
  # alternative calculation - same results
  dom_alt_df <- otu_table(bact_data_rank)
  dom_alt_vect <- apply(dom_alt_df, 2, max)/colSums(dom_alt_df)
  
  dom_alt_LCA_df <- data.frame(dom_alt_vect, factor(sample_data(bact_data_rank)$Hospital_Death, levels=c("0", "1")))
  colnames(dom_alt_LCA_df) <- c("dominant", "Hospital_Death")
  
  # beta diversity -> compare diversities bw samples
  # add centroids to plot
  ord_df <- ordinate(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank), method = "NMDS", distance = "bray")
  scrs <- scores(ord_df)
  
  centroids <- setNames(aggregate(scrs$sites, by=list(bgc_meta_data_filt[rownames(scrs$sites), "Hospital_Death"]), FUN=mean), c("Hospital_Death", "NMDS1", "NMDS2"))
  
  rownames(centroids) <- centroids$Hospital_Death
  centroids$Hospital_Death <- factor(centroids$Hospital_Death, levels=c("0", "1"))
  
  line_df <- data.frame(scores(ord_df)$sites)
  line_df$Hospital_Death <- bgc_meta_data_filt[rownames(scrs$sites),]$Hospital_Death
  
  line_df$NMDS1c <- centroids[line_df$Hospital_Death, "NMDS1"]
  line_df$NMDS2c <- centroids[line_df$Hospital_Death, "NMDS2"]
  
  # permanova
  bact_data_gen_beta <- phyloseq::distance(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank), 
                                           method="bray")
  
  # sig
  set.seed(123)
  permanova_res <- adonis2(formula=bact_data_gen_beta ~ as(sample_data(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank)), "data.frame")$Hospital_Death)
  print(permanova_res$`Pr(>F)`[1])
  
  #mod2 <- betadisper(bact_data_gen_beta, as(sample_data(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank)), "data.frame")$Hospital_Death) ## messages
  #print(permutest(mod2))
  
  p_nmds <- plot_ordination(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank), 
                            ord_df, color="Hospital_Death") +
    annotate("text", size=3, label=paste0("PERMANOVA P=", permanova_res$`Pr(>F)`[1]), x=max(scrs$sites[, 1])/3, y=max(scrs$sites[, 2])) +
    scale_color_manual(values=c("0"="steelblue", "1"="goldenrod")) + 
    theme_bw() + theme(text = element_text(family = "Helvetica", size=8),
                       axis.text.x = element_text(family = "Helvetica", size=8),
                       axis.text.y = element_text(family = "Helvetica", size=8),
                       axis.title.x = element_text(family = "Helvetica", size=8),
                       axis.title.y = element_text(family = "Helvetica", size=8),
                       legend.text = element_text(family = "Helvetica", size=8),
                       legend.title = element_text(family = "Helvetica", size=8),
                       strip.text.x = element_text(family = "Helvetica", size=8),
                       strip.text.y = element_text(family = "Helvetica", size=8),
                       plot.title = element_text(family = "Helvetica", size=8))
  
  p_nmds <- p_nmds + xlim(min(scrs$sites[,1]),max(scrs$sites[,1])) + ylim(min(scrs$sites[,2]),max(scrs$sites[,2]))
  
  p_nmds <- p_nmds + geom_segment(data=line_df, 
                                  aes(x=NMDS1c, y=NMDS2c, xend=NMDS1, yend=NMDS2, colour=Hospital_Death), 
                                  size=0.25, show.legend=FALSE)
  
  p_nmds <- p_nmds + geom_point(data=centroids, size=5, pch=21, color="black")
  
  pdf(paste0(results_path, paste0(paste(microb_type, paste(suffix_, rank_, sep="_"), sep="_"), "_nmds_boxplot.pdf")), height = 3, width = 4)
  print(p_nmds)
  dev.off()
  
  # family level
  microb_list <- microb_data[["bacteria"]]
  microb_data_rank <- microb_data_glom[["family"]]
  bact_data_rank <- prune_taxa(microb_list, microb_data_rank)
  
  bact_family_data <- otu_table(bact_data_rank)
  
  tax_data_df <- data.frame(tax_data)
  row_ind <- rownames(bact_family_data) %in% tax_data_df[tax_data_df$family=="Enterobacteriaceae" & !is.na(tax_data_df$family), "taxon"]
  enterobacteriaceae_prop <- bact_family_data[row_ind, ] / (colSums(bact_family_data[!(row_ind), ]) + bact_family_data[row_ind, ])
  sample_data(bact_data_rank)$enterobacteriaceae_prop <- t(enterobacteriaceae_prop)
  
  sample_data(bact_data_rank)$rel_dom <- dom_alt_LCA_df$dominant
  
  df_out <- data.frame(sample_data(bact_data_rank)[, c("Hospital_Death", "microbial_mass", "rel_dom", "enterobacteriaceae_prop")])
  
  write.csv(df_out, paste0(results_path, name_out))
}

###Add cardiac arrest to the No ICU output file
output<-read.csv("Outputs/Supplementary_Data17_ICUControls_microbialmass_and_pathdominance.csv")
output$Barcode<-substr(output$X, 7,11)
#Create a metadata column for cardiac arrest
metadata$cardiacarrest<-0
for (i in (1:nrow(metadata))) {
  if (metadata$Primary.Diagnosis[i] == "Cardiovascular:Cardiac arrest")
  {metadata$cardiacarrest[i] = 1}
}
metadata_short<-metadata %>% filter (Group == "4_NO_Sepsis")
metadata_short<-metadata_short[,c(1,38)]
output<-merge(metadata_short, output, by="Barcode", all.x=TRUE, all.y=FALSE)
names(output)
write.csv(output, "Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv")

###CREATING A DATABASE OF NT_RPM PER ESTABLISHED BACTERIAL PATHOGEN PER patient with sepsis - becomes Figures 3E and Supp Figure 1
######Read in files
metadata<-read.csv("Inputs/Cleaned_Metadata_070324.csv")
sepsis<-metadata %>% filter(Group == "1_Sepsis+BldCx+"| Group =="2_Sepsis+OtherCx+") 
rbm1<-sepsis[,c(1,32,33)]

rbm2<-sepsis[,c(1,34,35)]
colnames(rbm2)[2]<-"name"
colnames(rbm2)[3]<-"nt_rpm"

rbm3<-sepsis[,c(1,36,37)]
colnames(rbm3)[2]<-"name"
colnames(rbm3)[3]<-"nt_rpm"

rbm<-rbind(rbm1,rbm2, rbm3)

#Create an empty data frame 
q<- unique(sepsis$Barcode)
res = data.frame(Barcode = q)

#Pickout each pathogen and list their nt_rpm 
Bacteriodes <-rbm %>% dplyr::filter(rbm$name == "Bacteroides fragilis" | rbm$name == "Bacteroides thetaiotaomicron" |
                                      rbm$name == "Bacteroides vulgatus")
Bacteriodes_s<-Bacteriodes[c("Barcode", "nt_rpm")]
colnames(Bacteriodes_s)[2]<-"Bacteriodes"
res<-merge(res,Bacteriodes_s, all.x=TRUE, by = "Barcode")

#Citrobacter
Citrobacter <-rbm %>% dplyr::filter(rbm$name == "Citrobacter portucalensis" | rbm$name == "Citrobacter braakii")
Citrobacter_s<-Citrobacter[c("Barcode", "nt_rpm")]
colnames(Citrobacter_s)[2]<-"Citrobacter"
res<-merge(res,Citrobacter_s, all.x=TRUE)

#Ecoli
Escherichia1 <-rbm %>% dplyr::filter(rbm$name == "Escherichia coli")
Escherichia_s<-Escherichia1[c("Barcode", "nt_rpm")]
colnames(Escherichia_s)[2]<-"Escherichia"
res<-merge(res,Escherichia_s, all.x=TRUE)

#Enterobacter
Enterobacter <-rbm %>% dplyr::filter(rbm$name == "Enterobacter hormaechei" | rbm$name == "Enterobacter roggenkampii")
Enterobacter_s<-Enterobacter[c("Barcode", "nt_rpm")]
colnames(Enterobacter_s)[2]<-"Enterobacter"
res<-merge(res,Enterobacter_s, all.x=TRUE)

#Klebsiella
Klebsiella <-rbm %>% dplyr::filter(rbm$name == "Klebsiella michiganensis" | 
                                     rbm$name == "Klebsiella quasivariicola"
                                   |rbm$name == "Klebsiella pneumoniae" | 
                                     rbm$name == "Klebsiella oxytoca")
Klebsiella_s<-Klebsiella[c("Barcode", "nt_rpm")]
colnames(Klebsiella_s)[2]<-"Klebsiella"
res<-merge(res,Klebsiella_s, all.x=TRUE)

#Haemophilus
Haemophilus <-rbm %>% dplyr::filter(rbm$name == "Haemophilus influenzae")
Haemophilus_s<-Haemophilus[c("Barcode", "nt_rpm")]
colnames(Haemophilus_s)[2]<-"Haemophilus"
res<-merge(res,Haemophilus_s, all.x=TRUE)

#Lactobacillus
Lactobacillus <-rbm %>% dplyr::filter(rbm$name == "Lactobacillus crispatus" | 
                                        rbm$name == "Lactobacillus johnsonii"
                                      |rbm$name == "Lactobacillus fermentum")
Lactobacillus_s<-Lactobacillus[c("Barcode", "nt_rpm")]
colnames(Lactobacillus_s)[2]<-"Lactobacillus"
res<-merge(res,Lactobacillus_s, all.x=TRUE)

#Morganella
Morganella <-rbm %>% dplyr::filter(rbm$name == "Morganella morganii")
Morganella_s<-Morganella[c("Barcode", "nt_rpm")]
colnames(Morganella_s)[2]<-"Morganella"
res<-merge(res,Morganella_s, all.x=TRUE)

#PRevotella
Prevotella <-rbm %>% dplyr::filter(rbm$name == "Prevotella intermedia")
Prevotella_s<-Haemophilus[c("Barcode", "nt_rpm")]
colnames(Prevotella_s)[2]<-"Prevotella"
res<-merge(res, Prevotella_s, all.x=TRUE)

#Proteus
Proteus <-rbm %>% dplyr::filter(rbm$name == "Proteus mirabilis")
Proteus_s<-Proteus[c("Barcode", "nt_rpm")]
colnames(Proteus_s)[2]<-"Proteus"
res<-merge(res, Proteus_s, all.x=TRUE)

#Pseudomonas
Pseudomonas <-rbm %>% dplyr::filter(rbm$name == "Pseudomonas aeruginosa")
Pseudomonas_s<-Pseudomonas[c("Barcode", "nt_rpm")]
colnames(Pseudomonas_s)[2]<-"Pseudomonas"
res<-merge(res, Pseudomonas_s, all.x=TRUE)

#Serratia 
Serratia <-rbm %>% dplyr::filter(rbm$name == "Serratia rubidaea" |rbm$name == "Serratia marcescens")
Serratia_s<-Serratia[c("Barcode", "nt_rpm")]
colnames(Serratia_s)[2]<-"Serratia"
res<-merge(res, Serratia_s, all.x=TRUE)

#Staphylococcus 
Staphylococcus <-rbm %>% dplyr::filter(rbm$name == "Staphylococcus argenteus" |rbm$name == "Staphylococcus aureus")
Staphylococcus_s<-Staphylococcus[c("Barcode", "nt_rpm")]
colnames(Staphylococcus_s)[2]<-"Staphylococcus"
res<-merge(res, Staphylococcus_s, all.x=TRUE)

#Streptococcus 
Streptococcus <-rbm %>% dplyr::filter(rbm$name == "Streptococcus agalactiae" |rbm$name == "Streptococcus salivarius"|rbm$name == "Streptococcus mitis" |
                                        rbm$name == "Streptococcus oralis"|rbm$name == "Streptococcus pyogenes"|rbm$name == "Streptococcus pneumoniae"
                                      |rbm$name == "Streptococcus thermophilus"|rbm$name == "Streptococcus dysgalactiae")
Streptococcus_s<-Streptococcus[c("Barcode", "nt_rpm")]
colnames(Streptococcus_s)[2]<-"Streptococcus"
res<-merge(res, Streptococcus_s, all.x=TRUE)

#### now we have a matrix of nt_rpm per pathogen vs patient barcode 
bacteria <- merge(rbm, res, all.x=TRUE, all.y=TRUE)

###remove double rows
bacteria<-bacteria[!duplicated(bacteria$Barcode),]

###CREATING A DATABASE OF NT_RPM PER DNA VIRAL PATHOGEN PER PATIENT 
bgc_data1 <- read.csv("Inputs/bgccounts_mortality.csv")
virus<-bgc_data1

#Restrict to only samples with viruses 
virus<-virus[virus$category == "viruses",]

###restrict to viruses present at >0.1 nt_rpm
virus<-virus[virus$nt_rpm > 0.1,]
#Pull out the barcode_name
for (i in 1:nrow(virus)){
  virus$barcode_name<-substr(virus$sample_name, 7,11)
}

###use grepl to combine MULTIPLE different naems per virsu (in particular for annelovirus family)
virus$name[grepl("Torque|Anello|TTV|TTV-like|SEN virus|uncultured", virus$name)] <- "Anelloviridae"
virus$name[grepl("Human betaherpesvirus 6|Human betaherpesvirus 6B", virus$name)] <- "HHV-6"
virus$name[grepl("Human polyomavirus 1|Human polyomavirus 2", virus$name)] <- "Human polyomavirus"
virus$name[grepl("Human gammaherpesvirus 4", virus$name)] <- "HHV-4"
virus$name[grepl("Human betaherpesvirus 5", virus$name)] <- "HHV-5"
virus$name[grepl("Hepatitis B virus", virus$name)] <- "Hepatitis B"
virus$name[grepl("Human gammaherpesvirus 8", virus$name)] <- "HHV-8"
virus$name[grepl("Human alphaherpesvirus 1", virus$name)] <- "HSV-1"
virus$name[grepl("Human alphaherpesvirus 2", virus$name)] <- "HSV-2"
virus$name[grepl("Human alphaherpesvirus 3", virus$name)] <- "VZV"
table(virus$name)
q<-unique(virus$barcode_name)
res = data.frame(Barcode = q)

###for anellovirus only, summing the nt_rpm that match to each member of the anneloviridae
Anelloviridae <-virus %>% dplyr::filter(virus$name == "Anelloviridae")
Anelloviridae<-Anelloviridae[c("barcode_name", "nt_rpm")]
Anello<- Anelloviridae %>%  group_by(barcode_name) %>% 
  summarise_all(sum)
colnames(Anello)[2]<-"Anelloviridae"
res<-merge(res,Anello, by.x="Barcode", by.y="barcode_name", all.x=TRUE)


#HBV
HepatitisB <-virus %>% dplyr::filter(virus$name == "Hepatitis B")
HepatitisB_s<-HepatitisB[c("barcode_name", "nt_rpm")]
colnames(HepatitisB_s)[2]<-"HepatitisB"
res<-merge(res,HepatitisB_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#HSV1
HSV1 <-virus %>% dplyr::filter(virus$name == "HSV-1")
HSV1_s<-HSV1[c("barcode_name", "nt_rpm")]
colnames(HSV1_s)[2]<-"HSV1"
res<-merge(res,HSV1_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#HSV2
HSV2 <-virus %>% dplyr::filter(virus$name == "HSV-2")
HSV2_s<-HSV2[c("barcode_name", "nt_rpm")]
colnames(HSV2_s)[2]<-"HSV2"
res<-merge(res,HSV2_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#VZV
VZV <-virus %>% dplyr::filter(virus$name == "VZV")
VZV_s<-HSV2[c("barcode_name", "nt_rpm")]
colnames(VZV_s)[2]<-"VZV"
res<-merge(res,VZV_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#HHV4
HHV4 <-virus %>% dplyr::filter(virus$name == "HHV-4")
HHV4_s<-HHV4[c("barcode_name", "nt_rpm")]
colnames(HHV4_s)[2]<-"HHV4"
res<-merge(res,HHV4_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#HHV5
HHV5 <-virus %>% dplyr::filter(virus$name == "HHV-5")
HHV5_s<-HHV5[c("barcode_name", "nt_rpm")]
colnames(HHV5_s)[2]<-"HHV5"
res<-merge(res,HHV5_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#HHV6
HHV6 <-virus %>% dplyr::filter(virus$name == "HHV-6")
HHV6_s<-HHV6[c("barcode_name", "nt_rpm")]
colnames(HHV6_s)[2]<-"HHV6"
res<-merge(res,HHV6_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#HHV8
HHV8 <-virus %>% dplyr::filter(virus$name == "HHV-8")
HHV8_s<-HHV8[c("barcode_name", "nt_rpm")]
colnames(HHV8_s)[2]<-"HHV8"
res<-merge(res,HHV8_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#Boca virus
Bocavirus <-virus %>% dplyr::filter(virus$name == "Human bocavirus")
Bocavirus_s<-Bocavirus[c("barcode_name", "nt_rpm")]
colnames(Bocavirus_s)[2]<-"Bocavirus"
res<-merge(res,Bocavirus_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

#Polyomavirus
Polyomavirus <-virus %>% dplyr::filter(virus$name == "Human polyomavirus")
Polyomavirus_s<-Polyomavirus[c("barcode_name", "nt_rpm")]
colnames(Polyomavirus_s)[2]<-"Polyomavirus"
res<-merge(res,Polyomavirus_s, by.x="Barcode", by.y="barcode_name", all.x=TRUE)

######Combine bacteria and viruses
alltogether<-merge(bacteria, res, by.x = "Barcode", by.y="Barcode", all.x=TRUE)

#Remove NAs 
alltogether[is.na(alltogether)]<-0
#Remove duplicate rows
alltogether<- alltogether[!duplicated(alltogether$Barcode), ]

#Transform any positive number into a 1, meaning bug was detected
allbinary<-alltogether
for (i in (1:length(allbinary$Bacteriodes))) {if (allbinary$Bacteriodes[i] >0){allbinary$Bacteriodes[i] = 1}}
for (i in (1:length(allbinary$Citrobacter))) {if (allbinary$Citrobacter[i] >0){allbinary$Citrobacter[i] = 1}}
for (i in (1:length(allbinary$Escherichia))) {if (allbinary$Escherichia[i] >0){allbinary$Escherichia[i] = 1}}
for (i in (1:length(allbinary$Enterobacter))) {if (allbinary$Enterobacter[i] >0){allbinary$Enterobacter[i] = 1}}
for (i in (1:length(allbinary$Klebsiella))) {if (allbinary$Klebsiella[i] >0){allbinary$Klebsiella[i] = 1}}
for (i in (1:length(allbinary$Haemophilus))) {if (allbinary$Haemophilus[i] >0){allbinary$Haemophilus[i] = 1}}
for (i in (1:length(allbinary$Lactobacillus))) {if (allbinary$Lactobacillus[i] >0){allbinary$Lactobacillus[i] = 1}}
for (i in (1:length(allbinary$Morganella))) {if (allbinary$Morganella[i] >0){allbinary$Morganella[i] = 1}}
for (i in (1:length(allbinary$Prevotella))) {if (allbinary$Prevotella[i] >0){allbinary$Prevotella[i] = 1}}
for (i in (1:length(allbinary$Proteus))) {if (allbinary$Proteus[i] >0){allbinary$Proteus[i] = 1}}
for (i in (1:length(allbinary$Pseudomonas))) {if (allbinary$Pseudomonas[i] >0){allbinary$Pseudomonas[i] = 1}}
for (i in (1:length(allbinary$Serratia))) {if (allbinary$Serratia[i] >0){allbinary$Serratia[i] = 1}}
for (i in (1:length(allbinary$Staphylococcus))) {if (allbinary$Staphylococcus[i] >0){allbinary$Staphylococcus[i] = 1}}
for (i in (1:length(allbinary$Streptococcus))) {if (allbinary$Streptococcus[i] >0){allbinary$Streptococcus[i] = 1}}
for (i in (1:length(allbinary$Anelloviridae))) {if (allbinary$Anelloviridae[i] >0){allbinary$Anelloviridae[i] = 1}}
for (i in (1:length(allbinary$HepatitisB))) {if (allbinary$HepatitisB[i] >0){allbinary$HepatitisB[i] = 1}}
for (i in (1:length(allbinary$HSV1))) {if (allbinary$HSV1[i] >0){allbinary$HSV1[i] = 1}}
for (i in (1:length(allbinary$HSV2))) {if (allbinary$HSV2[i] >0){allbinary$HSV2[i] = 1}}
for (i in (1:length(allbinary$VZV))) {if (allbinary$VZV[i] >0){allbinary$VZV[i] = 1}}
for (i in (1:length(allbinary$HHV4))) {if (allbinary$HHV4[i] >0){allbinary$HHV4[i] = 1}}
for (i in (1:length(allbinary$HHV5))) {if (allbinary$HHV5[i] >0){allbinary$HHV5[i] = 1}}
for (i in (1:length(allbinary$HHV6))) {if (allbinary$HHV6[i] >0){allbinary$HHV6[i] = 1}}
for (i in (1:length(allbinary$HHV8))) {if (allbinary$HHV8[i] >0){allbinary$HHV8[i] = 1}}
for (i in (1:length(allbinary$Bocavirus))) {if (allbinary$Bocavirus[i] >0){allbinary$Bocavirus[i] = 1}}
for (i in (1:length(allbinary$Polyomavirus))) {if (allbinary$Polyomavirus[i] >0){allbinary$Polyomavirus[i] = 1}}

mortality<-metadata[,c(1,3)]
allbinary<-merge(allbinary, mortality, all.x=TRUE, all.y=FALSE)
write.csv(allbinary, "Outputs/Source_Data4_mortalitypresenceabsenceofpathogens.csv")

#Analyze each pathogen by: total numbers, % of surviving patients, pval of fisher test - output tabulated in Supplementary Table ****
Staphylococcus <- allbinary %>% dplyr::filter(Staphylococcus ==1)
table(Staphylococcus$Hospital_Death)
1 - mean(as.numeric(Staphylococcus$Hospital_Death))
mean(as.numeric(Staphylococcus$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Staphylococcus))

#p=0.0061 

HHV8 <- allbinary %>% dplyr::filter(HHV8 ==1)
table(HHV8$Hospital_Death)
1 - mean(as.numeric(HHV8$Hospital_Death))
mean(as.numeric(HHV8$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$HHV8))

Bacteriodes <- allbinary %>% dplyr::filter(Bacteriodes ==1)
table(Bacteriodes$Hospital_Death)
1 - mean(as.numeric(Bacteriodes$Hospital_Death))
mean(as.numeric(Bacteriodes$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Bacteriodes))

HHV5 <- allbinary %>% dplyr::filter(HHV5 ==1)
table(HHV5$Hospital_Death)
1 - mean(as.numeric(HHV5$Hospital_Death))
mean(as.numeric(HHV5$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$HHV5))
###p = 0,0096

HHV6 <- allbinary %>% dplyr::filter(HHV6 ==1)
table(HHV6$Hospital_Death)
1 - mean(as.numeric(HHV6$Hospital_Death))
mean(as.numeric(HHV6$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$HHV6))

Escherichia <- allbinary %>% dplyr::filter(Escherichia ==1)
table(Escherichia$Hospital_Death)
1 - mean(as.numeric(Escherichia$Hospital_Death))
mean(as.numeric(Escherichia$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Escherichia))

Klebsiella <- allbinary %>% dplyr::filter(Klebsiella ==1)
table(Klebsiella$Hospital_Death)
1 - mean(as.numeric(Klebsiella$Hospital_Death))
mean(as.numeric(Klebsiella$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Klebsiella))

Streptococcus <- allbinary %>% dplyr::filter(Streptococcus ==1)
table(Streptococcus$Hospital_Death)
1 - mean(as.numeric(Streptococcus$Hospital_Death))
mean(as.numeric(Streptococcus$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Streptococcus))

Anello <- allbinary %>% dplyr::filter(Anelloviridae ==1)
table(Anello$Hospital_Death)
1 - mean(as.numeric(Anello$Hospital_Death))
mean(as.numeric(Anello$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Anelloviridae))

Pseudo <- allbinary %>% dplyr::filter(Pseudomonas ==1)
table(Pseudo$Hospital_Death)
1 - mean(as.numeric(Pseudo$Hospital_Death))
mean(as.numeric(Pseudo$Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$Pseudomonas))

HHV4 <- allbinary %>% dplyr::filter(HHV4 ==1)
table(HHV4 $Hospital_Death)
1 - mean(as.numeric(HHV4 $Hospital_Death))
mean(as.numeric(HHV4 $Hospital_Death))
fisher.test(table(allbinary$Hospital_Death,allbinary$HHV4))

####the following pathogens are found in less than 3 patients and therefore fishers' tests cannot be done
Polyomavirus <- allbinary %>% dplyr::filter(Polyomavirus ==1)
table(Polyomavirus$Hospital_Death)
HepatitisB <- allbinary %>% dplyr::filter(HepatitisB ==1)
table(HepatitisB$Hospital_Death)
HSV <- allbinary %>% dplyr::filter(HSV1 ==1)
table(HSV$Hospital_Death)
Citrobacter <- allbinary %>% dplyr::filter(Citrobacter ==1)
table(Citrobacter$Hospital_Death)
Haemophilus <- allbinary %>% dplyr::filter(Haemophilus ==1)
table(Haemophilus$Hospital_Death)
Lactobacillus <- allbinary %>% dplyr::filter(Lactobacillus ==1)
table(Lactobacillus$Hospital_Death)
Morganella <- allbinary %>% dplyr::filter(Morganella ==1)
table(Morganella$Hospital_Death)
Prevotella <- allbinary %>% dplyr::filter(Prevotella ==1)
table(Prevotella$Hospital_Death)
Proteus <- allbinary %>% dplyr::filter(Proteus ==1)
table(Proteus$Hospital_Death)
Serratia <- allbinary %>% dplyr::filter(Serratia ==1)
table(Serratia$Hospital_Death)
HSV2 <- allbinary %>% dplyr::filter(HSV2 ==1)
table(HSV2$Hospital_Death)
VZV <- allbinary %>% dplyr::filter(VZV ==1)
table(VZV$Hospital_Death)



#########SECTION: Analyzing De genes by high and low pathogen dominance for Sepsis- becomes Figure 3D
# Read data
genecountspc<-read.csv("Inputs/earli_counts_kallisto_mortality.csv", header=TRUE)
genecountspc<-as.matrix(genecountspc)
rownames(genecountspc)<-genecountspc[,1]
genecountspc<-genecountspc[,-1]
genecountspc <- apply(genecountspc, c(1,2), function(x) { (as.integer(x))})
head(genecountspc)

####remove the 'earli' from the colnames 
colnames(genecountspc)<-substr(colnames(genecountspc), 7,11)

#unique(metadata$Group)
metadata<-read.csv("Inputs/Cleaned_Metadata_070324.csv")
output<-read.csv("Outputs/Supplementary_Data4_sepsis_microbialmass_and_pathdominance.csv")
output$Barcode<-substr(output$X, 7,11)
metadata_analysis<-metadata %>% filter(Group == "1_Sepsis+BldCx+"|Group == "2_Sepsis+OtherCx+")
merge<-merge(metadata_analysis, output, by.x = "Barcode", by.y="Barcode", all.x=TRUE, all.y=FALSE)
merge$rel_dom <- as.numeric(merge$rel_dom)
metadata_analysis <- merge 
metadata_analysis <- metadata_analysis[!is.na(metadata_analysis$rel_dom), ]

###scale losdeath Age
metadata_analysis$scaleAge <- scale(metadata_analysis$Age)

#Ensuring gender is a factor
metadata_analysis$Gender<-as.factor(metadata_analysis$Gender) 
rownames(metadata_analysis)<-metadata_analysis$Barcode

metadata<-metadata_analysis
counts<-genecountspc
x<- c(rownames(metadata))
y<-c(colnames(counts))
x
y
counts<-counts[,colnames(counts) %in% x]
metadata<-metadata[rownames(metadata) %in% y,]
counts<-counts[,order(noquote(colnames(counts)))]
metadata<-metadata[order(noquote(rownames(metadata))),]
identical(rownames(metadata), colnames(counts))

quantile(metadata$rel_dom, na.rm=T)

##top 1/2 is 0.6238125
#create a variable for top half
metadata$domhighlow<-NA

##assign pts with top half of microbial mass 
for (i in (1:nrow(metadata))){
  if(metadata$rel_dom[i] > 0.6238125){metadata$domhighlow[i]=1}
  else
  {metadata$domhighlow[i]=0}
}

metadata$domhighlow<-as.factor(metadata$domhighlow)

##Set up target category
target_cat <- "domhighlow" 

##tilda is formula ##plus is not interacting 
design <- model.matrix(as.formula(paste(paste("~", paste0(" 0 + ", target_cat)), paste0("+ scaleAge + Gender "))), data = metadata)
dds_paxgene_tmp <- DGEList(counts = counts)

#normalization factor, including library death and number of genes expressed. Linear model
dds_paxgene_tmp <- calcNormFactors(dds_paxgene_tmp)
dds_paxgene_tmp <- voom(dds_paxgene_tmp,design, plot=T)
dds_paxgene_tmp <- lmFit(dds_paxgene_tmp, design)

cont_mat <- makeContrasts(paste(paste0(target_cat, 1), paste0(target_cat, 0), sep=" - "), 
                          levels=design)

dds_paxgene_tmp <- contrasts.fit(dds_paxgene_tmp, cont_mat)
dds_paxgene_tmp <- eBayes(dds_paxgene_tmp, robust=TRUE)
res_paxgene <- topTable(dds_paxgene_tmp,
                        coef=paste(paste0(target_cat, 1), paste0(target_cat, 0), sep=" - "),
                        sort.by = "P", n = Inf)

# sort the genes from lowest to highest given adjusted p-values
res_paxgene <- res_paxgene[order(res_paxgene$adj.P.Val, decreasing = F), ]

# replace NA values with 1s and keep only significant genes
res_paxgene$adj.P.Val[is.na(res_paxgene$adj.P.Val)] <- 1

# filter
sig_res_paxgene <- data.frame(res_paxgene)
sig<-sig_res_paxgene %>% filter (adj.P.Val <0.1)
sigup<-sig %>% filter (logFC>0)
sigdown<-sig %>% filter (logFC<0)
dim(sig)

#Putting in the gene names to ensemble IDS
gene_attr<-read.csv("Inputs/gene_attr.csv", header=T)
head(gene_attr)
res_paxgene$Ensembl_ID<-rownames(res_paxgene)
res_table<- merge(gene_attr, res_paxgene, by.x="gene_id", by.y="Ensembl_ID", all.x=FALSE, all.y=TRUE)
head(res_table)

write.csv(res_table, "Outputs/Supplementary_Data3_DEGenes_pathdominance_highlow.csv")


######Making plots of micro abundance - Supplementary Figure 1
################this section makes the plots#####################
###define groups as part of the nt_rpm
tograph<-alltogether
tograph$Hospital_Death <- sapply(tograph$Barcode, function(x) metadata_full[metadata_full$Barcode==x, "Hospital_Death"])
tograph$Group <- sapply(tograph$Barcode, function(x) metadata_full[metadata_full$Barcode==x, "Group"])
tograph<-tograph %>% dplyr::filter(Hospital_Death =="1"|Hospital_Death == "0" )

tograph<-tograph[,c("Barcode", "Hospital_Death", "Group", "Bacteriodes", "Citrobacter", "Escherichia", "Enterobacter",
                    "Klebsiella", "Lactobacillus", "Morganella", "Prevotella", "Proteus", "Pseudomonas", "Serratia", "Staphylococcus",
                    "Streptococcus", "Anelloviridae", "HepatitisB", "HSV1", "HSV2", "VZV", "HHV4", "HHV5", "HHV6", 
                    "HHV8", "Bocavirus", "Polyomavirus")]

#Fix some colnames
colnames(tograph)[4]<-"Bacteriodes spp."
colnames(tograph)[5]<-"Citrobacter spp."
colnames(tograph)[6]<-"Escherichia coli"
colnames(tograph)[7]<-"Enterobacter spp."
colnames(tograph)[8]<-"Klebsiella spp."
colnames(tograph)[9]<-"Lactobacillus spp."
colnames(tograph)[10]<-"Morganella spp."
colnames(tograph)[11]<-"Prevotella spp."
colnames(tograph)[12]<-"Proteus spp."
colnames(tograph)[13]<-"Pseudomonas aeruginosa"
colnames(tograph)[14]<-"Serratia spp."
colnames(tograph)[15]<-"S. aureus complex"
colnames(tograph)[16]<-"Streptococcus spp."
colnames(tograph)[17]<-"Torque tenovirus"
colnames(tograph)[18]<-"Hepatitis B"
colnames(tograph)[19]<-"HSV-1"
colnames(tograph)[20]<-"HSV-2"
colnames(tograph)[21]<-"VZV"
colnames(tograph)[22]<-"HHV-4 (EBV)"
colnames(tograph)[23]<-"HHV-5 (CMV)"
colnames(tograph)[24]<-"HHV-6"
colnames(tograph)[25]<-"HHV-8"
colnames(tograph)[26]<-"Bocavirus"
colnames(tograph)[27]<-"Polyomavirus"

###choose patients with sepsis
group12<-tograph %>% dplyr::filter(Group =="1_Sepsis+BldCx+" |Group =="2_Sepsis+OtherCx+" )
group12_bin <- group12[, !(colnames(group12) %in% c("Barcode", "X", "Hospital_Death", "Group"))]
group12_bin <- group12_bin[ , hclust(dist(t(group12_bin), "euclidean"), "complete")$order]
group12_bin<-group12

rownames(group12_bin) <- group12_bin$Barcode

patho_levels <- colnames(group12_bin)[!(colnames(group12) %in% c("Barcode", "X", "Hospital_Death", "Group"))]
patho_levels <- rev(names(sort(colSums(group12_bin[, patho_levels]>0))))

group12_bin <- group12_bin[do.call(order, lapply(patho_levels, function(col) group12_bin[[col]] > 0)), ]

group12pivot <- group12_bin %>%
  tidyr::pivot_longer(cols=c("Bacteriodes spp.":"Polyomavirus"), names_to="pathogen", values_to="yes_no") %>%
  mutate(Barcode=as.character(Barcode))
group12pivot1 <- group12pivot %>%
  mutate(
    Outcome = dplyr::case_when(
      Hospital_Death == 0 ~ "Survived",
      Hospital_Death == 1 ~ "Deceased",
      TRUE ~ NA_character_))
rownames(group12) <- group12$Barcode
group12pivot2 <- group12pivot1

group12pivot2$pathogen <- factor(group12pivot2$pathogen, levels = names(sort(table(group12pivot2[group12pivot2$yes_no>0, ]$pathogen))))
colnames(group12_bin)
group12pivot2$Outcome <- factor(group12pivot2$Outcome, levels = c("Survived", "Deceased"))
survived_df <- group12pivot2[group12pivot2$Outcome=="Survived", ]
survived_hc <- hclust(dist(group12_bin[rownames(group12_bin) %in% unique(survived_df$Barcode), ], "euclidean"), "complete")$order
deceased_df <- group12pivot[group12pivot1$Outcome=="Deceased", ]
deceased_hc <- hclust(dist(group12_bin[rownames(group12_bin) %in% unique(deceased_df$Barcode), ], "euclidean"), "complete")$order

#group12pivot2$Barcode <- factor(group12pivot2$Barcode, 
#                                levels=c(rownames(group12_bin[rownames(group12_bin) %in% unique(survived_df$Barcode), ])[survived_hc], 
#                                         rownames(group12_bin[rownames(group12_bin) %in% unique(deceased_df$Barcode), ])[deceased_hc]))

group12pivot2$Barcode <- factor(group12pivot2$Barcode, levels=rev(unique(group12pivot2$Barcode)))
group12pivot2_1 <-group12pivot2 %>% filter(is.na(pathogen) == FALSE)

hm <- ggplot(group12pivot2_1,
             aes(x = Barcode,
                 y = pathogen,
                 fill = log(yes_no)))+ geom_tile() +
  scale_fill_gradient(low = "white", high ="purple", na.value="white")  + 
  facet_grid(~Outcome, scales = "free", space="free") +  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Pathogens") +xlab("Patients with sepsis")

gplot <- ggplotGrob(hm)
gplot$grobs[[10]]$grobs[[1]]$children[[1]]$gp$fill <- "steelblue"
gplot$grobs[[11]]$grobs[[1]]$children[[1]]$gp$fill<- "goldenrod"

require(gridExtra)
grid.arrange(gplot)
ggsave("Outputs/SuppFigure1.pdf")

#########################
# GENERATE ntRPM MATRIX
#########################
# G1+G2
# we only want samples with RNA-seq
# G1+G2
data_path <- "Inputs/"
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

metadata_bgc <- metadata_full[(metadata_full$Group %in% c("1_Sepsis+BldCx+", "2_Sepsis+OtherCx+")) & metadata_full$Barcode!="11942", ]

# kallisto filter
metadata_bgc <- metadata_bgc[rownames(metadata_bgc) %in% colnames(paxgene_cnt_data), ]

bgc_data <- bgc_data_full[bgc_data_full$EARLI_Barcode %in% metadata_bgc$EARLI_Barcode, ]

bgc_data_cast <- reshape2::dcast(bgc_data, tax_id~EARLI_Barcode, value.var="nt_rpm")
bgc_data_cast[is.na(bgc_data_cast)] <- 0
rownames(bgc_data_cast) <- bgc_data_cast$tax_id
bgc_data_cast <- bgc_data_cast[, colnames(bgc_data_cast)!="tax_id"]

# taxonomic data
tax_data <- tax_data_full[rownames(bgc_data_cast), ]

# with NHSN list
nhsn_list <-  read.csv("Inputs/nhsn_list.csv", header = T)

to_incl <- tax_data[(tax_data$species %in% nhsn_list[, 1]) | (tax_data$genus %in% nhsn_list[, 1]) | (tax_data$family %in% nhsn_list[, 1]) | (tax_data$kingdom %in% nhsn_list[, 1]), "taxon"]

# minus the new list
exclude_updated_list <-  read.csv("Outputs/Supplementary_Data18_EARLI_environmental_mNGS_exclude.csv", header = F)

to_incl_filt <- to_incl[!(to_incl %in% tax_data[tax_data$genus %in% exclude_updated_list$V1, "taxon"])]

ntrpm_bgc_data_cast_include_exclude_list <- bgc_data_cast[rownames(bgc_data_cast) %in% to_incl_filt, ]

write.csv(ntrpm_bgc_data_cast_include_exclude_list, "Outputs/ntrpm_bgc_data_cast_include_exclude_list.csv")

# repeat for G4
metadata_bgc <- metadata_full[(metadata_full$Group %in% c("4_NO_Sepsis")), ]

# kallisto filter
# we only want samples with RNA-seq
metadata_bgc <- metadata_bgc[rownames(metadata_bgc) %in% colnames(paxgene_cnt_data), ]

bgc_data <- bgc_data_full[bgc_data_full$EARLI_Barcode %in% metadata_bgc$EARLI_Barcode, ]

# duplicate
#bgc_data<- bgc_data[bgc_data$sample_name!="EARLI_11255_PLaD1_L6_OPS046_S92", ]
bgc_data<- bgc_data[bgc_data$sample_name!="EARLI_11255_PLaD1_A10_OPS046_S145", ]

bgc_data_cast <- dcast(bgc_data, tax_id~EARLI_Barcode, value.var="nt_rpm")
bgc_data_cast[is.na(bgc_data_cast)] <- 0
rownames(bgc_data_cast) <- bgc_data_cast$tax_id
bgc_data_cast <- bgc_data_cast[, colnames(bgc_data_cast)!="tax_id"]

# taxonomic data
tax_data <- tax_data_full[rownames(bgc_data_cast), ]

# with NHSN list
to_incl <- tax_data[(tax_data$species %in% nhsn_list[, 1]) | (tax_data$genus %in% nhsn_list[, 1]) | (tax_data$family %in% nhsn_list[, 1]) | (tax_data$kingdom %in% nhsn_list[, 1]), "taxon"]

# minus the new list
to_incl_filt <- to_incl[!(to_incl %in% tax_data[tax_data$genus %in% exclude_updated_list$V1, "taxon"])]

ntrpm_bgc_data_cast_include_exclude_list_G4 <- bgc_data_cast[rownames(bgc_data_cast) %in% to_incl_filt, ]

to_fill <- rownames(ntrpm_bgc_data_cast_include_exclude_list)[!(rownames(ntrpm_bgc_data_cast_include_exclude_list) %in% rownames(ntrpm_bgc_data_cast_include_exclude_list_G4))]

for (row_ in to_fill){
  ntrpm_bgc_data_cast_include_exclude_list_G4[row_, ] <- rep(0, ncol(ntrpm_bgc_data_cast_include_exclude_list_G4))
}

write.csv(ntrpm_bgc_data_cast_include_exclude_list_G4, "Outputs/ntrpm_bgc_data_cast_include_exclude_list_G4.csv")

