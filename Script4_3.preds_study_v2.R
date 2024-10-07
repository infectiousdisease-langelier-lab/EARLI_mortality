
library(rpart)
library(pROC)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(data.table)

data_path <- paste0(getwd(), "/Inputs/")
results_path <- paste0(getwd(), "/Outputs/")

group_folder <- "G1G2_Hospital_Death"

pred_table <- read.csv(paste0(results_path, paste0("pred_table_bsvm_dgea_0.1_100_0.1_1.0_50_", paste0(group_folder, ".csv"))), row.names = 1)

full_meta_data <- read.csv(paste0(data_path, "Cleaned_Metadata_070324.csv"),
                           row.names=1)

# SRSq values
srsq_data <-  read.csv(paste0(paste0(results_path, "Sepsis"), "_SRSq.csv"), row.names = 1)

# death scores
death_scores_data <-  read.csv(paste0(paste0(results_path, "Sepsis"), "_death_scores.csv"), row.names = 1)

auc_list <- list()

###############################
# functions
###############################
build_survival_curve <- function(meta_data, y_pred, death_day, file_name){

  meta_data$death_day <- meta_data[, death_day]
  meta_data$censored[is.na(meta_data$death_day)] <- 1
  meta_data$censored[is.na(meta_data$censored)] <- 2
  #meta_data$death_day[is.na(meta_data$death_day)] <- max(meta_data$death_day, na.rm = T)+1
  meta_data$death_day[is.na(meta_data$death_day)] <- 91
  
  meta_data[, "y_pred"] <- factor(meta_data[, y_pred], levels=c("1", "0"))
  
  meta_data <- meta_data[!is.na(meta_data$y_pred), ]
  
  p <- surv_fit(Surv(death_day, censored) ~ y_pred, data = meta_data)
  
  ggsurv <- ggsurvplot(p, xlim=c(0, 89), risk.table = T, pval = T, break.x.by=10,
                   palette = c("#c9993c", "#446ea0"))
  
  pdf(str_replace(paste0(results_path, file_name), ".csv", ".pdf"))
  print(ggsurv)
  dev.off()
  
  write.csv(ggsurv$data.survplot, paste0(results_path, file_name))
  
}

# add additional variables
meta_data_paxgene <- full_meta_data[rownames(full_meta_data) != "11942", ]
rownames(meta_data_paxgene) <- paste0("EARLI_", rownames(meta_data_paxgene))

###############################
# bsvm paxgene
###############################
# mean prob per sample
agg_table <- aggregate(y_pred_prob~samples, pred_table, mean)

rownames(agg_table) <- agg_table$samples
agg_table$y_pred[agg_table$y_pred_prob<0.5] <- 0
agg_table$y_pred[agg_table$y_pred_prob>=0.5] <- 1

meta_data_paxgene_all <- meta_data_paxgene
meta_data_paxgene <- meta_data_paxgene[rownames(agg_table), ]

meta_data_paxgene$y_pred <- agg_table$y_pred
meta_data_paxgene$y_pred_prob <- agg_table$y_pred_prob

# censor at 90?
meta_data_paxgene$deathday <- meta_data_paxgene$daysofsurvival

# allocation prob
agg_table_tmp <- agg_table
agg_table_tmp$y_pred_prob_max <- agg_table_tmp$y_pred_prob
agg_table_tmp$y_pred_prob_max[agg_table_tmp$y_pred_prob_max<0.5] <- 1-agg_table_tmp$y_pred_prob_max[agg_table_tmp$y_pred_prob_max<0.5]

meta_data_paxgene$y_pred_prob_max <- agg_table_tmp[rownames(meta_data_paxgene), "y_pred_prob_max"]

# what is unique about misclassified samples?
meta_data_paxgene$correctly_pred <- (meta_data_paxgene$y_pred == meta_data_paxgene$Hospital_Death)

# logistic regression
meta_data_paxgene$y_pred_prob <- agg_table_tmp[rownames(meta_data_paxgene), "y_pred_prob"]

# cart model
meta_data_paxgene$Hospital_Death <- factor(meta_data_paxgene$Hospital_Death)
meta_data_paxgene$Group <- factor(meta_data_paxgene$Group)

###############################
# LR clinical model - apacheiii
###############################
# APACHE III only
meta_data_clin_apacheiii_lr <- meta_data_paxgene
meta_data_clin_apacheiii_lr_G4 <- meta_data_paxgene_all
meta_data_clin_apacheiii_lr_G4 <- meta_data_clin_apacheiii_lr_G4[meta_data_clin_apacheiii_lr_G4$Group=="4_NO_Sepsis", ]

# 34 variables
cont_cols_toinclude <- c("APACHEIII")

extra_cols_toinclude <- c("Hospital_Death")
meta_data_clin_apacheiii_lr <- meta_data_clin_apacheiii_lr[, c(cont_cols_toinclude, extra_cols_toinclude)]
meta_data_clin_apacheiii_lr_G4 <- meta_data_clin_apacheiii_lr_G4[, c(cont_cols_toinclude, extra_cols_toinclude)]

meta_data_clin_apacheiii_lr$Hospital_Death <- as.factor(meta_data_clin_apacheiii_lr$Hospital_Death)
meta_data_clin_apacheiii_lr_G4$Hospital_Death <- as.factor(meta_data_clin_apacheiii_lr_G4$Hospital_Death)

# do it per fold
auc_vals <- c()
auc_vals_G4 <- c()
preds_clin_apacheiii_lr <- c()
preds_clin_apacheiii_lr_G4 <- c()
for (i in c(1:10)){
  # read test labels
  test_probs_file <- read.csv(paste0(results_path, paste0("test_labels_", paste0(paste(i, group_folder, sep="_"), ".csv"))), 
                              row.names = 1)
  test_samples <- test_probs_file$X0
  
  meta_data_train <- meta_data_clin_apacheiii_lr[!(rownames(meta_data_clin_apacheiii_lr) %in% test_samples), ]
  meta_data_train <- meta_data_train[, !(colnames(meta_data_train) %in% colnames(meta_data_train)[apply(meta_data_train, 2, function(x) (length(unique(x)))) < 2])]
  
  meta_data_test <- meta_data_clin_apacheiii_lr[rownames(meta_data_clin_apacheiii_lr) %in% test_samples, ]
  
  meta_data_train$APACHEIII <- scale(meta_data_train$APACHEIII)
  meta_data_test$APACHEIII <- scale(meta_data_test$APACHEIII, 
                              attr(meta_data_train$APACHEIII, "scaled:center"), 
                              attr(meta_data_train$APACHEIII, "scaled:scale"))
  
  meta_data_G4 <- meta_data_clin_apacheiii_lr_G4
  meta_data_G4$APACHEIII <- scale(meta_data_G4$APACHEIII, 
                                          attr(meta_data_train$APACHEIII, "scaled:center"), 
                                          attr(meta_data_train$APACHEIII, "scaled:scale"))
  
  # with pred rather than pred prob
  if (i==1){
    print(c("Hospital_Death", "APACHEIII"))
  }
  model <- glm(Hospital_Death ~ APACHEIII, data=meta_data_train, family="binomial")
  print(((-1*coef(model)[1]/coef(model)[2]) * attr(meta_data_train$APACHEIII, "scaled:scale")) + attr(meta_data_train$APACHEIII, "scaled:center"))
  
  # create roc curve
  preds_clin_apacheiii_lr <- c(preds_clin_apacheiii_lr, stats::predict(model, meta_data_test, type="response"))
  preds_clin_apacheiii_lr_G4 <- c(preds_clin_apacheiii_lr_G4, stats::predict(model, meta_data_G4, type="response"))
  
  roc_object <- roc(meta_data_test$Hospital_Death, stats::predict(model, meta_data_test, type="response"))
  roc_object_G4 <- roc(meta_data_G4$Hospital_Death, stats::predict(model, meta_data_G4, type="response"))
  
  # calculate area under curve
  auc_vals <- c(auc_vals, auc(roc_object))
  auc_vals_G4 <- c(auc_vals_G4, auc(roc_object_G4))
  
}

print("apacheiii")
print(mean(auc_vals))
print(sd(auc_vals))

print("apacheiii G4")
print(mean(auc_vals_G4))
print(sd(auc_vals_G4))

auc_list[['apacheiii']] <- auc_vals

###############################
# LR clinical model - SRSq
###############################
meta_data_paxgene$SRSq <- srsq_data[rownames(meta_data_paxgene), ]

# APACHE III only
meta_data_clin_apacheiii_lr <- meta_data_paxgene

# 34 variables
cont_cols_toinclude <- c("SRSq")

extra_cols_toinclude <- c("Hospital_Death")
meta_data_clin_apacheiii_lr <- meta_data_clin_apacheiii_lr[, c(cont_cols_toinclude, extra_cols_toinclude)]

meta_data_clin_apacheiii_lr$Hospital_Death <- as.factor(meta_data_clin_apacheiii_lr$Hospital_Death)

# do it per fold
auc_vals <- c()
preds_clin_apacheiii_lr <- c()
for (i in c(1:10)){
  # read test labels
  test_probs_file <- read.csv(paste0(results_path, paste0("/test_labels_", paste0(paste(i, group_folder, sep="_"), ".csv"))), 
                              row.names = 1)
  test_samples <- test_probs_file$X0
  
  meta_data_train <- meta_data_clin_apacheiii_lr[!(rownames(meta_data_clin_apacheiii_lr) %in% test_samples), ]
  meta_data_train <- meta_data_train[, !(colnames(meta_data_train) %in% colnames(meta_data_train)[apply(meta_data_train, 2, function(x) (length(unique(x)))) < 2])]
  
  meta_data_test <- meta_data_clin_apacheiii_lr[rownames(meta_data_clin_apacheiii_lr) %in% test_samples, ]
  
  # with pred rather than pred prob
  if (i==1){
    print(c("Hospital_Death", "SRSq"))
  }
  model <- glm(Hospital_Death ~ SRSq, data=meta_data_train, family="binomial")
  
  # create roc curve
  preds_clin_apacheiii_lr <- c(preds_clin_apacheiii_lr, stats::predict(model, meta_data_test, type="response"))
  
  roc_object <- roc(meta_data_test$Hospital_Death, stats::predict(model, meta_data_test, type="response"))
  
  # calculate area under curve
  auc_vals <- c(auc_vals, auc(roc_object))
  
}

print("SRSq")
print(mean(auc_vals))
print(sd(auc_vals))

auc_list[['SRSq']] <- auc_vals

###############################
# LR clinical model - death scores
###############################
meta_data_paxgene[, colnames(death_scores_data)] <- death_scores_data[rownames(meta_data_paxgene), ]

# 34 variables
for (death_score in c("death1_scores", "death2_scores","death3_scores","death4_scores")){
  print(death_score)
  
  # APACHE III only
  meta_data_clin_apacheiii_lr <- meta_data_paxgene
  
  cont_cols_toinclude <- death_score
  meta_data_clin_apacheiii_lr$death_scores <- meta_data_clin_apacheiii_lr[, cont_cols_toinclude]
  
  extra_cols_toinclude <- c("Hospital_Death")
  meta_data_clin_apacheiii_lr <- meta_data_clin_apacheiii_lr[, c("death_scores", extra_cols_toinclude)]
  
  meta_data_clin_apacheiii_lr$Hospital_Death <- as.factor(meta_data_clin_apacheiii_lr$Hospital_Death)
  
  # do it per fold
  auc_vals <- c()
  preds_clin_apacheiii_lr <- c()
  for (i in c(1:10)){
    # read test labels
    test_probs_file <- read.csv(paste0(results_path, paste0("/test_labels_", paste0(paste(i, group_folder, sep="_"), ".csv"))), 
                                row.names = 1)
    test_samples <- test_probs_file$X0
    
    meta_data_train <- meta_data_clin_apacheiii_lr[!(rownames(meta_data_clin_apacheiii_lr) %in% test_samples), ]
    meta_data_train <- meta_data_train[, !(colnames(meta_data_train) %in% colnames(meta_data_train)[apply(meta_data_train, 2, function(x) (length(unique(x)))) < 2])]
    
    meta_data_test <- meta_data_clin_apacheiii_lr[rownames(meta_data_clin_apacheiii_lr) %in% test_samples, ]
    
    # with pred rather than pred prob
    if (i==1){
      print(c("Hospital_Death", cont_cols_toinclude))
    }
    model <- glm(Hospital_Death ~ death_scores, data=meta_data_train, family="binomial")
    
    # create roc curve
    preds_clin_apacheiii_lr <- c(preds_clin_apacheiii_lr, stats::predict(model, meta_data_test, type="response"))
    
    roc_object <- roc(meta_data_test$Hospital_Death, stats::predict(model, meta_data_test, type="response"))
    
    # calculate area under curve
    auc_vals <- c(auc_vals, auc(roc_object))
    
  }
  
  print(mean(auc_vals))
  print(sd(auc_vals))
  
  auc_list[[death_score]] <- auc_vals
}

###############################
# BGC MODEL TOP FEATURES
###############################
bgc_pred_table <- read.csv("Outputs/bgc_pred_table_bsvm_dgea_1_100_0.9_1.0_50_incl_excl_filt_ntrpm_original_set_G1G2_Hospital_Death.csv", row.names = 1)

bgc_pred_table_full <- bgc_pred_table

samples_in_common <- intersect(bgc_pred_table$samples, pred_table$samples)

rownames(bgc_pred_table) <- bgc_pred_table$samples
rownames(pred_table) <- pred_table$samples

bgc_pred_table <- bgc_pred_table[samples_in_common, ]
pred_table <- pred_table[samples_in_common, ]

# best variables
# taxonomic data from NCBI
tax_data <- readRDS("Inputs/tax_data.RDS")

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

tax_data_df <- data.frame(tax_data)
rownames(tax_data_df) <- tax_data_df$taxon

top_10_i <- c()
ranks_i <- list()
for (i in c(1:10)){
  best_vars_i <- read.csv(paste0("Outputs/top_vars_bsvm_dgea_1_100_0.9_1.0_50_incl_excl_filt_ntrpm_original_set_", paste0(paste(i, group_folder, sep="_"), ".csv")), 
                          row.names = 1)

  best_vars_i_mean <- best_vars_i %>% 
        dplyr::group_by(best_vars)  %>% 
        dplyr::summarise(mean(best_coefs))
  
  print(tax_data_df[as.character(best_vars_i_mean[order(abs(best_vars_i_mean$`mean(best_coefs)`), decreasing = T), "best_vars"]$best_vars), "species"][1:10])
  
  top_10_i <- c(top_10_i, tax_data_df[as.character(best_vars_i_mean[order(abs(best_vars_i_mean$`mean(best_coefs)`), decreasing = T), "best_vars"]$best_vars), "species"][1:10])
  
  best_vars_i_mean_df <- as.data.frame(best_vars_i_mean)
  best_vars_i_mean_df <- best_vars_i_mean_df[order(abs(best_vars_i_mean_df$`mean(best_coefs)`), decreasing = T), ]
  best_vars_i_mean_df$species <- tax_data_df[as.character(best_vars_i_mean_df$best_vars), "species"]
  
  best_vars_i_mean_df$rank_inv <- rev(c(0:(nrow(best_vars_i_mean_df)-1)))
  
  ranks_i[[i]] <- best_vars_i_mean_df
  
}

table(top_10_i)

all_microbes <- unique(unlist(lapply(ranks_i, function(x) x$species)))
ranksum_microbes <- c()
means_microbes <- c()
for (microbe_ in all_microbes){
  ranksum_microbes <- c(ranksum_microbes, sum(unlist(lapply(ranks_i, function(x) x[x$species==microbe_, "rank_inv"]))))
  means_microbes <- c(means_microbes, mean(unlist(lapply(ranks_i, function(x) x[x$species==microbe_, "mean(best_coefs)"]))))
}
names(ranksum_microbes) <- all_microbes
names(means_microbes) <- all_microbes

means_microbes <- means_microbes[order(ranksum_microbes, decreasing = T)]
ranksum_microbes <- ranksum_microbes[order(ranksum_microbes, decreasing = T)]

write.csv(data.frame(ranksum_microbes, means_microbes), "Outputs/Supplementary_Data5_microbial_input_transcriptomic_classifier.csv")

###############################
# combine probs from host and bgc svm classifiers
###############################
# consider G4
bgc_pred_table_g4 <- read.csv("Outputs/bgc_pred_table_G4_bsvm_dgea_1_100_0.9_1.0_50_incl_excl_filt_ntrpm_original_set_G1G2_Hospital_Death.csv", row.names = 1)

agg_table_g4 <- aggregate(y_pred_prob~samples, bgc_pred_table_g4, mean)

rownames(agg_table_g4) <- agg_table_g4$samples
agg_table_g4$y_pred[agg_table_g4$y_pred_prob<0.5] <- 0
agg_table_g4$y_pred[agg_table_g4$y_pred_prob>=0.5] <- 1

meta_data_paxgene_g4 <- meta_data_paxgene_all
meta_data_paxgene_g4 <- meta_data_paxgene_g4[rownames(agg_table_g4), ]

meta_data_paxgene_g4$bgc_pred <- agg_table_g4$y_pred
meta_data_paxgene_g4$bgc_pred_prob <- agg_table_g4$y_pred_prob

pred_table_g4 <- read.csv("Outputs/pred_table_g4_bsvm_dgea_0.1_100_0.1_1.0_50_G1G2_Hospital_Death.csv", row.names = 1)

agg_table_g4 <- aggregate(y_pred_prob~samples, pred_table_g4, mean)

rownames(agg_table_g4) <- agg_table_g4$samples
agg_table_g4$y_pred[agg_table_g4$y_pred_prob<0.5] <- 0
agg_table_g4$y_pred[agg_table_g4$y_pred_prob>=0.5] <- 1

meta_data_paxgene_g4 <- meta_data_paxgene_g4[rownames(meta_data_paxgene_g4) %in% rownames(agg_table_g4), ]

meta_data_paxgene_g4$y_pred <- agg_table_g4$y_pred
meta_data_paxgene_g4$y_pred_prob <- agg_table_g4$y_pred_prob

bgc_pred_table <- bgc_pred_table_full
rownames(bgc_pred_table) <- bgc_pred_table$samples

meta_data_paxgene$bgc_pred_prob <- bgc_pred_table[rownames(meta_data_paxgene), "y_pred_prob"]

pred_cols <- c("bgc_pred_prob", "y_pred_prob")

meta_data_comb <- meta_data_paxgene
meta_data_comb <- meta_data_comb[, c(pred_cols, "Hospital_Death")]

meta_data_comb <- meta_data_comb[rownames(meta_data_comb)!="EARLI_11942", ]

# do it per fold
auc_vals <- c()
preds_clin_shannon_top_patho_lr <- c()

auc_vals_g4 <- c()

auc_vals_host <- c()
auc_vals_microbe <- c()

for (i in c(1:10)){
  # read test labels
  test_probs_file <- read.csv(paste0(results_path, paste0("test_labels_"), paste0(paste(i, group_folder, sep="_"), ".csv")), 
                              row.names = 1)
  test_samples <- test_probs_file$X0
  
  meta_data_train <- meta_data_comb[!(rownames(meta_data_comb) %in% test_samples), ]
  meta_data_test <- meta_data_comb[rownames(meta_data_comb) %in% test_samples, ]

  # create roc curve
  formula_ <- as.formula(paste0("Hospital_Death ~ ", paste(pred_cols, collapse =" + ")))
  
  model <- glm(formula_, data=meta_data_train, family="binomial")
  
  # create roc curve
  preds_clin_shannon_top_patho_lr <- c(preds_clin_shannon_top_patho_lr, stats::predict(model, meta_data_test, type="response"))
  
  roc_object <- roc(meta_data_test$Hospital_Death, stats::predict(model, meta_data_test, type="response"))
  roc_object_g4 <- roc(meta_data_paxgene_g4$Hospital_Death, stats::predict(model, meta_data_paxgene_g4, type="response"))

  # calculate area under curve
  auc_vals <- c(auc_vals, auc(roc_object))
  auc_vals_g4 <- c(auc_vals_g4, auc(roc_object_g4))
  
  # for individual modalities
  auc_vals_host <- c(auc_vals_host, auc(roc(meta_data_test$Hospital_Death, meta_data_test$y_pred_prob)))
  auc_vals_microbe <- c(auc_vals_microbe, auc(roc(meta_data_test$Hospital_Death, meta_data_test$bgc_pred_prob)))
  
}

print("mean")
print(mean(auc_vals))
print(sd(auc_vals))

auc_list[['host_microbe']] <- auc_vals
auc_list[['host']] <- auc_vals_host
auc_list[['microbe']] <- auc_vals_microbe

print(mean(auc_vals_g4))
print(sd(auc_vals_g4))

meta_data_comb$preds_comb_lr <- sapply(preds_clin_shannon_top_patho_lr[rownames(meta_data_comb)],
                                                              function(x){
                                                                if (!is.na(x)){
                                                                  if (x>=0.5){
                                                                    "1"
                                                                  }else{
                                                                    "0"
                                                                  }
                                                                }else{
                                                                  NA
                                                                }
                                                              })
meta_data_comb$deathday <- meta_data_paxgene[rownames(meta_data_comb), ]$deathday

file_name <- "Source_Data6_KaplanMeier_hostpathogen.csv"
build_survival_curve(meta_data_comb,
                     "preds_comb_lr",
                     "deathday",
                     file_name)

###############################
# SUMMARY AUC FIGURE
################################
write.csv(data.frame(mean=unlist(lapply(auc_list, mean)), sd=unlist(lapply(auc_list, sd))), row.names = names(auc_list),
          paste0(results_path, "Source_Data5_AUCforDifferentClassifiers.csv"))


