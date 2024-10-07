    # EARLI Mortality - Host and Microbe Prognostic Biomarkers

This repository contains the code for the paper **Molecular profiling identifies sepsis prognostic biomarkers based on host and microbe**.

# Table of contents

1. [Code](#code)
2. [Code usage](#code-usage)
3. [Figures](#figures)
    1. [Main Figures](#main-figures)
    2. [Supplementary Figures](#supplementary-figures)
4. [Supplementary Data](#supplementary-tables)
5. [Data files](#data-files)
6. [Required hardware and software](#required-hardware-and-software)

## Code
[Script1_EARLISurvivalvsMortality_HostAnalyses.R](Script1_EARLISurvivalvsMortality_HostAnalyses.R) is used to perform host analyses at the gene expression and protein level between patients who survived vs. those who did not. 
* [earli_counts_kallisto_mortality.csv](Inputs/earli_counts_kallisto_mortality.csv): Gene expression counts for patients in this cohort 
* [Cleaned_Metadata_070324.csv](Inputs_Final/Cleaned_Metadata_070324.csv): Metadata and cytokine values for patients in this cohort 
* [gene_attr.csv](Inputs_Final/gene_attr.csv): list of gene names with ENSEMBL IDs  

[Script2_de_external_cohorts.R](Script2_de_external_cohorts.R) script is designed to perform differential expression analysis between deceased and surviving patients in the GAINS cohort.
* [geneLevel_data.csv](Inputs_Final/geneLevel_data.csv): Gene expression counts for patients in the GAINS cohort
* [class_labels.csv](Inputs_Final/class_labels.csv): Metadata for patients in the GAINS cohort

[Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R](Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R) script is designed to perform microbial analyses between deceased and surviving patients.
* [Cleaned_Metadata_070324.csv](Inputs_Final/Cleaned_Metadata_070324.csv): Metadata for patients in this cohort 
* [bgccounts_mortality.csv](Inputs_Final/bgccounts_mortality.csv): Microbial counts for patients in this cohort 
* [tax_data.RDS](Inputs_Final/tax_data.RDS): NCBI taxonomic information for detected microbes

[Script4_0.build_splits_grp_paxgene.py](Script4_0.build_splits_grp_paxgene.py) script building the folds for the classifiers
[Script4_1.de_script_grp_paxgene.R](Script4_1.de_script_grp_paxgene.R) script running differential gene expression analyses within each fold
[Script4_2.build_classifiers_bsvm_grp_paxgene.py](Script4_2.build_classifiers_bsvm_grp_paxgene.py) builds host classifiers using folds and lists of differentially expressed genes
[bgc_classifiers.py](bgc_classifiers.py) builds microbe classifiers using folds generated in Script4_0.build_splits_grp_paxgene.py and filtered microbial data generated from Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R
[Script4_3.preds_study_v2.R](Script4_3.preds_study_v2.R) script combines outputs from Script4_2.build_classifiers_bsvm_grp_paxgene.py and bgc_classifiers.py to build a combined host and microbe model
All Script4 scripts use:
* [earli_counts_kallisto_mortality.csv](Inputs/earli_counts_kallisto_mortality.csv): Gene expression counts for patients in this cohort 
* [Cleaned_Metadata_070324.csv](Inputs_Final/Cleaned_Metadata_070324.csv): Metadata and cytokine values for patients in this cohort 

[custom_classes.py](custom_classes.py) utility functions used in Script4_2.build_classifiers_bsvm_grp_paxgene.py and bgc_classifiers.py scripts

## Code usage
The full data sets are available in the folder [Inputs](Inputs). To reproduce the results in the paper, the uploaded codes could be run directly. The codes' output are located in the folder [Outputs](Outputs).

## Figures

### Main Figures

Figure 1. Figure 1 created in Powerpoint. Source data is available in Cleaned_Metadata_070324.csv

Figure 2: Created using script Script1_EARLISurvivalvsMortality_HostAnalyses.R. Figure 2A was created in R and then edited in Prism (font size and other minor changes for appearance). Data for Figure 2B (3 panels) was exported as CSV and graphs made in Prism. Gene expression data for Figure 2C is exported from Script1_EARLISurvivalvsMortality_HostAnalyses.R in as Supplementary_Data1_Sepsis_DEGenes.csv; and then ORA pathways identified in WebGestalt with output in Supplementary_Data2_Sepsis_ORAPathways.csv; graphs made in Prism. Raw data for Figure 2D is in Source_Data2_CytokineandOutcomeData.csv and input data is in Table 4, then graphed in Prism. 
* A: [Sepsis_volcano_plot.pdf](Outputs/Sepsis_volcano_plot.csv). 
* B: [Source_Data1_SpecificGenesSepsis.csv](Outputs/Source_Data1_SpecificGenesSepsis.csv)
* C: [Supplementary_Data2_Sepsis_ORAPathways.csv](Outputs/Supplementary_Data2_Sepsis_ORAPathways.csv)
* D: [Table4.docx](Outputs/Table4.docx)

Figure 3: Created using script Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R and graphed in Prism with the following exceptions. Figure 3C was created as .pdf by R and edited in Prism (font size and other minor changes for appearance). Gene expression data for Figure 3D is exported from Script3_EARLISurvivalvsMortality_MicrobialAnalyses.csv in as Supplementary_Data3_DEGenes_pathdominance_highlow.csv; and then ORA pathways identified in WebGestalt with output in Supplementary_Data4_ORAPathways_SepsisPatients_bypathogendominance_highlow.csv; graphs made in Prism. Figure 3E had data created in R and tabulated in excel. Final edits made in Illustrator. 
* A: [Source_Data3_sepsis_microbialmass_and_pathdominance](Outputs/Source_Data3_sepsis_microbialmass_and_pathdominance.csv). 
* B: [Source_Data3_sepsis_microbialmass_and_pathdominance](Outputs/Source_Data3_sepsis_microbialmass_and_pathdominance.csv). 
* C: [bacteria_Sepsis_genus_nmds_boxplot.pdf](Outputs/bacteria_Sepsis_genus_nmds_boxplot.pdf) 
* D: [Supplementary_Data4_ORAPathways_SepsisPatients_bypathogendominance_highlow.csv](Outputs/Supplementary_Data4_ORAPathways_SepsisPatients_bypathogendominance_highlow.csv)
* E: [Source_Data4_mortalitypresenceabsenceofpathogens.csv](Outputs/Source_Data4_mortalitypresenceabsenceofpathogens.csv)

Figure 4: Generated using script Script4_3.preds_study_v2.R. Figure 4A data exported as a .csv and graph made in Prism. Figure 4B exported as a R graph and edits made in Illustrator (including moving patients at risk on timeline to graph, cosmetic adjustments to font)
* A: [Source_Data5_AUCforDifferentClassifiers.csv](Outputs/Source_Data5_AUCforDifferentClassifiers.csv). 
* B: [Figure5B_KaplanMeier_hostpathogen.pdf](Outputs/Figure5B_KaplanMeier_hostpathogen.pdf). Raw numbers are available in [Source_Data6_KaplanMeier_hostpathogen.csv](Outputs/Source_Data6_KaplanMeier_hostpathogen.csv)

Figure 5: A was created using script Script1_EARLISurvivalvsMortality_HostAnalyses.R and  script Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R. Gene expression data for Figure 5A is exported from Script1_EARLISurvivalvsMortality_HostAnalyses.R in as Supplementary_Data8_ICUControls_DEGenes.csv; and then ORA pathways identified in WebGestalt with output in Supplementary_Data9_ICUControls_ORAPathways.csv; graphs made in Prism. Figure 5B data exported as Source_Data7_DEGene_Overlap.csv and graph made in Powerpoint. Raw data for Figure 5C is in Source_Data2_CytokineandOutcomeData.csv and input data is in Table 4, then graphed in Prism. Figures D-F are all created by script Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R and output is available in Source_Data8_ICUControls_microbialmass_and_pathdominance.csv 
* A: [Supplementary_Data9_ICUControls_ORAPathways.csv](Outputs/Supplementary_Data9_ICUControls_ORAPathways.csv)
* B: [Source_Data7_DEGene_Overlap.csv](Outputs/Source_Data7_DEGene_Overlap.csv). 
* C: [Table4.docx](Outputs/Table4.docx)
* D: [Source_Data8_ICUControls_microbialmass_and_pathdominance.csv](Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv)
* E: [Source_Data8_ICUControls_microbialmass_and_pathdominance.csv](Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv)
* F: [Source_Data8_ICUControls_microbialmass_and_pathdominance.csv](Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv) 


### Supplementary Figures
Supplementary Figure 1. Raw detection of microbes. Generated by code Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R and exported as SuppFigure1.pdf. Cosmetic edits (font changes) performed in Illustrator.


Supplementary Figure 2 (comparison to outside data sets) Generated using scripts Script2_de_external_cohorts.R for Supplementary Figure 2A and Script4_3.preds_study_v2.R for Supplementary 2B. Data exported as .csv files and graphs made in Prism. 
* A: [Supplementary_Data6_ORAPathways_GAINS.csv](Outputs/Supplementary_Data6_ORAPathways_GAINS.csv). 
* B: [Source_Data5_AUCforDifferentClassifiers.csv](Outputs/Source_Data5_AUCforDifferentClassifiers.csv)


Supplementary Figure 3. ICU control gene expression data with cardiac arrest correction. Generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R. Gene expression data for Supp Fig 3 is exported from Script1_EARLISurvivalvsMortality_HostAnalyses.csv in as Supplementary_Data10_ICUControls_DEGenes_CardiacArrestasCovariate.csv; and then ORA pathways identified in WebGestalt with output in Supplementary_Data11_ICUControls_ArrestAsCofactor_ORAPathways.csv; graphs made in Prism.
* [Supplementary_Data11_ICUControls_ArrestAsCofactor_ORAPathways.csv](Outputs_Final/Supplementary_Data11_ICUControls_ArrestAsCofactor_ORAPathways.csv) 

Supplementary Figure 4 (additional ICU control microbiology data) Generated using script Script1_EARLISurvivalvsMortality_HostAnalyses.R. Graphs generated in Prism.
* A: [Source_Data8_ICUControls_microbialmass_and_pathdominance.csv](Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv). 
* B: [Source_Data8_ICUControls_microbialmass_and_pathdominance.csv](Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv) 
* C: [Source_Data8_ICUControls_microbialmass_and_pathdominance.csv](Outputs/Source_Data8_ICUControls_microbialmass_and_pathdominance.csv) 


## Figures
### Main Tables
Table 1. Created in Word. Source data is available in Cleaned_Metadata_070324.csv

## Supplementary Tables 
Supp Table 1. Created in Word. Source data is available in Cleaned_Metadata_070324.csv
Supp Table 2. Created in Word. Source data is available in Cleaned_Metadata_070324.csv
Supp Table 3. Created in Word. Source data is available in Supplementary_Data3_CytokineandOutcomeData.csv

## Supplementary Data 
Supp Data 1: generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R, saved as Supplementary_Data1_Sepsis_DEGenes.csv 

Supp Data 2: generated by inputting the DE genes from Supplementary_Data1_Sepsis_DEGenes.csv into WebGestalt, saved as Supplementary_Data2_Sepsis_ORAPathways.csv.

Supp Data 3: generated by the code Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R, saved as Supplementary_Data3_DEGenes_pathdominance_highlow.csv. 

Supp Data 4: generated by inputting the DE genes from Supplementary_Data3_DEGenes_pathdominance_highlow.csv into WebGestalt, saved as Supplementary_Data4_ORAPathways_SepsisPatients_bypathogendominance_highlow.csv. 

Supp Data 5: generated by the code Script4_3.preds_study_v2.R, saved as Supplementary_Data5_microbial_input_transcriptomic_classifier. 

Supp Data 6: generated by the code Script2_de_external_cohorts-2, saved as Supplementary_Data6_GAINS_DEGenes.csv. 

Supp Data 7: generated by inputting the DE genes from Supplementary_Data6_GAINS_DEGene.csv into WebGestalt, saved as Supplementary_Data7_GAINS_DEGene_ORAPathways.csv. 

Supp Data 8: generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R, saved as Supplementary_Data8_ICUControls_DEGenes.csv. 

Supp Data 9: generated by inputting the DE genes from Supplementary_Data13_ICUControls_DEGenes.csv into WebGestalt, saved as Supplementary_Data9_ICUControls_ORAPathways.csv.

Supp Data 10: generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R, saved as Supplementary_Data10_ICUControls_DEGenes_CardiacArrestasCovariate.csv. 

Supp Data 11: generated by inputting the DE genes from Supplementary_Data10_ICUControls_DEGenes_CardiacArrestasCovariate.csv into WebGestalt, saved as Supplementary_Data11_ICUControls_ArrestAsCofactor_ORAPathways.csv.

Supp Data 12: List of additional metagenomic contaminants excluded from analysis. 

## Source Data File 
###All Source Data files 1-8 combined into single xlsx file; "source_data.xlsx"
Source Data 1: generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R, saved as Source_Data1_SpecificGenesSepsis.csv 
Source Data 2: generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R, saved as Source_Data2_CytokineandOutcomeData.csv 
Source Data 3: generated by the code Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R, saved as Source_Data3_sepsis_microbialmass_and_pathdominance.csv. 
Source Data 4: generated by the code Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R, saved as Source_Data4_mortalitypresenceabsenceofpathogens.csv. 
Source Data 5: generated by the code Script4_3.preds_study_v2.R, saved as Source_Data5_AUCforDifferentClassifiers.csv. 
Source Data 6: generated by the code Script4_3.preds_study_v2.R, saved as Source_Data6_KaplanMeier_hostpathogen.csv. 
Source Data 7: generated by the code Script1_EARLISurvivalvsMortality_HostAnalyses.R, saved as  Source_Data7_DEGene_Overlap.csv 
Source Data 8: generated by the code Script3_EARLISurvivalvsMortality_MicrobialAnalyses.R, saved as Source_Data8_ICUControls_microbialmass_and_pathdominance.csv. 


## Required hardware and software

The codes were run on a Mac laptop. The required R packages (see below for more details) could be installed with the commands `install.package()` and `BiocManager::install()`. Each package can take up to 1 minute to install. Please refer to each package's website for more information on the installation.

```
R version 4.3.1 (2023-06-16)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Sonoma 14.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Los_Angeles
tzcode source: internal

attached base packages:
 [1] grid      parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pROC_1.18.4                    rpart_4.1.19                   SepstratifieR_1.0             
 [4] VennDiagram_1.7.3              futile.logger_1.4.3            hablar_0.3.2                  
 [7] magrittr_2.0.3                 mia_1.8.0                      TreeSummarizedExperiment_2.8.0
[10] Biostrings_2.68.1              XVector_0.40.0                 lubridate_1.9.3               
[13] forcats_1.0.0                  purrr_1.0.2                    readr_2.1.4                   
[16] tibble_3.2.1                   tidyverse_2.0.0                MAST_1.25.1                   
[19] SingleCellExperiment_1.22.0    org.Hs.eg.db_3.17.0            AnnotationDbi_1.62.2          
[22] pheatmap_1.0.12                factoextra_1.0.7               readxl_1.4.3                  
[25] ggstats_0.6.0                  mclust_6.0.0                   clv_0.3-2.4                   
[28] class_7.3-22                   cluster_2.1.4                  doRNG_1.8.6                   
[31] rngtools_1.5.2                 foreach_1.5.2                  microViz_0.12.1               
[34] ANCOMBC_2.2.2                  ggstance_0.3.6                 nnet_7.3-19                   
[37] scales_1.3.0                   singscore_1.20.0               vegan_2.6-4                   
[40] lattice_0.21-8                 permute_0.9-7                  infotheo_1.2.0.1              
[43] rstatix_0.7.2                  binom_1.1-1.1                  patchwork_1.2.0               
[46] e1071_1.7-13                   biomaRt_2.56.1                 DescTools_0.99.50             
[49] fgsea_1.26.0                   ggrepel_0.9.4                  plyr_1.8.8                    
[52] survival_3.5-7                 survminer_0.4.9                phyloseq_1.44.0               
[55] data.table_1.14.10             tidyr_1.3.0                    DESeq2_1.40.2                 
[58] edgeR_3.42.4                   limma_3.56.2                   msigdbr_7.5.1                 
[61] bestNormalize_1.9.1            table1_1.4.3                   ggsankey_0.0.99999            
[64] dplyr_1.1.4                    ggfortify_0.4.16               poolr_1.1-1                   
[67] mice_3.16.0                    corrplot_0.92                  tidyLPA_1.1.0                 
[70] haven_2.5.4                    lmerTest_3.1-3                 ggpubr_0.6.0                  
[73] ggplot2_3.5.0                  simr_1.0.7                     lme4_1.1-34                   
[76] Matrix_1.6-1.1                 pwr_1.3-0                      reshape2_1.4.4                
[79] stringr_1.5.1                  SeuratDisk_0.0.0.9021          reticulate_1.37.0             
[82] import_1.3.2                   scp_1.9.0                      QFeatures_1.10.0              
[85] MultiAssayExperiment_1.26.0    SummarizedExperiment_1.30.2    Biobase_2.60.0                
[88] GenomicRanges_1.52.0           GenomeInfoDb_1.36.3            IRanges_2.34.1                
[91] S4Vectors_0.38.1               BiocGenerics_0.46.0            MatrixGenerics_1.12.3         
[94] matrixStats_1.0.0              SeuratObject_4.1.3             Seurat_4.3.0.1                
[97] anndata_0.7.5.6               
```

### Python scripts
```
Python 3.12.3 with the following libraries:

Package          Version
---------------- -----------
joblib           1.4.2
matplotlib       3.9.1
numpy            1.26.4
pandas           2.2.2
rpy2             3.5.16
seaborn          0.13.2
scikit-learn     1.4.2
contourpy        1.2.1
cycler           0.12.1
fonttools        4.53.1
kiwisolver       1.4.5
packaging        24.1
pillow           10.4.0
pyparsing        3.1.2
python-dateutil  2.9.0.post0
six              1.16.0
pytz             2024.1
tzdata           2024.1
cffi             1.16.0
pycparser        2.22
Jinja2           3.1.4
MarkupSafe       2.1.5
tzlocal          5.2
scipy            1.14.0
threadpoolctl    3.5.0
```




