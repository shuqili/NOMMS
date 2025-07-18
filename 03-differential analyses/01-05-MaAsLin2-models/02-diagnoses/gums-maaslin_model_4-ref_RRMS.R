library(Maaslin2)
library(tidyverse)

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/")


# read_tsv("metadata-CL04-NOMMS-v5.txt") %>%
#   filter(Anatomical_Location=="G") %>%
#   rename(sample_ID = `#Sample ID`) %>%
#   write_tsv("09-MaAsLin/01-input/metadata-v5-gums.tsv")
list.files("09-MaAsLin/01-input/")

df_input_metadata = read.table(file             = "09-MaAsLin/01-input/metadata-v5-gums.tsv", 
                               header           = TRUE, 
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE)
df_input_metadata[1:5, ]

df_input_data = read.table(file             = "09-MaAsLin/01-input/input_genus_rel_abund.tsv",
                           header           = TRUE,
                           sep              = "\t", 
                           row.names        = 1,
                           stringsAsFactors = FALSE)

df_input_data[1:5, 1:5]

# ASV relative abundance as data has already been normalized.
# Fit 1: LM-AST
# Fit 2: LM-LOG
# Fit 3: LM-LOGIT
# Fit 4: CPLM-NONE
getwd()
dir.create("09-MaAsLin/03-output/02-gums/l6taxa/diagnoses_ref_RRMS/")

fit_data_1 = Maaslin2(output = "09-MaAsLin/03-output/02-gums/l6taxa/diagnoses_ref_RRMS/lmfit-ast-transformed-relab-asv-to-diseasestat",
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                       
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "AST",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c("DiseaseStatusOracle_abbrev",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DiseaseStatusOracle_abbrev,RRMS;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))


fit_data_2 = Maaslin2(output = "09-MaAsLin/03-output/02-gums/l6taxa/diagnoses_ref_RRMS/lmfit-log-transformed-relab-asv-to-diseasestat", 
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOG",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c("DiseaseStatusOracle_abbrev",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DiseaseStatusOracle_abbrev,RRMS;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

fit_data_3 = Maaslin2(output = "09-MaAsLin/03-output/02-gums/l6taxa/diagnoses_ref_RRMS/lmfit-logit-transformed-relab-asv-to-diseasestat", 
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "LM",
                      normalization  = "NONE",
                      transform = "LOGIT",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c("DiseaseStatusOracle_abbrev",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DiseaseStatusOracle_abbrev,RRMS;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))

fit_data_4 = Maaslin2(output = "09-MaAsLin/03-output/02-gums/l6taxa/diagnoses_ref_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat", 
                      input_data     = df_input_data, 
                      input_metadata = df_input_metadata, 
                      analysis_method = "CPLM",
                      normalization  = "NONE",
                      transform = "NONE",
                      min_prevalence = 0.1,
                      standardize = TRUE,
                      cores = 4,
                      plot_heatmap = F,
                      plot_scatter = F,
                      fixed_effects  = c("DiseaseStatusOracle_abbrev",
                                         "age","sex","bmi",
                                         "smoking",
                                         "dental_procedure_simplified",
                                         "oral_disease_simplified",
                                         "oral_hygiene_cumsum_score"),
                      reference      = c("DiseaseStatusOracle_abbrev,RRMS;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))
