
# install.packages("devtools")
# library("devtools")
# install_github("biobakery/Maaslin2")

# library(maaslin2)
# call dev installed pkg by --

# Maaslin2::Maaslin2()


getwd()


# library(Maaslin2)
library(tidyverse)

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/")

# 
# read_tsv("metadata-CL04-NOMMS-v5.txt") %>%
#   filter(Anatomical_Location=="G") %>%
#   rename(sample_ID = `#Sample ID`) %>%
#   write_tsv("09-MaAsLin/01-input/metadata-v5-gums.tsv")

df_input_metadata = read.table(file             = "09-MaAsLin/01-input/metadata-v5-nostril.tsv", 
                               header           = TRUE, 
                               sep              = "\t", 
                               row.names        = 1,
                               stringsAsFactors = FALSE) %>%
  mutate(smoking_YN = if_else(smoking_YN == 0, "no", "yes"))
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

# fit_data_1 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/lmfit-ast-transformed-relab-asv-to-diseasestat",
#                       input_data     = df_input_data, 
#                       input_metadata = df_input_metadata, 
#                       
#                       analysis_method = "LM",
#                       normalization  = "NONE",
#                       transform = "AST",
#                       min_prevalence = 0.1,
#                       standardize = TRUE,
#                       cores = 4,
#                       fixed_effects  = c("age","sex","bmi",
#                                          "smoking"))
# 
# 
# fit_data_2 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/lmfit-log-transformed-relab-asv-to-diseasestat", 
#                       input_data     = df_input_data, 
#                       input_metadata = df_input_metadata, 
#                       analysis_method = "LM",
#                       normalization  = "NONE",
#                       transform = "LOG",
#                       min_prevalence = 0.1,
#                       standardize = TRUE,
#                       cores = 4,
#                       fixed_effects  = c("age","sex","bmi",
#                                          "smoking"))
# 
# fit_data_3 = Maaslin2(output = "09-MaAsLin/03-output/01-nostril/lmfit-logit-transformed-relab-asv-to-diseasestat", 
#                       input_data     = df_input_data, 
#                       input_metadata = df_input_metadata, 
#                       analysis_method = "LM",
#                       normalization  = "NONE",
#                       transform = "LOGIT",
#                       min_prevalence = 0.1,
#                       standardize = TRUE,
#                       cores = 4,
#                       fixed_effects  = c("age","sex","bmi",
#                                          "smoking"))
dir.create("09-MaAsLin/03-output/01-nostril/l6taxa/Smoking_YN + [age + sex + bmi ]/")
fit_data_4 = Maaslin2::Maaslin2(output = "09-MaAsLin/03-output/01-nostril/l6taxa/Smoking_YN + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/", 
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
                      fixed_effects  = c("age","sex","bmi",
                                         "smoking_YN"),
                      reference = c("smoking_YN","no"))

fixed_effects  = c("DiseaseStatusOracle_abbrev",
                   "age","sex","bmi",
                   "smoking",
                   "dental_procedure_simplified",
                   "oral_disease_simplified",
                   "oral_hygiene_cumsum_score"),
reference      = c("DiseaseStatusOracle_abbrev,RRMS;smoking,NONE;dental_procedure_simplified,none;oral_disease_simplified,none"))




getwd()

































--------------------------------------------------------------------
  sessionInfo()
# 
# R version 4.4.2 (2024-10-31)
# Platform: x86_64-apple-darwin20
# Running under: macOS Monterey 12.7.5
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] devtools_2.4.5      usethis_3.0.0       BiocManager_1.30.25
# 
# loaded via a namespace (and not attached):
#   [1] miniUI_0.1.1.1    compiler_4.4.2    promises_1.3.0    Rcpp_1.0.13-1     callr_3.7.6      
# [6] later_1.3.2       fastmap_1.2.0     mime_0.12         R6_2.5.1          curl_6.0.1       
# [11] htmlwidgets_1.6.4 desc_1.4.3        profvis_0.4.0     shiny_1.9.1       rlang_1.1.4      
# [16] cachem_1.1.0      httpuv_1.6.15     fs_1.6.5          pkgload_1.4.0     memoise_2.0.1    
# [21] cli_3.6.3         magrittr_2.0.3    ps_1.8.1          processx_3.8.4    digest_0.6.37    
# [26] rstudioapi_0.17.1 xtable_1.8-4      remotes_2.5.0     lifecycle_1.0.4   vctrs_0.6.5      
# [31] glue_1.8.0        urlchecker_1.0.1  sessioninfo_1.2.2 pkgbuild_1.4.5    purrr_1.0.2      
# [36] tools_4.4.2       ellipsis_0.3.2    htmltools_0.5.8.1
