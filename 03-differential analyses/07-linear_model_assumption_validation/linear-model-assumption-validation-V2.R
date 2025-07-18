library(tidyverse)
library(ggplot2)
library(ggtext)
# https://stats.stackexchange.com/questions/99052/residuals-in-poisson-regression
# https://rpubs.com/Julian_Sampedro/1047952
getwd()
setwd("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")
# use gingival data to illustrate
setwd("03-output/02-gums/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/")
list.files()

#1. Import the residual.rds and fitted.rds from CPLM model and the log.normal model.
gingival.cplm.residuals <- file.choose(new = "cplmfit-no-transformation-relab-asv-to-diseasestat/residuals.rds")
gingival.cplm.fitted <- file.choose(new = "cplmfit-no-transformation-relab-asv-to-diseasestat/fitted.rds")
gingival.log.normal.residuals <- file.choose("../lmfit-log-transformed-relab-asv-to-diseasestat/residuals.rds")
gingival.log.normal.fitted <- file.choose("../lmfit-log-transformed-relab-asv-to-diseasestat/fitted.rds")

#2. Read-in .rds as matrix.
cplm.residuals.rds <- read_rds(gingival.cplm.residuals)
cplm.fitted.rds <- read_rds(gingival.cplm.fitted) 
log.normal.residuals.rds <- read_rds(gingival.log.normal.residuals)
log.normal.fitted.rds <- read_rds(gingival.log.normal.fitted)

#3. Plot residual against predicted value for individual taxon. rep1 = ASV_08, rep2 = ASV_29.
# representative feature 1 Fusobacterium ASV_08
##cplm
g.fuso.residual.cplm <- cplm.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>% 
  filter(str_detect(feature, "g__Fusobacterium.s__uncultured_bacterium.ASV__bf7ec0c65c8f126e7b866def613f7145")) %>%
  as.matrix()
g.fuso.predicted.cplm <- cplm.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "g__Fusobacterium.s__uncultured_bacterium.ASV__bf7ec0c65c8f126e7b866def613f7145")) %>%
  as.matrix()
plot(y=g.fuso.residual.cplm, x=g.fuso.predicted.cplm,
     ylab = "CPLM residuals",xlab = "CPLM fitted values",
     main = substitute(paste(italic('Fusobacterium ASV_08'))))


##log.normal
g.fuso.residual.log <-log.normal.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "g__Fusobacterium.s__uncultured_bacterium.ASV__bf7ec0c65c8f126e7b866def613f7145")) %>%
  as.matrix()
g.fuso.predicted.log <- log.normal.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "g__Fusobacterium.s__uncultured_bacterium.ASV__bf7ec0c65c8f126e7b866def613f7145")) %>%
  as.matrix()
plot(y = g.fuso.residual.log, x = g.fuso.predicted.log,
     ylab = "Log-transformed GLM residuals",xlab = "Log-transformed GLM fitted values",
     main = substitute(paste(italic('Fusobacterium ASV_08'))))
ggsave("Fusobaacterium ASV_08 log-transformed GLM residual plot.pdf")
getwd()
# representative feature 2 C. showae ASV_29
##cplm
g.showae.residual.cplm <- cplm.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(.,"feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
g.showae.predicted.cplm <-cplm.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
plot(y = g.showae.residual.cplm, x=g.showae.predicted.cplm)

##log.normal
g.showae.residual.log <- log.normal.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
g.showae.predicted.log <- log.normal.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
plot(y=g.showae.residual.log,x=g.showae.predicted.log)  







#4. Plot Pearson residuals against fitted values for Poisson model
# representative feature 1 Fusobacterium ASV_08
##cplm
g.fuso.residual.cplm.pearson <- cplm.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>% 
  filter(str_detect(feature, "g__Fusobacterium.s__uncultured_bacterium.ASV__bf7ec0c65c8f126e7b866def613f7145")) %>% 
  rowwise() %>%
  mutate(mean = mean(c_across(c(2:106)))) %>% 
  mutate(sd = sd(c_across(c(2:106)))) %>%
  ungroup() %>%
  mutate(across(c(2:106), ~ .x/sd - mean/sd)) %>%
  select(-mean, -sd) %>%
  as.matrix() 
g.fuso.predicted.cplm <- cplm.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "g__Fusobacterium.s__uncultured_bacterium.ASV__bf7ec0c65c8f126e7b866def613f7145")) %>%
  as.matrix()
plot(y=g.fuso.residual.cplm.pearson, x=g.fuso.predicted.cplm)
# representative feature 2 C. showae ASV_29
##cplm
g.showae.residual.cplm.pearson <-
  cplm.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(.,"feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  rowwise() %>%
  mutate(mean = mean(c_across(c(2:106)))) %>%
  mutate(sd = sd(c_across(c(2:106)))) %>%
  ungroup() %>%
  mutate(across(c(2:106), ~ .x/sd - mean/sd)) %>%
  select(-mean, -sd) %>%
  as.matrix()
  
g.showae.predicted.cplm <- cplm.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
plot(y = g.showae.residual.cplm.pearson, x=g.showae.predicted.cplm)

##log.normal
g.showae.residual.log <- log.normal.residuals.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
g.showae.predicted.log <- log.normal.fitted.rds %>%
  as.data.frame() %>%
  rownames_to_column(., "feature") %>%
  filter(str_detect(feature, "s__Campylobacter_showae.ASV__b21bbc630376fdeea3d592a807d5f0d8")) %>%
  as.matrix()
plot(y=g.showae.residual.log,x=g.showae.predicted.log)  

#5. Plot residual histograms - pearson residuals of CPLM model
# rep1 ASV_08 Fuso
g.fuso.residual.cplm.pearson %>%
  as.t
  ggplot(aes(x))








