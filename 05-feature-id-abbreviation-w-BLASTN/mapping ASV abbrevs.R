library(tidyverse)
library(svglite)
library(seqinr)
# getwd()

setwd("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin")

list.files()

# Set critical variables (uniform for all)
p = 0.05

# Read in results with original ASV
# DS (disease status)
ds.vs.HC.nasal <-
  read_tsv("03-output/01-nostril/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(pval < p) %>% 
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs HC")) %>%
  mutate(AL = "nasal")

ds.vs.HC.gingival <-
  read_tsv("03-output/02-gums/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs HC")) %>%
  mutate(AL = "gingival")

ds.vs.HC.oro <-
  read_tsv("03-output/03-throat/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs HC")) %>%
  mutate(AL = "oropharyngeal")

ds.vs.HC <- bind_rows(ds.vs.HC.nasal, ds.vs.HC.gingival, ds.vs.HC.oro)

ds.vs.RR.nasal <-
  read_tsv("03-output/01-nostril/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/PMS_vs_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(pval < p) %>% 
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs RR")) %>%
  mutate(AL = "nasal")

ds.vs.RR.gingival <-
  read_tsv("03-output/02-gums/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/PMS_vs_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs RR")) %>%
  mutate(AL = "gingival")

ds.vs.RR.oro <-
  read_tsv("03-output/03-throat/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/PMS_vs_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs RR")) %>%
  mutate(AL = "oropharyngeal")

ds.vs.RR <- bind_rows(ds.vs.RR.nasal, ds.vs.RR.gingival, ds.vs.RR.oro)

ds.sig.res <- bind_rows(ds.vs.HC, ds.vs.RR)

# DMT nasal

dmt.vs.untreated <-
  read_tsv("03-output/01-nostril/ASVlevel/DMT-N_over_5 vs Untreated/DMT-model4/cplmfit-no-transformation-relab-asv-to-DMT/all_results.tsv") %>%
  filter(pval < p) %>% 
  filter(metadata == "DMT_v2") %>%
  mutate(value = str_replace(value, "^(.*)$", "\\1 vs Untreated")) %>%
  mutate(AL = "nasal")


# EDSS change over 2 years.
edss.delta.nasal <-
read_tsv("03-output/01-nostril/ASVlevel/edss_change + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "edss_change") %>%
  mutate(AL = "nasal")

edss.delta.gingival <-
  read_tsv("03-output/02-gums/ASVlevel/edss_change + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "edss_change") %>%
  mutate(AL = "gingival")

edss.delta.oro <-
  read_tsv("03-output/03-throat/ASVlevel/edss_change + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(pval < p) %>%
  filter(metadata == "edss_change") %>%
  mutate(AL = "oro")

edss.delta <- bind_rows(edss.delta.nasal,edss.delta.gingival, edss.delta.oro)

# combine all ASV ID in one form

all.sig.ASVs <- bind_rows(ds.sig.res, dmt.vs.untreated, edss.delta) %>%
  arrange(desc(coef)) %>%
  select(AL, metadata,feature)

asv.abbrev.list <-
  all.sig.ASVs %>%
  mutate(ASV.old = str_replace(feature, "^(.*)ASV__", "")) %>%
  distinct(ASV.old) %>%
  mutate(ASV.abbrev = seq.int(1, 92, 1)) %>%
  mutate(ASV.abbrev= as.character(ASV.abbrev)) %>%
  mutate(ASV.abbrev = str_replace(ASV.abbrev,"^(\\d)$","0\\1")) %>%
  mutate(ASV.abbrev = str_replace(ASV.abbrev,"^(\\d+)$","ASV_\\1")) 

# Export asv.abbrev.list
write_tsv(x = asv.abbrev.list, file = "ASV.abbreviation.list.tsv")

# Expand ASV abbreviation list by adding female/male specific taxonomic alteration results based on two-year edss change.
# 2025-03-31 eBioMedicine R1

# read-in new data
# Read in data to plot
female.nasal <- read_tsv("03-output/01-nostril/ASVlevel/female-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  mutate(anatomical_site = "nasal")
female.gin <- read_tsv("03-output/02-gums/ASVlevel/female-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")  %>%
  mutate(anatomical_site = "gingival")
female.oro <- read_tsv("03-output/03-throat/ASVlevel/female-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  mutate(anatomical_site = "oropharyngeal")

male.nasal <- read_tsv("03-output/01-nostril/ASVlevel/male-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  mutate(anatomical_site = "nasal")
male.gin <- read_tsv("03-output/02-gums/ASVlevel/males-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")  %>%
  mutate(anatomical_site = "gingival")
male.oro <- read_tsv("03-output/03-throat/ASVlevel/male-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")   %>%
  mutate(anatomical_site = "oropharyngeal")

females.res <-
  bind_rows(female.nasal, female.gin, female.oro) %>%
  mutate(sex = "female")

males.res <-
  bind_rows(male.nasal, male.gin, male.oro) %>%
  mutate(sex = "male")


bind_rows(females.res, males.res) %>%
  filter(metadata == "edss_change") %>%
  filter(pval < p) %>%
  select(feature) %>%
  distinct() %>%
  mutate(feature = str_remove(feature, ".*ASV__")) %>%
  mutate(feature.copy = feature) %>%
  full_join(., asv.abbrev.list, by = c("feature" = "ASV.old")) %>% 
  arrange(ASV.abbrev) %>%
  rename(added.feature = feature.copy) %>%
  mutate( ASV_ = "ASV_", ASV.number = c(1:109)) %>%
  mutate(ASV.number= as.character(ASV.number)) %>%
  mutate(ASV.number = str_replace(ASV.number,"^(\\d)$","0\\1")) %>%
  mutate(ASV.number = str_replace(ASV.number,"^(\\d+)$","ASV_\\1"))  %>% 
  select(feature, ASV.number) %>%
  write_tsv("ASV.abbreviation.list.v2.tsv")

  
getwd()







# taxonomy <- read_tsv("../07-taxa/02-mapped_taxonomy/taxonomy-mapped-w-human-oral-weighted-silva138_1.tsv") %>%
#   rename( ASV = `Feature ID`)
# 
# left_join(asv.abbrev.list, taxonomy, by = c("ASV.old"="ASV") )  %>%
#   rename(taxon.oral.silva.138.1 = Taxon) %>%
#   # filter(ASV.abbrev == "ASV_28")
#   # filter(str_detect())
#   filter(str_detect(taxon.oral.silva.138.1, "Rothia"))


# read_tsv("../03-features/oral-silva138_1/sequences.fasta")
# read.fasta("../03-features/oral-silva138_1/sequences.fasta")
# smallFastaFile <- system.file("sequences/smallAA.fasta", package = "seqinr")
# mySmallProtein <- read.fasta(file = smallFastaFile, as.string = TRUE, seqtype = "AA")[[1]]
