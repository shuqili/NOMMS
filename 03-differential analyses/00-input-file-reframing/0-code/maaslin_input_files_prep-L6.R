library(tidyverse)

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/07-taxa/05-count-tables")
list.files()
read_tsv("../../metadata-CL04-NOMMS-v5.txt")
read_tsv("filtered-genus-count-table-mapped-w-human-oral-weighted-silva138_1.txt")#, skip = 1) %>% colnames()


metadata <-
  read_tsv("../../metadata-CL04-NOMMS-v5.txt") %>%
  select(c(sample_name,"SubjectID", "Anatomical_Location", "DiseaseStatusOracle_abbrev"), everything()) %>%
  rename(sample_id = sample_name) %>%
  mutate(DiseaseStatusOracle_abbrev= factor(DiseaseStatusOracle_abbrev, 
                                            levels=c("HC", "RRMS", "Prog MS"),
                                            labels=c("HC", "RRMS", "Prog MS"))) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, 
                                      levels=c("N", "G", "T"),
                                      labels=c("nasal", "gingival", "oropharyngeal"))) %>%
  mutate(sex = factor(sex, 
                      levels = c("M","F"),
                      labels = c("Males","Females"))) 

# read_tsv("../02-mapped_taxonomy/taxonomy-mapped-w-human-oral-weighted-silva138_1.tsv")
# taxonomy-mapped-w-human-oral-weighted-silva138_1.tsv

taxonomy <-
  read_tsv("../02-mapped_taxonomy/taxonomy-mapped-w-human-oral-weighted-silva138_1.tsv") %>%
  rename_all(tolower) %>%
  rename(feature_ID = `feature id`) %>%
  separate(taxon, sep = "; ", fill = "warn", extra = "warn",
           into = c("kingdom", "phylum", "class", "order", 
                    "family", "genus", "species")) %>%
  mutate(order = str_replace_na(order, "o__unclassified")) %>% #summarize(count = sum(is.na(order)))
  mutate(family = str_replace_na(family, "f__unclassified")) %>% #summarize(count = sum(is.na(family)))
  mutate(genus = str_replace_na(genus, "g_unclassified")) %>%
  mutate(species = str_replace_na(species, "s__unclassified")) %>%
  mutate(ASV = str_c("ASV__",feature_ID)) 

genus_counts <-
  read_tsv("filtered-genus-count-table-mapped-w-human-oral-weighted-silva138_1.txt") %>%
  rename(sample_id = index) %>% #colnames()
  select(c(1:137)) %>%
  pivot_longer(-sample_id, names_to = "genus", values_to = "count")

genus_rel_abund <-
  genus_counts %>%  
  group_by(sample_id) %>%
  mutate(rel_abund = 100 * count/sum(count)) %>%
  ungroup() %>%
  select(-count) 

n_sample <-
  inner_join(metadata, genus_rel_abund, by = c("sample_id")) %>% 
  select(Anatomical_Location, DiseaseStatusOracle_abbrev, sample_id, excluded) %>%
  filter(excluded == "no") %>%
  select(-excluded) %>%
  distinct() %>%
  unite("Groups", c("Anatomical_Location", "DiseaseStatusOracle_abbrev"), remove = T) %>%
  group_by(Groups) %>%
  summarise(n_sample_per_group = n(), .groups = "drop")

genus_prevalence_remove <-
  inner_join(metadata, genus_counts, by = c("sample_id")) %>%
  filter(str_detect(excluded, "no")) %>%
  select(-excluded) %>%
  select(sample_id, Anatomical_Location, DiseaseStatusOracle_abbrev, genus, count) %>%
  unite("Groups", c("Anatomical_Location", "DiseaseStatusOracle_abbrev"), remove = T) %>%
  group_by(Groups, sample_id, genus) %>%
  mutate(presence = if_else(count > 0, 1, 0)) %>% 
  ungroup() %>%
  group_by(Groups, genus) %>%
  mutate(n_sample_w_feature_present_per_group = sum(presence)) %>%
  ungroup() %>%
  left_join(., n_sample, by = "Groups") %>%
  group_by(Groups, genus) %>%
  reframe(prevalence = 100 * (n_sample_w_feature_present_per_group/n_sample_per_group)) %>%
  group_by(genus) %>%
  reframe(remove = max(prevalence) < 10)

maaslin_input_l6 <-
  inner_join(genus_rel_abund, genus_prevalence_remove, by = "genus") %>%
  filter(remove=="FALSE") %>%
  select(genus, sample_id, rel_abund) %>%
  pivot_wider(names_from = sample_id, values_from = rel_abund, values_fill = 0) 

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin")
write_tsv(x=maaslin_input_l6,
          file = "input_genus_rel_abund.tsv")
