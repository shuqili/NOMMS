library(tidyverse)
library(readxl)

getwd()
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/03b-EbioMedicine-Revise/Adjusted_Figures_Tables/")
list.files()
#1 oral disease, dental procedures, and oral hygiene related practice.
a<-read_excel("../../02-NOMMS_SRA_submission/NOMMS.MIMARKS.survey.human-associated.6.0.xlsx", skip = 11) %>%
  rename(sample_id = `*sample_name`) %>%
  select(sample_id, anatomical_location,host_subject_id, nose_throat_disord) %>%
  select(host_subject_id, nose_throat_disord) %>%
  # replace_na(.,"")
  mutate(Procedure.bridges = if_else(str_detect(nose_throat_disord, "bridges"), "1","0")) %>%
  mutate(Procedure.crowns = if_else(str_detect(nose_throat_disord, "crowns"), "1","0")) %>%
  mutate(Procedure.dentures = if_else(str_detect(nose_throat_disord, "dentures"), "1","0")) %>%
  mutate(Procedure.extractions = if_else(str_detect(nose_throat_disord, "extractions"), "1","0")) %>%
  mutate(Procedure.fillings = if_else(str_detect(nose_throat_disord, "fillings"), "1","0")) %>%
  mutate(Procedure.gum_surgery = if_else(str_detect(nose_throat_disord, "gum_surgery"), "1","0")) %>%
  mutate(Procedure.oral_surgery = if_else(str_detect(nose_throat_disord, "oral_surgery"), "1","0")) %>%
  mutate(Procedure.root_canals = if_else(str_detect(nose_throat_disord, "root_canals"), "1","0")) %>%
  mutate(Procedure.sealants = if_else(str_detect(nose_throat_disord, "sealants"), "1","0")) %>%
  mutate(Procedure.teeth_whitening = if_else(str_detect(nose_throat_disord, "teeth_whitening"), "1","0")) %>%
  mutate(Procedure.veneers = if_else(str_detect(nose_throat_disord, "veneers"), "1","0")) %>%
  mutate(Procedure.apicoectomy = if_else(str_detect(nose_throat_disord, "apicoectomy"), "1","0")) %>%
  mutate(Procedure.implant = if_else(str_detect(nose_throat_disord, "implant"), "1","0")) %>%
  mutate(Disease.dental_caries = if_else(str_detect(nose_throat_disord, "dental_caries"), "1","0")) %>%
  mutate(Disease.gingivitis = if_else(str_detect(nose_throat_disord, "gingivitis"), "1","0")) %>%
  mutate(Disease.periodontal_disease = if_else(str_detect(nose_throat_disord, "periodontal_disease"), "1","0")) %>%
  mutate(Disease.oral_cancer = if_else(str_detect(nose_throat_disord, "oral_cancer"), "1","0")) %>%
  mutate(Disease.dental_abscess = if_else(str_detect(nose_throat_disord, "dental_abscess"), "1","0")) %>%
  mutate(Disease.reoccurring_ulcers = if_else(str_detect(nose_throat_disord, "reoccurring_ulcers"), "1","0")) %>%
  distinct()
  
b <- read_tsv("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/metadata-CL04-NOMMS-v4.txt") %>%
   select(SubjectID, DiseaseStatusOracle_abbrev, brush_freq, floss_freq, mouthwash_freq, dentist_freq) %>%
   filter(!is.na(DiseaseStatusOracle_abbrev)) %>%
   distinct()
  
oral.health <- full_join(a, b, by = c("host_subject_id"="SubjectID")) %>% 
  select(-nose_throat_disord) 

oral.health %>% colnames()


table.2.d.p <-
  oral.health %>%
  select(host_subject_id, DiseaseStatusOracle_abbrev, everything()) %>%
  select(-brush_freq,-floss_freq, -mouthwash_freq,-dentist_freq) %>%
  pivot_longer(-c(1:2), names_to = "variable", values_to = "value") %>%
  mutate(value = if_else(value == "1", 1, 0))  %>%
  group_by(DiseaseStatusOracle_abbrev) %>%
  mutate(total.number = n_distinct(host_subject_id)) %>%
  ungroup() %>% 
  replace_na(list(value = 0)) %>%
  group_by(DiseaseStatusOracle_abbrev, variable) %>%
  mutate(yes.number = sum(value))  %>%
  ungroup() %>%
  select(-host_subject_id, -value) %>%
  distinct() %>%
  select(DiseaseStatusOracle_abbrev, total.number, everything()) %>%
  unite("group", DiseaseStatusOracle_abbrev:total.number) %>%
  pivot_wider(names_from = group, values_from = yes.number)
  
brush.freq.count.table <-
oral.health %>% #colnames()
  select(host_subject_id,DiseaseStatusOracle_abbrev, brush_freq) %>%
  replace_na(list(brush_freq="no.data")) %>%
  # mutate(brush_freq = factor(brush_freq, levels=c("Never","Once a month","Once a week",
  #                                                 "Once a day","Twice a day",
  #                                                 "More than twice a day", "no.data"))) %>%
  group_by(DiseaseStatusOracle_abbrev, brush_freq) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(-host_subject_id) %>% distinct() %>%
  pivot_wider(names_from = DiseaseStatusOracle_abbrev, values_from = count,
              values_fill = 0) %>%
  add_row(brush_freq="Never", HC = 0, RRMS = 0, `Prog MS` = 0) 
  
brush.freq.count.table <-
  oral.health %>% #colnames()
  select(host_subject_id,DiseaseStatusOracle_abbrev, brush_freq) %>%
  replace_na(list(brush_freq="no.data")) %>%
  group_by(DiseaseStatusOracle_abbrev, brush_freq) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(-host_subject_id) %>% distinct() %>%
  pivot_wider(names_from = DiseaseStatusOracle_abbrev, values_from = count,
              values_fill = 0) %>%
  add_row(brush_freq="Never", HC = 0, RRMS = 0, `Prog MS` = 0) %>%
  mutate(brush_freq = str_replace(brush_freq,"(.*)","brush(\\1)")) %>%
  rename(var = brush_freq)

floss.freq.count.table <-
  oral.health %>% 
  select(host_subject_id,DiseaseStatusOracle_abbrev, floss_freq) %>%
  replace_na(list(floss_freq="no.data")) %>%
  group_by(DiseaseStatusOracle_abbrev, floss_freq) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(-host_subject_id) %>% distinct() %>%
  pivot_wider(names_from = DiseaseStatusOracle_abbrev, values_from = count,
              values_fill = 0) %>%
  mutate(floss_freq = str_replace(floss_freq,"(.*)","floss(\\1)")) %>%
  rename(var = floss_freq)

mouthwash_freq.count.table <-
  oral.health %>% #colnames()
  select(host_subject_id,DiseaseStatusOracle_abbrev, mouthwash_freq) %>%
  replace_na(list(mouthwash_freq="no.data")) %>%
  group_by(DiseaseStatusOracle_abbrev, mouthwash_freq) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(-host_subject_id) %>% distinct() %>%
  pivot_wider(names_from = DiseaseStatusOracle_abbrev, values_from = count,
              values_fill = 0) %>%
  # add_row(mouthwash_freq="Never", HC = 0, RRMS = 0, `Prog MS` = 0) %>%
  mutate(mouthwash_freq = str_replace(mouthwash_freq,"(.*)","mouthwash(\\1)")) %>%
  rename(var = mouthwash_freq)


dentist_freq.count.table <-
  oral.health %>% #colnames()
  select(host_subject_id,DiseaseStatusOracle_abbrev, dentist_freq) %>%
  replace_na(list(dentist_freq="no.data")) %>%
  group_by(DiseaseStatusOracle_abbrev, dentist_freq) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(-host_subject_id) %>% distinct() %>%
  pivot_wider(names_from = DiseaseStatusOracle_abbrev, values_from = count,
              values_fill = 0) %>%
  mutate(dentist_freq = str_replace(dentist_freq,"(.*)","dentist(\\1)")) %>%
  rename(var = dentist_freq)

table2.v3 <-
bind_rows(brush.freq.count.table, floss.freq.count.table,mouthwash_freq.count.table,
          dentist_freq.count.table) %>%
  rename(variable=var, HC_40=HC, RRMS_53=RRMS, `Prog MS_16`=`Prog MS`) %>%
  bind_rows(., table.2.d.p) %>%
  mutate(HC.freq = HC_40/40,  RRMS.freq = RRMS_53/53, ProgMS.freq = `Prog MS_16`/16)
  
# Export
getwd()
write_tsv(x = oral.health, file = "oral-health-related-variables.tsv")
write_tsv(x = table.2, file = "Table2.tsv")
write_tsv(x = table2.v3, file = "Table2_v3.tsv")
# read_xlsx("oral-health-related-variables.xlsx")




#2. Smoking types
table3 <- read_excel("../../02-NOMMS_SRA_submission/NOMMS.MIMARKS.survey.human-associated.6.0.xlsx", skip = 11) %>%
  rename(sample_id = `*sample_name`) %>%
  select(host_subject_id, smoker,host_disease) %>%
  distinct() %>%
  group_by(host_disease, smoker) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(-host_subject_id) %>%
  distinct() %>%
  pivot_wider(names_from = host_disease, values_from = count, values_fill = 0) %>%
  mutate(smoker = factor(smoker, levels = c("Cigarettes_only", "Cigarettes_Tobacco", "Cigarettes_Marijuana",
                                            "Cigars_only",
                                            "Marijuana_only",
                                            "NONE","NA"),
                         labels = c("Cigarettes", "Cigarettes & Tobacco", "Cigarettes & Marijuana",
                                    "Cigars",
                                    "Marijuana",
                                    "Not a smoker","no data"))) %>%
  arrange(smoker) %>%
  mutate(HC = HC/40, RRMS = RRMS/53, `Prog MS` = `Prog MS`/16)

# Export
write_tsv(x=table3, file = "Table3_v1.tsv")


#3. Upper respiratory infection last attack, asthma, and allergy symptoms on the day of visit

table4.b <-
read_excel("../../02-NOMMS_SRA_submission/NOMMS.MIMARKS.survey.human-associated.6.0.xlsx", skip = 11) %>%
  rename(sample_id = `*sample_name`) %>% # colnames()
  select(sample_id, host_subject_id, host_disease, 
         allergy_symptoms_today, upper_respiratory_infection_last_attack,	asthma_last_attack) %>% 
  mutate(allergy_symptoms_today = if_else(allergy_symptoms_today==1,"Yes","No"))  %>%
  mutate(upper_respiratory_infection_last_attack = case_when(upper_respiratory_infection_last_attack == 0 ~ "Unsure or Never",
                                                             upper_respiratory_infection_last_attack == 1 ~ "Over a year ago",
                                                             upper_respiratory_infection_last_attack == 2 ~ "Within the last year",
                                                             upper_respiratory_infection_last_attack == 3 ~ "Within the last month",
                                                             upper_respiratory_infection_last_attack == 4 ~ "Within the last week"))  %>%
  mutate(asthma_last_attack = case_when(asthma_last_attack == 0 ~ "Unsure or Never",
                                        asthma_last_attack == 1 ~ "Over a year ago",
                                        asthma_last_attack == 2 ~ "Within the last year",
                                        asthma_last_attack == 3 ~ "Within the last month",
                                        asthma_last_attack == 4 ~ "Within the last week"))  %>%
  pivot_longer(-c(1:3), names_to = "var", values_to = "value") %>%
  select(-sample_id) %>%
  distinct() %>%
  unite("var_value", var:value,sep = ":", remove = T) %>%
  group_by(host_disease, var_value) %>%
  mutate(count = n_distinct(host_subject_id)) %>%
  ungroup() %>%
  select(host_disease, var_value, count) %>%
  distinct()%>%
  pivot_wider(names_from = host_disease, values_from = count, values_fill = 0) %>%
  arrange(var_value) %>%
  separate(var_value, into = c("var","value"), sep = ":") %>%
  filter(value != "No") %>%
  add_row(var = "asthma_last_attack", value = "Within the last year", HC = 0, RRMS = 0, `Prog MS` = 0) %>%
  add_row(var = "asthma_last_attack", value = "Within the last week", HC = 0, RRMS = 0, `Prog MS` = 0) %>%
  mutate(value = factor(value, 
                        levels = c("Yes", "Unsure or Never", "Over a year ago", "Within the last year", 
                                   "Within the last month","Within the last week"))) %>%
  arrange(var, value) %>%
  mutate(HC = HC/40, RRMS = RRMS/53, `Prog MS` = `Prog MS`/16)

# Export
write_tsv(x=table4.b, file = "../../03b-EbioMedicine-Revise/Adjusted_Figures_Tables/Table4b_v1.tsv")


