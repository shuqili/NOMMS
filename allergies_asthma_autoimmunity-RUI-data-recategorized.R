library(tidyverse)
library(lubridate)
library(fs)
library(readxl)
library(ggplot2)
library(ggtext)

getwd()
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/10-phenotypes")
metadata <- read_tsv("../metadata-CL04-NOMMS-v4.txt")
metadata_clinical <- read_xlsx('01-descriptive-data/NommsSamples030223-COMPILATION-FM_LMC_SL_reviewed_v4.xlsx')

metadata_MS <- metadata %>%
  select(SubjectID, DiseaseStatusOracle_abbrev, excluded) %>%
  filter(excluded == "no") %>%
  distinct()%>%
  # select(-excluded) %>%
  mutate(DiseaseStatusOracle_abbrev = factor(DiseaseStatusOracle_abbrev,
                                                               levels = c("HC","RRMS","Prog MS"))) %>%
  drop_na()


metadata_sample_date <-
  metadata_clinical %>%
  select(NOMMS_ID, SAMPLE_DATE) %>%
  distinct() %>%
  mutate(sample_date = as.Date(SAMPLE_DATE)) %>%
  rename(subject_ID = NOMMS_ID) %>%
  select(subject_ID,sample_date)
  



pretty_names = c("HC" = "Healthy<br>Subjects<br>(N=40)",
                 "RRMS"= "Relapsing<br>Remitting<br>MS Subjects<br>(N=53)",
                 "Prog MS" = "Progressive<br>MS Subjects<br>(N=10)")


var <- c("allergy_symptoms_today",
         "envir_allergy_yn",
         "envir_allergy_type___1",	"envir_allergy_type___2",	"envir_allergy_type___3",
         "envir_allergy_type___4",	"envir_allergy_type___5",	"envir_allergy_other",
         "food_allergy_yn",
         "food_allergy_desc___1",	"food_allergy_desc___2",	"food_allergy_desc___3",
         "food_allergy_desc___4",	"food_allergy_desc___5",	"food_allergy_desc___6",
         "food_allergy_desc___7",	"food_allergy_desc___8",	"food_allergy_desc___9",
         "food_allergy_other",	
         "last_resp_infection",	 "resp_infection_type",	
         "chronic_sinusitis_yn",
         "sinusitis_medication_yn",	"sinusitis_med",
         "asthma_yn",	 "asthma_last_attack",	
         "autoimm_disease_yn",
         "autoimm_desc___1",	"autoimm_desc___2",	"autoimm_desc___3",
         "autoimm_desc___4",	"autoimm_desc___5",	"autoimm_other",
         "antibiotics_lasttime",	
         "probiotics_yn",	"probiotic_freq",	"probiotic_type")

var1 <- c("allergy_symptoms_today",
          "envir_allergy_yn",
          "food_allergy_yn",
          "last_resp_infection",
          "chronic_sinusitis_yn",
          "asthma_yn",	 "asthma_last_attack",	
          "autoimm_disease_yn",
          "autoimm_desc___1",	"autoimm_desc___2",	"autoimm_desc___3",
          "autoimm_desc___4",	"autoimm_desc___5",	"autoimm_other")

var2 <- c("allergy_symptoms_today",
          "envir_allergy_yn",
          "food_allergy_yn",
          "chronic_sinusitis_yn",
          "asthma_yn",
          "autoimm_disease_yn")


var3 <- c("autoimm_disease_yn",
          "autoimm_desc___1",	"autoimm_desc___2",	"autoimm_desc___3",
          "autoimm_desc___4",	"autoimm_desc___5",	"autoimm_other")

var4 <- c("last_resp_infection","resp_infection_type",
          "asthma_last_attack")

# read_xlsx("01-descriptive-data/RedCap-NOMMS-20230510_1840_SL_compilation-bind_rows.xlsx") %>% #colnames()
#   select("subject_code", any_of(var2), any_of(var4)) %>%View
#   
#   left_join(., metadata_MS, by = c("subject_code" = "SubjectID")) %>%
#   select(subject_code, DiseaseStatusOracle_abbrev, everything())


allergies_sinusitis_asthma_autoimmune_YN <-
  read_xlsx("01-descriptive-data/RedCap-NOMMS-20230510_1840_SL_compilation-bind_rows.xlsx") %>% #colnames()
  select("subject_code", any_of(var2)) %>%
  filter(subject_code != "J100104")  %>% #subject didnt file questionnaire
  mutate(autoimm_disease_yn = if_else(subject_code == "J200089", 0, autoimm_disease_yn)) %>% 
  pivot_longer(-subject_code, names_to = "names", values_to = "value") %>%
  mutate(value = if_else(names != "chronic_sinusitis_yn",
                         case_when(value == 1 ~ "Yes",
                                   value == 2 ~ "No",
                                   value == 3 ~ "Unsure"),
                         case_when(value == 1 ~"Yes",
                                   value == 0 ~ "No"))) %>%
  mutate(value = if_else(value =="Yes", 1, 0)) %>%
  mutate(value = replace_na(value, 0)) %>% # fill missing values (NA) using zero indicates a "unsure/unknown" state.
  pivot_wider(id_cols = "subject_code", names_from = "names", values_from = "value", values_fill = 0) %>%
  left_join(., metadata_MS, by = c("subject_code" = "SubjectID")) %>% 
  drop_na()


allergies_sinusitis_asthma_autoimmune_ratio <- 
  allergies_sinusitis_asthma_autoimmune_YN %>%
  group_by(DiseaseStatusOracle_abbrev) %>%
  summarise(N_total = n(),
            N_allergies_today = sum(allergy_symptoms_today),
            N_envir_y = sum(envir_allergy_yn),
            N_food_y = sum(food_allergy_yn),
            N_chronic_sinusitis_y = sum(chronic_sinusitis_yn),
            N_asthma_y = sum(asthma_yn),  
            N_autoimm_disease_y = sum(autoimm_disease_yn),
            .groups = "drop") %>%
  group_by(DiseaseStatusOracle_abbrev) %>%
  summarise(ratio_allergies_today = 100 * N_allergies_today/N_total,
           ratio_envir_y = 100 * N_envir_y/N_total,
           ratio_food_y = 100 * N_food_y/N_total,
           ratio_chronic_sinusitis_y = 100 * N_chronic_sinusitis_y/N_total,
           ratio_asthma_y = 100 * N_asthma_y/N_total,
           ratio_autoimm_disease_y = 100 * N_autoimm_disease_y/N_total,
           .groups = "drop") 

allergies_sinusitis_asthma_autoimmune_ratio %>%
  pivot_longer(-DiseaseStatusOracle_abbrev, names_to = "name", values_to = "value") %>%
  group_by(name) %>%
  mutate(mean_value = mean(value)) %>%
  ungroup() %>%
    mutate(name= factor(name,
                       levels=c("ratio_allergies_today",
                                "ratio_envir_y",
                                "ratio_food_y",
                                "ratio_chronic_sinusitis_y",
                                "ratio_asthma_y",
                                "ratio_autoimm_disease_y"),
                       labels=c("Allergic Symptoms on Visit Day",
                                "Environmental Allergies",
                                "Food Allergies",
                                "Chronic Sinusitis",
                                "Asthma",
                                "Autoimmune Disease besides MS"))) %>%
    # mutate(name =fct_reorder(name, mean_value,.desc = T)) %>%
    ggplot(aes(x=DiseaseStatusOracle_abbrev, y = value, 
               fill=name)) +
    geom_bar(position="dodge", stat="identity") +
    labs(x=NULL, y = "Percentage of Subjects in Group (%)") +
    scale_x_discrete() +
    scale_y_continuous(n.breaks = 10, limits = c(0,100), expand = c(0,0)) +
    scale_fill_manual(name=NULL,
                      values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c',
                                 '#fb9a99','#e31a1c','#fdbf6f','#ff7f00',
                                 '#cab2d6','#6a3d9a','#ffff99','#b15928',
                                 'darkgrey')) +
    theme_classic() +
    theme(legend.text = element_markdown(face = "bold"),
          axis.title =  element_markdown(face = "bold", color = "black"),
          axis.text = element_markdown(face = "bold", color = "black"),
          strip.background = element_blank(),
          strip.text = element_markdown(face = "bold", vjust = 1),
          plot.background = element_blank())
  
  
ggsave("03-outputs/allergies_sinusitis_asthma_autoimmune_ratio_barchart.pdf",width = 9, height = 4)
write_tsv(x = allergies_sinusitis_asthma_autoimmune_YN,
          file = "03-outputs/allergies_sinusitis_asthma_autoimmune_YN_table.tsv")



autoimmune_types <-
read_xlsx("01-descriptive-data/RedCap-NOMMS-20230510_1840_SL_compilation-bind_rows.xlsx") %>% #colnames()
  select("subject_code", any_of(var3)) %>%
  filter(subject_code != "J100104") %>%
  rename(type_I_diabetes = autoimm_desc___1,
         lupus= autoimm_desc___2,
         rheumatoid_arthritis= autoimm_desc___3,
         thyroid_disease= autoimm_desc___4,
         other= autoimm_desc___5) %>%
  mutate(autoimm_disease_yn = if_else(subject_code == "J200089", 0, autoimm_disease_yn)) %>% 
  pivot_wider(id_cols = c(1:6), names_from = "autoimm_other", values_from = "other", values_fill = 0) %>%# colnames()
  select(-Shingles) %>% # not considered autoimmune related. Shingles are caused by the reactivation of dormant varicella zoster and are a complication of having had chickenpox.
  rename(psoriatic_arthritis = `Psoriatic arthritis`) %>%
  rename(retinal_vasculitis = `Panuveitis, retil vasculitis, maybe cogans?`) %>% 
  # Panuveitis (removed), also called Total uveitis, is an eye disease affecting the internal structures of the eye. 
  # cogans removed due to uncertainty.
  # retinal vasculitis (kept) for it is majorly explained by Autoimmune mechanisms.
  select(-`Multiple Sclerosis`) %>%
  select(-`Fibromyalgia, borderline positive thyroid antibiodies`) %>% # Fibromyalgia NOT considered as autoimmune.
  select(-Eczema) %>%
  rename(cns_vasculitis = `CNS vasculitis`) %>% # somewhat likely autoimmune
  rename(autoimmune_hepatitis = `Auto Immune Hepatitis`) %>%
  select(-`NA`, -autoimm_disease_yn) %>%
  replace(is.na(.), 0)

autoimmune_types



autoimmune_numbers<-
  autoimmune_types %>%
  left_join(., metadata_MS, by = c("subject_code" = "SubjectID")) %>%
  
  pivot_longer(-c(subject_code,DiseaseStatusOracle_abbrev),
               names_to = "name", values_to = "value") %>%
  group_by(DiseaseStatusOracle_abbrev, name) %>%
  summarise(N_total = n(), 
            value = sum(value),
            .groups = "drop") %>%
  drop_na() 

autoimmune.table <- autoimmune_numbers %>%
  select(-N_total) %>%
  pivot_wider(names_from = DiseaseStatusOracle_abbrev, values_from = value) %>%
  mutate(HC = HC/40, RRMS = RRMS/53, `Prog MS` = `Prog MS`/16)

autoimmune_numbers%>%
  mutate(DiseaseStatusOracle_abbrev = factor(DiseaseStatusOracle_abbrev,
                                             levels = c("HC", "RRMS", "Prog MS")),
         value = as.factor(value)) %>% 
  mutate(name= factor(name,
                   levels=c("type_I_diabetes",
                            "thyroid_disease", 
                            "rheumatoid_arthritis",
                            "retinal_vasculitis",
                            "psoriatic_arthritis", 
                            "lupus",
                            "cns_vasculitis", 
                            "autoimmune_hepatitis"),
                   labels=c("Type I Diabetes",
                            "Thyroid disease", 
                            "Rheumatoid Arthritis",
                            "Retinal Vasculitis",
                            "Psoriatic Arthritis",
                            "Lupus",
                            "CNS Vasculitis", 
                            "Autoimmune Hepatitis"))) %>%
  ggplot(aes(x=DiseaseStatusOracle_abbrev, y = name, fill=value)) +
  geom_tile(color = "black", linewidth = 0.4) +
  labs(x=NULL, y = NULL) +
  scale_x_discrete() +
  scale_y_discrete(position = "right") +
  scale_fill_manual(name=NULL,
                    
                    breaks = c(0,1), values = c("white","darkgreen")) +
  theme_classic() +
  theme(legend.text = element_markdown(face = "bold"),
        axis.title =  element_markdown(face = "bold", color = "black"),
        axis.text.y.right = element_markdown(face = "bold", color = "black"),
        axis.text.x = element_markdown(face = "bold", color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none")

ggsave("03-outputs/autoimmune_numbers_tile_plot.pdf",width = 3.2, height = 4)
write_tsv(x=autoimmune_types, file = "03-outputs/autoimmune_types_table.tsv")
write_tsv(x=autoimmune_numbers, file = "03-outputs/autoimmune_numbers_of_groups_table.tsv")

asthma_rui_score <-
read_xlsx("01-descriptive-data/RedCap-NOMMS-20230510_1840_SL_compilation-bind_rows.xlsx") %>% #colnames()
  select("subject_code", any_of(var4)) %>%
  filter(subject_code != "J100104") %>%
  left_join(., metadata_sample_date, by = c("subject_code" = "subject_ID")) %>% #colnames()
  mutate(last_resp_infection_date = as.Date(last_resp_infection)) %>%
  select(subject_code, sample_date, last_resp_infection_date, asthma_last_attack, resp_infection_type) %>%
  # drop_na() %>%
  mutate(time_diff_unit_month = as.numeric(sample_date-last_resp_infection_date)/ 30) %>%
  mutate(rui_last_attack = case_when(time_diff_unit_month <0.24 ~1,
                                               0.24<time_diff_unit_month & time_diff_unit_month <1 ~2,
                                               1<time_diff_unit_month & time_diff_unit_month  <12 ~3,
                                               12<time_diff_unit_month  ~4),
         rui_last_attack = replace_na(rui_last_attack,5)) %>%
  select(subject_code, rui_last_attack, asthma_last_attack) %>%
  mutate(asthma_last_attack = replace_na(asthma_last_attack, 5)) %>%
  left_join(., metadata_MS, by = c("subject_code"= "SubjectID")) %>%
  drop_na(DiseaseStatusOracle_abbrev) %>%
  pivot_longer(-c(DiseaseStatusOracle_abbrev, subject_code), 
               names_to = "name", values_to = "value") %>%
  mutate(meaning = case_when(value == 1 ~ "Within the last week",
                             value == 2 ~ "Within the last month",
                             value == 3 ~ "Within the last year",
                             value == 4 ~ "Over a year ago",
                             value == 5 ~ "Unsure or Never")) %>%
  mutate(value_new = case_when(meaning == "Unsure or Never" ~ 0,
                               meaning == "Over a year ago" ~ 1,
                               meaning == "Within the last year" ~ 2,
                               meaning == "Within the last month" ~ 3,
                               meaning == "Within the last week" ~ 4)) %>%
  select(subject_code, DiseaseStatusOracle_abbrev, name, value_new) %>%
  pivot_wider(id_cols = c("subject_code", "DiseaseStatusOracle_abbrev"), 
              names_from = "name", values_from = "value_new") 

allergy_sinusitis_asthma_rui_score <-
  inner_join(allergies_sinusitis_asthma_autoimmune_YN,asthma_rui_score,
           by = c("subject_code","DiseaseStatusOracle_abbrev")) %>%
  select(DiseaseStatusOracle_abbrev, subject_code,
         allergy_symptoms_today,
         chronic_sinusitis_yn,
         asthma_last_attack,
         rui_last_attack)


allergy_sinusitis_asthma_rui_score %>% 
  rowwise() %>%
  mutate(pertubation_subject_score = allergy_symptoms_today+ chronic_sinusitis_yn+ asthma_last_attack+ rui_last_attack) %>%
  ungroup() %>% #View
  mutate(pertubation_cutoff = 4) %>%
  mutate(pertubated = pertubation_subject_score > pertubation_cutoff) %>%
  filter(pertubated == TRUE)
  
# write_tsv(x = allergy_sinusitis_asthma_rui_score, file = "03-outputs/allergy_sinusitis_asthma_rui_score.tsv")  












# Export
write_tsv(x=autoimmune.table, file = "../../03b-EbioMedicine-Revise/Adjusted_Figures_Tables/Table4-autoimmune-table_v1.tsv")
