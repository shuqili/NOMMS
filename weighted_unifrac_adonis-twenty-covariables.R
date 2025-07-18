# (qiime2-amplicon-2023.9) ➜  3-sites-d1000 qiime diversity adonis --i-distance-matrix 01-distance-matrix/weighted_unifrac_distance_matrix-metadata_v5_filtered.qza --m-metadata-file ../../../metadata-CL04-NOMMS-v5-qiime2.txt --p-formula Anatomical_Location+DiseaseStatusOracle_abbrev+DMT_v2+disease_duration_years+age+sex+bmi+race+ethnicity+smoking+dental_procedure_simplified+oral_disease_simplified+oral_hygiene_cumsum_score+vit_d_yn+autoimm_disease_yn+asthma_last_attack+rui_last_attack+allergy_symptoms_today+chronic_sinusitis_yn --p-n-jobs 4 --o-visualization 04-adonis/weighted_unifrac_adonis-twenty-covariables.qzv
# Saved Visualization to: 04-adonis/weighted_unifrac_adonis-twenty-covariables.qzv
# Downloaded weighted_unifrac_adonis-twenty-covariables.tsv file from view.qiime2.org online visualization site.
# Adjust headers as follows: "covariate	Df	SumsOfSqs	MeanSqs	F.Model	R2	Pr(>F)"


library(tidyverse)
library(ggplot2)
library(ggtext)
library(ggbreak)

getwd()
setwd("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/06-beta/oral-silva138_1/3-sites-d1000/04-adonis/weighted_unifrac_adonis-twenty-covariables/")
list.files()

# threshold

p = 0.05

read_tsv("weighted_unifrac_adonis-twenty-covariables.tsv") %>%
  mutate(label = if_else(`Pr(>F)` < p, "*", "")) %>%
  select(covariate, R2, label) %>%
  filter(covariate != "Residuals") %>%
  filter(covariate != "Total") %>% #arrange(R2)
  # mutate(log2R2 = log2(R2),
  #        log10R2 = log10(R2) )
  mutate(covariate = factor(covariate, # change the label visualization as designed
                            levels = c("Anatomical_Location","DiseaseStatusOracle_abbrev","DMT_v2",
                                       "disease_duration_years","age","sex", "bmi","race","ethnicity",
                                       "smoking","dental_procedure_simplified","oral_disease_simplified",
                                       "oral_hygiene_cumsum_score","vit_d_yn","autoimm_disease_yn",
                                       "asthma_last_attack","rui_last_attack","allergy_symptoms_today",
                                       "chronic_sinusitis_yn"),
                            labels = c("Anatomical<br>Location","MS Diagnosis","MS Treatment",
                                       "MS Duration","Age","Sex", "BMI","Race","Ethnicity",
                                       "Smoking<br>Preference","Dental","Oral Disease",
                                       "Oral Hygiene<br>Cumulative Score","Vitamin D<br>Supplementation","Autoimmune<br>Disease not MS",
                                       "Asthma Attack","Upper Respiratory<br>Infection","Allergy Symptoms",
                                       "Chronic Sinusitis"))) %>%
  mutate(covariate = fct_reorder(covariate, R2, .desc = F)) %>%
  ggplot(aes(x=R2, y=covariate)) +
  geom_col(fill="lightblue") +
  geom_line()+
  geom_text(aes(label=label),hjust = 0.6, vjust = 0.75, size=12) +
  geom_vline(xintercept = c(0.02, 0.5), linetype = 'dotted')+
  labs(title = NULL, 
       x = "R<sup>2</sup> Contribution to<br>Microbiome Variation",
       y = NULL) +
  scale_x_continuous(expand = c(0,0),
                     breaks = c(0, 0.01, 0.02, 0.5),
                     limits = c(0,0.5)) +
  scale_x_cut(breaks = c(0.02), which=c(1, 2), scales=c(0.6, 0.5), space = 0.25) +
  scale_y_discrete() +
  theme_classic() +
  theme(axis.title = element_markdown(face = "bold", colour = "black", hjust = 0.65),
        axis.text.x.top = element_blank(),
        axis.line.x.top = element_blank(),
        axis.ticks.x.top =element_blank(),
        # axis.ticks.x.bottom = element_blan,
        axis.text.x.bottom = element_markdown(face = "bold", colour = "black",
                                              angle = 45, hjust = 1),
        axis.text.y.left = element_markdown(face = "bold", colour = "black"),
        panel.background = element_rect(fill='transparent'))
# getwd()

ggsave("weighted_unifrac_adonis-twenty-covariables-v3.pdf",
       width = 4.6, height = 8, dpi = 300, bg='transparent')





