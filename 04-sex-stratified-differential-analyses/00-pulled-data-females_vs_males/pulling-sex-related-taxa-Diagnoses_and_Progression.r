library(tidyverse)
library(ggplot2)
library(ggtext)
library(stringr)
getwd()
setwd("Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")


# Define statistical thresholds.
p = 0.05
q = 0.05


# l6-genera; diagnosis
# Import files
nasal.diagnose.sex.ref.female <-
  read_tsv("03-output/01-nostril/l6taxa/diagnoses_ref_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "nasal")
gin.diagnose.sex.ref.female <-
  read_tsv("03-output/02-gums/l6taxa/diagnoses_ref_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "gingival")
oro.diagnose.sex.ref.female <-
  read_tsv("03-output/03-throat/l6taxa/diagnoses_ref_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "oropharyngeal")

diagnose.sex.ref.female <-
  bind_rows(nasal.diagnose.sex.ref.female, gin.diagnose.sex.ref.female,oro.diagnose.sex.ref.female) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, levels = c("nasal", "gingival", "oropharyngeal"))) %>%
  filter(pval < p) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  mutate(feature = str_remove(feature, "d__.*o__")) %>% 
  mutate(feature = str_replace(feature, "Bacteroidales.f__Paludibacteraceae.g__F0058", "***Paludibacteraceae*** **F0058**")) %>%
  mutate(feature = str_replace(feature, "Lactobacillales.f__P5D1.392.g__P5D1.392", "***Lactobacillales*** **P5D1-392**")) %>%
  mutate(feature = str_replace(feature, "(.*g__)(.*)", "***\\2***"))


diagnose.sex.ref.female %>%
  arrange(Anatomical_Location, coef) %>%
  mutate(feature.pos = c(1:15)) %>%
  complete(Anatomical_Location, feature) %>% 
  mutate(feature = fct_reorder(feature, feature.pos, .desc = F)) %>%
  ggplot(aes(x=Anatomical_Location, y= feature, fill = coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with sex",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-3,3),
                       breaks = c(-3, 0, 3),
                       labels = c("females","0","males")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0.2))
ggsave("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/02-diagnoses/l6taxa_more_abundant_in_males_independent_of_diagnosis.pdf",
       width = 2.5, height = 4.2, dpi = 300)


# l6-genera; progression
# Import files
nasal.progression.sex.ref.female <-
  read_tsv("03-output/01-nostril/l6taxa/edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "nasal")
gin.progression.sex.ref.female <-
  read_tsv("03-output/02-gums/l6taxa/edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "gingival")
oro.progression.sex.ref.female <-
  read_tsv("03-output/03-throat/l6taxa/edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "oropharyngeal")

progression.sex.ref.female <-
  bind_rows(nasal.progression.sex.ref.female, gin.progression.sex.ref.female,oro.progression.sex.ref.female) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, levels = c("nasal", "gingival", "oropharyngeal"))) %>%
  filter(pval < p) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  mutate(feature = str_remove(feature, "d__.*o__")) %>% 
  mutate(feature = str_replace(feature, "Actinomycetales.f__Actinomycetaceae.g__F0332", "***Actinomycetales*** **F0332**")) %>%
  mutate(feature = str_replace(feature, "(.*g__)(.*)", "***\\2***"))


progression.sex.ref.female %>%
  arrange(Anatomical_Location, coef) %>%
  mutate(feature.pos = c(1:10)) %>%
  complete(Anatomical_Location, feature) %>%
  mutate(feature = fct_reorder(feature, feature.pos, .desc = F)) %>%
  ggplot(aes(x=Anatomical_Location, y= feature, fill = coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with sex",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-3,3),
                       breaks = c(-3, 0, 3),
                       labels = c("females","0","males")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0.2))
ggsave("../../03b-EbioMedicine-Revise/l6taxa_more_abundant_in_males_independent_of_progression.pdf",
       width = 2.8, height = 3.4, dpi = 300)


# l8 - ASV



# Diagnoses
# Import files

asv.abbrev <-
read_tsv("../09-MaAsLin/04-blastn-significant-ASVs/ASV.abbreviation.list.v3.tsv") %>% 
  select(number, feature)

nasal.diagnose.sex.ref.female <-
  read_tsv("03-output/01-nostril/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "nasal")
gin.diagnose.sex.ref.female <-
  read_tsv("03-output/02-gums/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "gingival")
oro.diagnose.sex.ref.female <-
  read_tsv("03-output/03-throat/ASVlevel/DiseaseStat + [age + sex + bmi + smoking + dental_procedure_simplified + oral_disease_simplified + oral_hygiene_cumsum_score]/RRMS_vs_HC_and_PMS_vs_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "oropharyngeal")

diagnose.sex.ref.female <-
  bind_rows(nasal.diagnose.sex.ref.female, gin.diagnose.sex.ref.female,oro.diagnose.sex.ref.female) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, levels = c("nasal", "gingival", "oropharyngeal"))) %>%
  filter(pval < p) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>% 
  mutate(feature = str_replace(feature, "(d__.*)(g__.*)(s__.*)(ASV__.*)", "\\2\\3\\4")) %>% 
  mutate(feature = if_else(str_detect(feature, "s__un"), str_replace(feature, "(g__)(.*)(s__.*)(ASV__.*)", "\\2 \\4"), 
                           str_replace(feature, "(g__)(.*)(s__)(.*)_(.*)(ASV__.*)", "\\4 \\5 \\6" ))) %>% 
  mutate(feature = str_replace(feature, "F0058. ASV__20f887a0f58476d52027e0798805c1d2", 
                               "Paludibacteraceae F0058 ASV__20f887a0f58476d52027e0798805c1d2")) %>%
  mutate(feature = str_replace(feature, "P5D1.392. ASV__386a607b118090dd71bfc8b4b4081447", 
                               "Lactobacillales P5D1-392 ASV__386a607b118090dd71bfc8b4b4081447")) %>% 
  mutate(feature = str_remove(feature, "\\.")) %>%
    separate(feature, into = c("species", "ASV"), sep = " ASV__", remove = F) %>% #select(species, ASV) %>%
    left_join(., asv.abbrev, by = c("ASV" = "feature")) %>% #  %>% view
    select(species, number, feature, ASV, everything()) %>%
    unite("feature.abbrev", species:number, sep = " ") %>%
    select(-ASV) %>% #select(feature.abbrev) %>%
    mutate(feature.abbrev = str_replace(feature.abbrev, "(.*) (ASV_.*)", "***\\1*** **\\2**")) 
  
  # mutate(feature = str_replace(feature, "(.*) (ASV__.*)", "***\\1*** **\\2**")) 
  
  


diagnose.sex.ref.female %>%
  arrange(Anatomical_Location, coef) %>%
  mutate(feature.pos = c(1:45)) %>% #view()
  complete(Anatomical_Location, feature.abbrev) %>%
  mutate(feature.abbrev = fct_reorder(feature.abbrev, feature.pos, .desc = F)) %>%
  ggplot(aes(x=Anatomical_Location, y= feature.abbrev, fill = coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with sex",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-4,4),
                       breaks = c(-4, 0, 4),
                       labels = c("females","0","males")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", size = 12,
                                        angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black", size = 10),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0.2))
ggsave("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/02-diagnoses/l8_ASV_more_abundant_in_males_independent_of_diagnosis.pdf",
       width = 4, height = 12, dpi = 300)


# l8-ASV; progression
# Import files
nasal.progression.sex.ref.female <-
  read_tsv("03-output/01-nostril/ASVlevel/edss_change + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "nasal")
gin.progression.sex.ref.female <-
  read_tsv("03-output/02-gums/ASVlevel/edss_change + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "gingival")
oro.progression.sex.ref.female <-
  read_tsv("03-output/03-throat/ASVlevel/edss_change + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  filter(metadata == "sex") %>%
  mutate(Anatomical_Location = "oropharyngeal")

progression.sex.ref.female <-
  bind_rows(nasal.progression.sex.ref.female, gin.progression.sex.ref.female,oro.progression.sex.ref.female) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, levels = c("nasal", "gingival", "oropharyngeal"))) %>%
  filter(pval < p) %>% 
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  mutate(feature = str_replace(feature, "(d__.*)(g__.*)(s__.*)(ASV__.*)", "\\2\\3\\4")) %>% 
  mutate(feature = if_else(str_detect(feature, "s__un"), str_replace(feature, "(g__)(.*)(s__.*)(ASV__.*)", "\\2 \\4"), 
                           str_replace(feature, "(g__)(.*)(s__)(.*)_(.*)(ASV__.*)", "\\4 \\5 \\6" ))) %>% 
  mutate(feature = str_replace(feature, "F0332. ASV__cad533d4ec4ec33683ca5652ff1b32b2",
                               "Actinomycetaceae F0332. ASV__cad533d4ec4ec33683ca5652ff1b32b2")) %>% 
  mutate(feature = str_remove(feature, "\\.")) %>% 
  mutate(feature = str_replace(feature, "(.*) (ASV__.*)", "***\\1*** **\\2**")) 
  
  
progression.sex.ref.female %>%
  arrange(Anatomical_Location, coef) %>%
  mutate(feature.pos = c(1:31)) %>% 
  complete(Anatomical_Location, feature) %>%
  mutate(feature = fct_reorder(feature, feature.pos, .desc = F)) %>%
  ggplot(aes(x=Anatomical_Location, y= feature, fill = coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with sex",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-6,6),
                       breaks = c(-6, 0, 6),
                       labels = c("females","0","males")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0.2))

ggsave("../../03b-EbioMedicine-Revise/l8_ASV_more_abundant_in_males_independent_of_progression.pdf",
       width = 6, height = 8, dpi = 300)
