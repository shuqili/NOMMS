library(tidyverse)
library(glue)
library(ggtext)
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")

# Import Map for abbreviated ASV
asv.map <- read_tsv("04-blastn-significant-ASVs/ASV.abbreviation.list.v3.tsv")

# set up variables
p = 0.05
q = 0.05
x_var <- "DiseaseStat"
title <- glue("rel_ab ~ {x_var} + Covariates (p < {p}, q < {q})")
title

# read in data to plot

nasal.l6.ref.HC <- read_tsv("03-output/01-nostril/l6taxa/diagnoses_ref_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>% 
  mutate(anatomical_site = "nasal") %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(reference = "HC") 
nasal.l6.ref.RR <- read_tsv("03-output/01-nostril/l6taxa/diagnoses_ref_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_site = "nasal") %>%
  filter(value == "Prog MS") %>%
  mutate(reference = "RRMS") 
  
gingival.l6.ref.HC <- read_tsv("03-output/02-gums/l6taxa/diagnoses_ref_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>% 
  mutate(anatomical_site = "gingival") %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(reference = "HC") 
gingival.l6.ref.RR <- read_tsv("03-output/02-gums/l6taxa/diagnoses_ref_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_site = "gingival") %>%
  filter(value == "Prog MS") %>%
  mutate(reference = "RRMS") 

oro.l6.ref.HC <- read_tsv("03-output/03-throat/l6taxa/diagnoses_ref_HC/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>% 
  mutate(anatomical_site = "oropharyngeal") %>%
  filter(metadata == "DiseaseStatusOracle_abbrev") %>%
  mutate(reference = "HC") 
oro.l6.ref.RR <- read_tsv("03-output/03-throat/l6taxa/diagnoses_ref_RRMS/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_site = "oropharyngeal") %>%
  filter(value == "Prog MS") %>%
  mutate(reference = "RRMS") 

l6.taxa <-
  bind_rows(nasal.l6.ref.HC,nasal.l6.ref.RR, gingival.l6.ref.HC, gingival.l6.ref.RR, 
          oro.l6.ref.HC, oro.l6.ref.RR) %>%
  select(feature, anatomical_site, metadata, value, reference, everything()) %>%
  unite("value", value:reference, sep = " vs ", remove = T) %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_site, levels = c("nasal", "gingival", "oropharyngeal"))) %>%
  mutate(feature = str_remove(feature, "d__.*g__")) %>%
  mutate(feature = str_replace(feature, "RF39", "Bacilli RF39")) %>%
  mutate(feature = str_replace(feature, "F0058", "Paludibacteraceae F0058")) %>%
  mutate(feature = str_replace(feature, "uncultured", "Lachnospiraceae")) %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***")) 
  
  
# "RRMS vs HC" Starts  

l6.taxa %>%  
  filter(value == "RRMS vs HC") %>%
  # unite("feature.anatomical.site", feature:anatomical_site, remove = F, sep = ":" ) %>%
  group_by(anatomical_site) %>%
  arrange(coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:11)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  complete(anatomical_site, feature) %>%
  mutate(feature = fct_reorder(feature, feature_pos, .desc = F)) %>% 
  # view()
  ggplot(aes(x=anatomical_site,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  # labs(title = glue("{title}")) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0), position = "top",
                   breaks = c("nasal","gingival","oropharyngeal"),
                   labels = c("N","G","O")) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#0096FF","#0096FF","white","#F8D448","#F8D448" ),
                       na.value = "transparent",
                       breaks=c(-17, -4,0,4, 17),
                       labels=c("HC","","","","RRMS"),
                       limits=c(-17,17)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x.top  = element_markdown(face = "bold", color = "black", angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold"))

ggsave("../../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/02-diagnoses/RR_vs_HC_diagnoses_related_genus.pdf", 
       width = 2.1, height = 4.2, limitsize = FALSE)

#"RRMS vs HC" END



# "Prog MS vs HC" Start
l6.taxa %>%  
  filter(value == "Prog MS vs HC") %>%
  unite("feature.anatomical.site", feature:anatomical_site, remove = F, sep = ":" ) %>%
  group_by(anatomical_site) %>%
  arrange(coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:11)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  complete(anatomical_site, feature) %>%
  mutate(feature = fct_reorder(feature, feature_pos, .desc = F)) %>%
  ggplot(aes(x=anatomical_site,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0), position = "top",
                   breaks = c("nasal","gingival","oropharyngeal"),
                   labels = c("N","G","O")) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#0096FF","#0096FF","white","#EA4025","#EA4025"),
                       na.value = "transparent",
                       breaks=c(-5, -4,0,4, 5),
                       labels=c("HC","","","","Prog MS"),
                       limits=c(-5,5)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x.top  = element_markdown(face = "bold", color = "black", angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold"))

ggsave("../../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/02-diagnoses/PMS_vs_HC_diagnoses_related_genus.pdf", 
       width = 2.7, height = 4.2, limitsize = FALSE)

# "Prog MS vs HC" Ends



# "Prog MS vs RR" Starts
l6.taxa %>%  
  filter(value == "Prog MS vs RRMS") %>%
  # unite("feature.anatomical.site", feature:anatomical_site, remove = F, sep = ":" ) %>%
  group_by(anatomical_site) %>%
  arrange(coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:6)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  complete(anatomical_site, feature) %>%
  mutate(feature = fct_reorder(feature, feature_pos, .desc = F)) %>%
  ggplot(aes(x=anatomical_site,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0), position = "top",
                   breaks = c("nasal","gingival","oropharyngeal"),
                   labels = c("N","G","O")) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#F8D448","#F8D448", "white","#EA4025","#EA4025"),
                       na.value = "transparent",
                       breaks=c(-17, -4,0,4, 17),
                       labels=c("RRMS","","","","Prog MS"),
                       limits=c(-17,17)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x.top  = element_markdown(face = "bold", color = "black", angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0))

ggsave("../../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/02-diagnoses/PMS_vs_RR_diagnoses_related_genus.pdf", 
       width = 2.5, height = 2.8, limitsize = FALSE)

# "Prog MS vs RR" Ends.

    
