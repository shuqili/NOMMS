library(tidyverse)
library(glue)
library(ggtext)
library(ggplot2)
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")

# set up variables
p = 0.05
q = 0.05
x_var <- "edss change"
title <- glue("rel_ab ~ {x_var} + Covariates (p < {p}, q < {q})")
title

# Read in data to plot
female.nasal <- read_tsv("03-output/01-nostril/l6taxa/female-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  mutate(anatomical_site = "nasal")
female.gin <- read_tsv("03-output/02-gums/l6taxa/female-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")  %>%
  mutate(anatomical_site = "gingival")
female.oro <- read_tsv("03-output/03-throat/l6taxa/female-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")   %>%
  mutate(anatomical_site = "oropharyngeal")

male.nasal <- read_tsv("03-output/01-nostril/l6taxa/male-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  mutate(anatomical_site = "nasal")
male.gin <- read_tsv("03-output/02-gums/l6taxa/male-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")  %>%
  mutate(anatomical_site = "gingival")
male.oro <- read_tsv("03-output/03-throat/l6taxa/male-edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")   %>%
  mutate(anatomical_site = "oropharyngeal")

females.res <-
  bind_rows(female.nasal, female.gin, female.oro) %>%
  mutate(sex = "female") 

males.res <-
  bind_rows(male.nasal, male.gin, male.oro) %>%
  mutate(sex = "male")


l6.taxa <-
  bind_rows(females.res, males.res) %>%
  filter(metadata == "edss_change") %>%
  select(feature, anatomical_site, metadata, value, everything()) %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_site, levels = c("nasal", "gingival", "oropharyngeal"))) %>% 
  mutate(sex = factor(sex, levels = c("female", "male"))) %>% 
  mutate(feature = str_remove(feature, "d__.*g__"))  %>% 
  mutate(feature = str_replace(feature, ".Eubacterium._brachy_group", "Eubacterium brachy group")) %>%
  mutate(feature = str_replace(feature, "d__Bacteria.p__Firmicutes.c__Negativicutes.o__Veillonellales.Selenomonadales.f__Selenomonadaceae.__",
                               "Selenomonadaceae sp.")) %>%
  mutate(feature = str_replace(feature, "F0058", "Paludibacteraceae F0058")) %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***"))

#

l6.taxa %>%  
  unite("feature.site", feature:anatomical_site, remove = F, sep = ":") %>% 
  filter(sex == "female") %>%
  group_by(anatomical_site) %>%
  arrange(-coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:6)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>% 
  # group_by(sex) %>%
  complete(anatomical_site, feature.site) %>%
  # ungroup() %>%
  mutate(feature.site = fct_reorder(feature.site, feature_pos,
                                    .desc = F, .na_rm = TRUE)) %>%
  ggplot(aes(x=anatomical_site,y=feature.site, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  # facet_wrap(vars(sex), nrow = 2) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-6,6),
                       breaks = c(-6, 0, 6),
                       labels = c("-6","0","6")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0.5))
getwd()
ggsave("../../03b-EbioMedicine-Revise/00-github repository/04-sex-stratified-differential-analyses/01-female-sub_cohort/l6-female-edss-change.pdf",
       width = 3.4, height = 4, dpi = 300)


l6.taxa %>%  
  unite("feature.site", feature:anatomical_site, remove = F, sep = ":") %>% 
  filter(sex == "male") %>%
  group_by(anatomical_site) %>%
  arrange(-coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:6)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>% 
  # group_by(sex) %>%
  complete(anatomical_site, feature.site) %>%
  # ungroup() %>%
  mutate(feature.site = fct_reorder(feature.site, feature_pos,
                                    .desc = F, .na_rm = TRUE)) %>%
  ggplot(aes(x=anatomical_site,y=feature.site, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  # facet_wrap(vars(sex), nrow = 2) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-6,6),
                       breaks = c(-6, 0, 6),
                       labels = c("-6","0","6")) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold", hjust = 0.5))
getwd()
ggsave("../../03b-EbioMedicine-Revise/00-github repository/04-sex-stratified-differential-analyses/02-male-sbu_cohort/l6-male-edss-change.pdf",
       width = 3.4, height = 2.4, dpi = 300)
