library(tidyverse)
library(glue)
library(ggtext)
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")

# set up variables
p = 0.05
q = 0.05
x_var <- "edss change"
title <- glue("rel_ab ~ {x_var} + Covariates (p < {p}, q < {q})")
title

# Read in data to plot
nasal <- read_tsv("03-output/01-nostril/l6taxa/edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv") %>%
  mutate(anatomical_site = "nasal")
gin <- read_tsv("03-output/02-gums/l6taxa/edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")  %>%
  mutate(anatomical_site = "gingival")
oro <- read_tsv("03-output/03-throat/l6taxa/edss_change/cplmfit-no-transformation-relab-asv-to-edss_change/all_results.tsv")   %>%
  mutate(anatomical_site = "oropharyngeal")

l6.taxa <-
  bind_rows(nasal, gin, oro) %>%
  filter(metadata == "edss_change") %>%
  select(feature, anatomical_site, metadata, value, everything()) %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_site, levels = c("nasal", "gingival", "oropharyngeal"))) %>% 
  mutate(feature = str_remove(feature, "d__.*g__"))  %>%
  mutate(feature = str_replace(feature, ".Eubacterium._brachy_group", "Eubacterium brachy group")) %>%
  mutate(feature = str_replace(feature, "F0058", "Paludibacteraceae F0058")) %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***"))

#
l6.taxa %>%  
  unite("feature.site", feature:anatomical_site, remove = F, sep = ":") %>%
  group_by(anatomical_site) %>%
  arrange(-coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:10)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>% 
  mutate(feature_pos = if_else(feature.site == "***Paludibacteraceae F0058***:oropharyngeal",13, feature_pos)) %>%
  complete(anatomical_site, feature) %>%
  mutate(feature = fct_reorder(feature, feature_pos, .desc = F)) %>%
  ggplot(aes(x=anatomical_site,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   labels = c("N", "G", "O")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-3,3),
                       breaks = c(-3, 0, 3),
                       labels = c("-3","0","3")) +
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
ggsave("../../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/04-EDSS-change/EDSS-change-l6.pdf",
       width = 2.8, height = 4, dpi = 300)

