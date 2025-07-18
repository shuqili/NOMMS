library(tidyverse)
library(glue)
library(ggtext)
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")



# set up variables
p = 0.05
q = 0.05
x_var <- "DMT_v2"
title <- glue("rel_ab ~ {x_var} + Covariates (p < {p}, q < {q})")
title

# read in data to plot

nasal.l6.ref.Untreated <-
  read_tsv("03-output/01-nostril/l6taxa/DMT_N_over_5_ref_Untreated/cplmfit-no-transformation-relab-asv-to-DMT/all_results.tsv") %>%
  mutate(anatomical_site = "nasal") %>%
  filter(metadata == "DMT_v2") %>%
  mutate(reference = "Untreated") 

gin.l6.ref.Untreated <-
  read_tsv("03-output/02-gums/l6taxa/DMT_N_over_5_ref_Untreated/cplmfit-no-transformation-relab-asv-to-DMT/all_results.tsv") %>%
  mutate(anatomical_site = "gingival") %>%
  filter(metadata == "DMT_v2") %>%
  mutate(reference = "Untreated") 

oro.l6.ref.Untreated <-
  read_tsv("03-output/03-throat/l6taxa/DMT_N_over_5_ref_Untreated/cplmfit-no-transformation-relab-asv-to-DMT/all_results.tsv") %>%
  mutate(anatomical_site = "oropharyngeal") %>%
  filter(metadata == "DMT_v2") %>%
  mutate(reference = "Untreated") 

l6.taxa <-
  bind_rows(nasal.l6.ref.Untreated, gin.l6.ref.Untreated, oro.l6.ref.Untreated) %>% 
  select(feature, anatomical_site, metadata, value, reference, everything()) %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_site, levels = c("nasal", "gingival", "oropharyngeal"))) %>% 
  mutate(value = factor(value, levels = c("HC", "anti-CD20","DMF", "Fingolimod",
                                          "glatiramer acetate", "Interferon", "Natalizumab"))) %>%
  mutate(feature = str_remove(feature, "d__.*g__")) %>%
  mutate(feature = str_replace(feature, "P5D1.392", "Lactobacillales P5D1-392")) %>%
  mutate(feature = str_replace(feature, "F0332", "Actinomycetaceae F0332")) %>%
  mutate(feature = str_replace(feature, "F0058", "Paludibacteraceae F0058")) %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***"))

# "nasal" begins
l6.taxa %>%  
  filter(anatomical_site == "nasal") %>% 
  group_by(feature) %>%
  mutate(mean.coef = mean(coef)) %>%
  ungroup() %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  complete(feature, value) %>%
  mutate(feature = fct_reorder(feature, mean.coef, .desc = F)) %>%
  ggplot(aes(x=value,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0), position = "right") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#0096FF","white","#EA4025"),
                       na.value = "transparent",
                       breaks=c(-12,0,12),
                       limits=c(-12,12)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 30, vjust = 1, hjust=1),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold"))
ggsave("../../03b-EbioMedicine-Revise/github repository/03-differential analyses/06-visualization/nasal-DMT_ref_Untreated_l6.pdf",
       width = 4.5, height = 5.5, dpi = 300)

#"nasal" END

#"gingival" begins

l6.taxa %>%  
  filter(anatomical_site == "gingival") %>% 
  group_by(feature) %>%
  mutate(mean.coef = mean(coef)) %>%
  ungroup() %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>% # view
  complete(feature, value) %>%
  mutate(feature = fct_reorder(feature, mean.coef, .desc = F)) %>%
  ggplot(aes(x=value,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0), position = "right") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#0096FF","white","#EA4025"),
                       na.value = "transparent",
                       breaks=c(-12,0,12),
                       limits=c(-12,12)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 30, vjust = 1, hjust=1),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold"))
ggsave("../../03b-EbioMedicine-Revise/github repository/03-differential analyses/06-visualization/gingival-DMT_ref_Untreated_l6.pdf",
       width = 4.5, height = 5.5, dpi = 300)

#"gingival" END

#"oropharyngeal" begins

l6.taxa %>%  
  filter(anatomical_site == "oropharyngeal") %>% 
  group_by(feature) %>%
  mutate(mean.coef = mean(coef)) %>%
  ungroup() %>% 
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>%
  complete(feature, value) %>%
  mutate(feature = fct_reorder(feature, mean.coef, .desc = F)) %>%
  ggplot(aes(x=value,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0), position = "right") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#0096FF","white","#EA4025"),
                       na.value = "transparent",
                       breaks=c(-12,0,12),
                       limits=c(-12,12)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x  = element_markdown(face = "bold", color = "black", 
                                        angle = 30, vjust = 1, hjust=1),
        axis.text.y = element_markdown(color = "black"),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.text = element_text(face = "bold"))
ggsave("../../03b-EbioMedicine-Revise/github repository/03-differential analyses/06-visualization/oropharyngeal-DMT_ref_Untreated_l6.pdf",
       width = 4.5, height = 5.5, dpi = 300)

#"oropharyngeal" END




