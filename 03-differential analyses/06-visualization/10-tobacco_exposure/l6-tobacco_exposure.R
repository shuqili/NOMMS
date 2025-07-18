library(tidyverse)
library(ggtext)

getwd()
setwd("09-MaAsLin/03-output/")

nasal <-
  read_tsv("01-nostril/l6taxa/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_location = "N")

gin <-
  read_tsv("02-gums/l6taxa/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_location = "G")

oro <-
  read_tsv("03-throat/l6taxa/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_location = "O")


# variables
p = 0.05
q = 0.05

l6.taxa <-
bind_rows(nasal, gin, oro) %>%
  filter(metadata == "tobacco_exposure") %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_location, levels = c("N", "G", "O"))) %>% 
  mutate(feature = str_remove(feature, "d__.*g__"))  %>% 
  mutate(feature = str_replace(feature, "d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Pasteurellaceae.__", "Pasteurellaceae spp.")) %>%
  mutate(feature = str_replace(feature, ".Eubacterium._brachy_group", "Eubacterium brachy group")) %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***"))

#
l6.taxa %>%  
  select(feature, anatomical_site, everything()) %>%
  unite("feature.site", feature:anatomical_site, remove = F, sep = ":") %>%
  group_by(anatomical_site) %>%
  arrange(-coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:12)) %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") %>% 
  complete(anatomical_site, feature) %>%
  mutate(feature = fct_reorder(feature, feature_pos, .desc = F)) %>%
  ggplot(aes(x=anatomical_site,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0),position = "top",
                   breaks = c("N", "G", "O"),
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
ggsave("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/10-tobacco_exposure/tobacco_exposure-l6.pdf",
       width = 2.5, height = 4, dpi = 300)

