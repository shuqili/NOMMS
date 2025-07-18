library(tidyverse)
library(glue)
library(ggtext)
library(reshape2)
getwd()
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin/")

list.files(path = "03-output/00-anatomical_location/l6-ref-oropharyngeal/")
list.files(path = "03-output/00-anatomical_location/l6-ref_gingival/")


p <- 0.05
q <- 0.05
x_var <- "Anatomical_Location"

l6.ref.oro <-
  read_tsv("03-output/00-anatomical_location/l6-ref-oropharyngeal/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(reference = "T")

l6.ref.gin <-
  read_tsv("03-output/00-anatomical_location/l6-ref_gingival/cplmfit-no-transformation-relab-asv-to-site/all_results.tsv") %>%
  filter(value != "T") %>%
  mutate(reference ="G")

l6.t <-
  bind_rows(l6.ref.gin, l6.ref.oro) %>%
  filter(metadata=="Anatomical_Location") %>%
  select(feature, value, reference, pval, qval, coef, everything()) %>% 
  filter(pval < p) %>% #View()
  # filter(qval < q) %>% 
  unite("value", value:reference, na.rm = TRUE, remove = T, sep = " vs ") %>% #works for continuous and discrete
  mutate(value = factor(value, levels=c("N vs G", 
                                        "N vs T", 
                                        "G vs T"))) %>%
  mutate(feature=str_replace(feature, "d__Bacteria.*f__","f__"))

# bold 4 joint taxa in nasal (NvsG and NvsO), all genus level annotation
bold.set.1 <-
  l6.t %>%
  filter(value != "G vs T") %>%
  filter(coef > 0) %>%
  group_by(feature) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n>1) %>%
  distinct(feature) 

# bold 7 joint taxa in gingiva (NvsG and GvsO), all genus level annotation
bold.set.2 <-
  l6.t %>%
  filter(value != "N vs T") %>%
  mutate(coef=if_else(value == "N vs G", -coef, coef)) %>%
  filter(coef > 0) %>%
  group_by(feature) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1) %>% 
  distinct(feature) 

# bold12 joint taxa in in oropharyngeal (NvsO and GvsO), all genus level annotation.
bold.set.3 <-
  l6.t %>%
  filter(value != "N vs G") %>%
  filter(coef < 0) %>%
  group_by(feature) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n>1) %>% 
  distinct(feature) 

bold.set<-
  bind_rows(bold.set.1,bold.set.2, bold.set.3) %>%
  distinct() %>% pull

l6.t.final <-
  l6.t %>%
  mutate(feature = str_replace(feature, "f__Neisseriaceae.g__uncultured", "Uncultured ***Neisseriaceae***")) %>%
  mutate(feature = str_replace(feature, "f__P5D1.392.g__P5D1.392", "***Lactobacillales P5D1-392***")) %>%
  mutate(feature = str_replace(feature, "f__.*g__(.*)","*\\1*")) %>%
  mutate(feature = if_else(feature %in% bold.set, str_replace(feature, "(.*)","**\\1**"), feature)) %>%
  mutate(feature = str_replace(feature, "\\*.Eubacterium._nodatum_group\\*", "*Eubacterium nodatum group*")) %>%
  mutate(feature = str_replace(feature, "\\*.Eubacterium._brachy_group\\*", "*Eubacterium brachy group*"))          
  


# pairsise heatmaps #1.nasal vs gingival
value <- "N vs G"
l6.x <-
  l6.t.final %>%  
  filter(value =="N vs G") %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") 
  # 
l6.x %>%
  mutate(feature=fct_reorder(feature,coef,.desc = F)) %>%
  ggplot(aes(x=value,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),
                   position = "right") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(name = NULL,
                       colours = c("#7f3b08","#b35806","#e08214",
                                   # "#fdb863",
                                   "#f7f7f7",
                                   "#a6dba0","#1b7837","#00441b"),
                       na.value = "transparent",
                       breaks=c(-9,-4,-1, 0, 1, 4,9),
                       labels=c(-9,-4,"", 0, "", 4,9),
                       limits=c(-9,9)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_markdown(color = "black", face = "italic", size = rel(2)),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = 0,
        legend.text = element_text(face = "bold", size = rel(1)),
        legend.key.width = unit(0.6,"cm"))

getwd()
# spell out nasal, gingival, and oropharyngeal in legends, only bold consistent elements in taxa.
ggsave("../../01-NOMMS_Final_qiime_v202309/12-WriteUps-and-Figures/Figures/nNasalVsGingival-l6-bold.check.discrete.color.bar.pdf",
       dpi = 300, width = 4.5, height = 8)

# pairsise heatmaps #1.nasal vs oropharyngeal
value <- "N vs T"
l6.y<-
  l6.t.final %>%  
  filter(value =="N vs T") %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") # %>%
  
l6.y %>%
  mutate(feature = fct_reorder(feature, coef, .desc = F)) %>%
  ggplot(aes(x=value,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),
                   position = "right") +
  scale_x_discrete(expand = c(0,0)) +
  
  scale_fill_gradientn(name = NULL,
                       colours = c("#40004b","#762a83","#9970ab","#f7f7f7", "#5aae61","#1b7837","#00441b"),
                       na.value = "transparent",
                       breaks=c(-9,-4,0.03, 0, 0.03, 4, 9),
                       labels=c(-9,-4,  "", 0,   "", 4, 9),
                       limits=c(-9, 9)) + 
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_markdown(color = "black", face = "italic", size = rel(2)),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = 0,
        legend.text = element_text(face = "bold", size = rel(1)),
        legend.key.width = unit(0.6,"cm"))


ggsave("../../01-NOMMS_Final_qiime_v202309/12-WriteUps-and-Figures/Figures/nNasalVsOropharyngeal-l6-bold.check.discrete.color.bar.pdf",
       dpi = 300, width = 4.5, height = 8)


# pairsise heatmaps #3. gingival vs oropharyngeal
value <- "G vs T"
l6.z <-
  l6.t.final %>%  
  filter(value =="G vs T") %>%
  mutate(plabel= if_else(pval<0.05,"*","")) %>%
  mutate(qlabel=if_else(qval<0.05,"+","")) %>%
  unite("label", plabel:qlabel, sep = " ") 
  

l6.z %>% 
  mutate(feature=fct_reorder(feature, coef, .desc = F)) %>%
  ggplot(aes(x=value,y=feature, fill=coef)) +
  geom_tile(color="black") +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),
                   position = "right") +
  scale_x_discrete(expand = c(0,0)) +

  scale_fill_gradientn(name = NULL,
                       colours = c("#2d004b","#40004b","#542788","#762a83", "#f7f7f7","#e08214","#b35806","#7f3b08","#662506"),
                       na.value = "transparent",
                       breaks=c(-9, -4, -1, -0.1, 0, 0.1, 1, 4, 9),
                       labels=c(-9, -4, "",   "", 0, "", "", 4, 9),
                       limits=c(-9, 9))+  
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_markdown(color = "black", face = "italic", size = rel(2)),
        plot.title = element_markdown(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.justification = 0,
        legend.text = element_text(face = "bold", size = rel(1)),
        legend.key.width = unit(0.6,"cm"))



ggsave("../../01-NOMMS_Final_qiime_v202309/12-WriteUps-and-Figures/Figures/nGingivalVsOropharyngeal-l6-bold.check.discrete.color.bar.pdf",
       dpi = 300, width = 4.5, height = 8)

