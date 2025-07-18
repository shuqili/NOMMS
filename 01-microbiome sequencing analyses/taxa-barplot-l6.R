library(tidyverse)
library(ggtext)

getwd()
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/07-taxa/05-count-tables/")
list.files()

l6.count <-
  read_tsv("filtered-genus-count-table-mapped-w-human-oral-weighted-silva138_1.txt") %>%
  rename(sample_id = index) %>% #colnames()
  select(c(1:137), SubjectID, Anatomical_Location, DiseaseStatusOracle_abbrev) 
metadata <-
  l6.count %>%
  select(sample_id, SubjectID, Anatomical_Location, DiseaseStatusOracle_abbrev) %>%
  mutate(DiseaseStatusOracle_abbrev = factor(DiseaseStatusOracle_abbrev, 
                                             levels = c("HC", "RRMS", "Prog MS"))) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, 
                                      levels = c("N", "G", "T"),
                                      labels = c("N", "G", "O")))

l6.rel.abund <-
  l6.count %>%
  select(c(1:137)) %>%
  pivot_longer(-sample_id, names_to = "taxa", values_to = "count") %>%
  group_by(sample_id) %>%
  mutate(rel_ab = count/sum(count) * 100) %>%
  ungroup() %>% 
  inner_join(., metadata, by = "sample_id") %>%
  select(Anatomical_Location, DiseaseStatusOracle_abbrev, sample_id, taxa, rel_ab) %>%
  filter(DiseaseStatusOracle_abbrev != "NA") %>%
  group_by(Anatomical_Location, DiseaseStatusOracle_abbrev, taxa) %>%
  reframe(mean_rel_ab = mean(rel_ab)) %>%
  mutate(taxa = str_replace(taxa, "d__.*o__", "o__")) %>% 
  mutate(taxa = str_replace(taxa, "(.*)f__(.*);__$", "Unclassified \\2")) %>% 
  mutate(taxa = str_replace(taxa, "(.*)o__(.*);g__uncultured$", "Uncultured \\2")) %>% # 
  mutate(taxa = if_else(str_detect(taxa, "f__uncultured"), 
                        str_remove(taxa, ";f__uncultured"), 
                        str_remove(taxa, "Uncultured.*;f__"))) %>% 
  mutate(taxa = str_remove(taxa, "o__.*g__")) %>%
  mutate(taxa = str_replace(taxa, "(.*)", "***\\1***"))
  

genus_pool <-
  l6.rel.abund %>%
  group_by(taxa) %>% #summarise(max = max(mean_rel_abund)) %>% arrange(desc(max)) %>% pull
  summarise(pool = max(mean_rel_ab) < 3, .groups = "drop")

l6.rel.abund.pooled <-
  l6.rel.abund %>%
  inner_join(., genus_pool, by = "taxa") %>%
  mutate(taxa = if_else(pool, "**Other (< 3%)**", taxa)) %>%
  group_by(Anatomical_Location, DiseaseStatusOracle_abbrev, taxa) %>% 
  reframe(mean_rel_abund = sum(mean_rel_ab)) 


#####

l6.taxa.bar.plot <-
l6.rel.abund.pooled %>%
  mutate(taxa =fct_reorder(taxa, mean_rel_abund,.desc = T)) %>%
  mutate(taxa =fct_shift(taxa,n = 2)) %>%
  ggplot(aes(x=DiseaseStatusOracle_abbrev, y=mean_rel_abund, fill = taxa)) +
  geom_col() +
  facet_wrap(~Anatomical_Location, nrow = 3) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(name="Mean Relative Abundance (%)<br>of Genera", #"Human Oral<br>Weighted<br>SILVA 138.1",
                      values= c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#1f78f4','#ffff99','#b15928',
                                '#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69', "darkgray",'#ffed6f','#bf812d',
                                '#543005','#8c510a','#bc80bd','#fccde5','#d9d9d9','#ccebc5','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30',
                                '#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a'))+
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0),n.breaks = 10)+
    theme_classic() +
    guides(fill = guide_legend(ncol = 1)) +
    theme(legend.title = element_markdown(face = "bold"),
          legend.text = element_markdown(),
          legend.justification = "top",
          # legend.key.size = unit(10,"pt"),
          strip.text.x.top = element_text(face = "bold", size = 9),
          axis.text.y = element_markdown(face = "bold", size = 12, color = "black"),
          axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                     angle = 45, hjust = 1, vjust = 1),
          axis.title.y.left = element_text(face = "bold", size = 16))

getwd()

l6.taxa.bar.plot
ggsave("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/03b-EbioMedicine-Revise/02-Adjusted_Figures_Tables/l6_rel_abundance_stacked_barplot-v4.pdf",
       dpi = 300, width = 4, height = 10)
