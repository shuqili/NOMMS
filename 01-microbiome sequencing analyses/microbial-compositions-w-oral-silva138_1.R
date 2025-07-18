library(readxl)
library(tidyverse)
library(RColorBrewer)
library(ggtext)
getwd()
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/07-taxa/05-count-tables")
list.files()
read_tsv("filtered-ASV-count-table-w-human-oral-weighted-silva138_1.tsv", skip = 1) %>% colnames()
# read_csv("filtered-count-table-mapped-w-human-oral-weighted-silva138_1.csv") %>% colnames()



metadata <-
  read_tsv("../../metadata-CL04-NOMMS-v4.txt") %>%
  select(c(`#Sample ID`,"SubjectID", "Anatomical_Location", "DiseaseStatusOracle_abbrev" , "sex", "excluded")) %>%
  rename(sample_ID = `#Sample ID`) %>%
  mutate(DiseaseStatusOracle_abbrev= factor(DiseaseStatusOracle_abbrev,  
                                            levels=c("HC", "RRMS", "Prog MS"),
                                            labels=c("HC", "RRMS", "Prog MS"))) %>%
  mutate(Anatomical_Location = factor(Anatomical_Location, 
                                      levels=c("N", "G", "T"),
                                      labels=c("N", "G", "O"))) %>%
  mutate(sex = factor(sex, 
                      levels = c("M","F"),
                      labels = c("Males","Females"))) 


asv_counts <-
  read_tsv("filtered-ASV-count-table-w-human-oral-weighted-silva138_1.tsv", skip = 1) %>%
  rename(feature_ID = `#OTU ID`) %>%
  pivot_longer(-feature_ID, names_to = "sample_ID", values_to = "count")

taxonomy <-
  read_tsv("../02-mapped_taxonomy/taxonomy-mapped-w-human-oral-weighted-silva138_1.tsv") %>%
  rename_all(tolower) %>%
  rename(feature_ID = `feature id`) %>%
  separate(taxon, sep = "; ", fill = "warn", extra = "warn",
           into = c("kingdom", "phylum", "class", "order", 
                           "family", "genus", "species")) %>%
  mutate(order = str_replace_na(order, "o__unclassified")) %>% #summarize(count = sum(is.na(order)))
  mutate(family = str_replace_na(family, "f__unclassified")) %>% #summarize(count = sum(is.na(family)))
  mutate(genus = str_replace_na(genus, "g_unclassified")) %>%
  mutate(species = str_replace_na(species, "s__unclassified")) %>%
  mutate(ASV = str_c("ASV__",feature_ID)) 
  
    
# metadata
# asv_counts
# taxonomy

asv_rel_abund <-
  inner_join(metadata, asv_counts, by = c("sample_ID"="sample_ID")) %>%
  inner_join(.,taxonomy, by =c("feature_ID"="feature_ID")) %>%
  group_by(sample_ID) %>%
  mutate(rel_abund = 100 * count/sum(count)) %>%
  mutate(sum_rel_abund = sum(rel_abund))%>%
  ungroup()%>%
  select(-count) %>%
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus", "species","ASV"),
               names_to = "level",
               values_to = "taxon")

asv_rel_abund %>% arrange(sum_rel_abund)
asv_rel_abund %>% arrange(desc(sum_rel_abund))

##1. phylum composition facet by anatomical location and sex
phylum_plot <-
  asv_rel_abund %>%
  filter(level=="phylum") %>%
  filter(excluded == "no") %>%
  group_by(sex,Anatomical_Location,DiseaseStatusOracle_abbrev,sample_ID,taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(#sex,
           Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
  summarise(mean_rel_abund = mean(rel_abund), .groups = "drop") %>%
  # mutate(Anatomical_Location = factor(Anatomical_Location))
  mutate(taxon =str_remove_all(taxon,"^\\w__")) %>%
  mutate(taxon = str_replace_all(taxon, "^(\\S*)$", "*\\1*")) %>%
  mutate(taxon =fct_reorder(taxon, mean_rel_abund,.desc = T)) %>%
  ggplot(aes(x=DiseaseStatusOracle_abbrev, y=mean_rel_abund, fill=taxon)) + 
  geom_col() +
  # facet_wrap(~ Anatomical_Location + sex, nrow = 3,
  facet_wrap(~ Anatomical_Location, nrow = 3,
             strip.position = "top") +
  labs(x=NULL,y=NULL
         # "Mean Relative Abundance (%)"
       ) +
  scale_y_continuous(expand = c(0,0),n.breaks = 10)+
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_manual(name="Mean Relative Abundance (%)<br>of Phylum", #"Human Oral<br>Weighted<br>SILVA 138.1",
                    values= c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#8dd3c7','#fccde5'))+
  
  theme_classic() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_markdown(face = "bold"),
        legend.text = element_markdown(face = "bold"),
        legend.justification = "top",
        strip.text.x.top = element_text(face = "bold", size = 9),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.y.left = element_text(face = "bold", size = 16))

phylum_plot


# genus: see separate genus level bar plot.

# taxonomy_genus_uncultured <-
#   taxonomy %>%
#   filter(str_detect(genus, "uncultured")) %>%
#   unite("genus", c("phylum", "order","family"), sep = "<br>", remove = F) %>%
#   mutate(genus = str_c("uncultured_",genus))
# 
# 
# taxonomy_genus_curated <-
#   taxonomy %>%
#   filter(!str_detect(genus, "uncultured")) %>%
#   bind_rows(.,taxonomy_genus_uncultured) 
# 
# genus_rel_abund <-
#   inner_join(metadata, asv_counts, by = c("sample_ID"="sample_ID")) %>%
#   inner_join(.,taxonomy_genus_curated, by =c("feature_ID"="feature_ID")) %>%
#   group_by(sample_ID) %>%
#   mutate(rel_abund = 100 * count/sum(count)) %>%
#   ungroup()%>%
#   select(-count) %>%
#   pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
#                         "family", "genus", "species","ASV"),
#                names_to = "level",
#                values_to = "taxon") %>%
#   filter(level=="genus") %>%
#   filter(excluded == "no") %>%
#   group_by(sex,Anatomical_Location,DiseaseStatusOracle_abbrev,sample_ID,taxon) %>%
#   summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
#   group_by(sex,Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
#   summarise(mean_rel_abund = mean(rel_abund), .groups = "drop") %>%
#   mutate(taxon =str_remove_all(taxon,"^\\w__")) %>%
#   mutate(taxon = str_replace_all(taxon, "^(\\S*)$", "*\\1*")) 
# 
# genus_pool <-
#   genus_rel_abund %>%
#   group_by(taxon) %>% #summarize(max = max(mean_rel_abund)) %>% arrange(max) %>% filter(max>10)
#   summarise(pool = max(mean_rel_abund) < 3, .groups = "drop")
# 
# genus_plot <-
#   inner_join(genus_rel_abund, genus_pool, by = "taxon") %>%
#   mutate(taxon = if_else(pool, "Other (< 3%)", taxon)) %>%
#   group_by(sex,Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
#   summarise(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>%
#   # mutate(taxon =fct_reorder(taxon, mean_rel_abund,.desc = T)) %>% #arrange(desc(taxon))
#   # mutate(taxon =fct_shift(taxon,n = 2)) %>%
#   ggplot(aes(x=DiseaseStatusOracle_abbrev, y=mean_rel_abund, fill=taxon,
#              color=NULL)) +
#   geom_col() +
#   facet_wrap(~Anatomical_Location, nrow = 3) +
#   labs(x=NULL,y="Mean Relative Abundance (%)") +
#   scale_fill_manual(name="Human Oral<br>Weighted<br>SILVA 138.1",
#                     values= c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
#                               '#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#ccebc5','#ffed6f','#bc80bd',
#                               # '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30',
#                               # '#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a',
#                               "darkgray"))+
#   
#   scale_y_continuous(expand = c(0,0),n.breaks = 10)+
#   theme_classic() +
#   guides(fill = guide_legend(ncol = 1)) +
#   theme(legend.title = element_markdown(face = "bold"),
#         legend.text = element_markdown(face = "bold"),
#         legend.justification = "top",
#         strip.text.x.top = element_text(face = "bold", size = 12),
#         axis.text.y = element_text(face = "bold", size = 12, color = "black"),
#         axis.text.x = element_text(face = "bold", size = 12, color = "black",
#                                    angle = 45, hjust = 1, vjust = 1),
#         axis.title.y.left = element_text(face = "bold", size = 16))
# 
# 
# genus_plot


##3. species composition facet by anatomical location and sex 

taxonomy_species_unclassified <-
  taxonomy %>%
  filter(str_detect(species,"s__unclassified$")) %>%
  unite("species",c("phylum","class","order","family","genus"), sep = "<br>",remove = F) %>%
  mutate(species = str_c("unclassified_",species))

taxonomy_species_uncultured <-
  taxonomy %>%
  filter(str_detect(species, "uncultured")) %>%
  unite("speices", c("phylum", "class","order","family","genus"), sep = "<br>", remove = F) %>%
  mutate(species = str_c("uncultured_",speices))

taxonomy_species_unidentified <-
  taxonomy %>%
  filter(str_detect(species, "unidentified")) %>%
  unite("speices", c("phylum","class","order","family", "genus"), sep = "<br>", remove = F) %>%
  mutate(species = str_c("unidentified",speices))

taxonomy_species_curated <-
  taxonomy %>%
  filter(!str_detect(species,"s__unclassified$")) %>%
  filter(!str_detect(species, "uncultured")) %>%
  filter(!str_detect(species, "unidentified")) %>%
  bind_rows(.,
            taxonomy_species_unclassified, 
            taxonomy_species_uncultured,
            taxonomy_species_unidentified) 

species_rel_abund <-
  inner_join(metadata, asv_counts, by = c("sample_ID"="sample_ID")) %>%
  inner_join(.,taxonomy_species_curated, by =c("feature_ID"="feature_ID")) %>%
  group_by(sample_ID) %>%
  mutate(rel_abund = 100 * count/sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(cols = c("kingdom", "phylum", "class", "order", 
                        "family", "genus", "species","ASV"),
               names_to = "level",
               values_to = "taxon") %>%
  filter(level=="species") %>%
  filter(excluded == "no") %>% 
  group_by(#sex,
           Anatomical_Location,DiseaseStatusOracle_abbrev,sample_ID,taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(#sex,
           Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
  summarise(mean_rel_abund = mean(rel_abund), .groups = "drop") %>%
  mutate(taxon =str_remove_all(taxon,"^\\w__")) %>%
  mutate(taxon = str_replace_all(taxon, "^(\\S*)$", "*\\1*")) %>%
    mutate(taxon = str_replace_all(taxon,"(.*)_(.*)", "\\1 \\2"))

species_pool <-
  species_rel_abund %>%
  group_by(taxon) %>% #summarise(max = max(mean_rel_abund)) %>% arrange(desc(max)) %>% pull
  summarise(pool = max(mean_rel_abund) < 3, .groups = "drop")




species_plot <-
  inner_join(species_rel_abund, species_pool, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other (< 3%)", taxon)) %>%
  group_by(#sex,
    Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>% #distinct(taxon) %>% pull
  mutate(taxon=case_when(taxon=="*unclassified_p__Bacteroidota<br>c__Bacteroidia<br>o__Bacteroidales<br>f__Prevotellaceae<br>g__Prevotella 7*"~"*Prevotella (unclassified)*",
                         taxon=="*unclassified_p__Firmicutes<br>c__Bacilli<br>o__Lactobacillales<br>f__Streptococcaceae<br>g_ Streptococcus*"~"*Streptococcus (unclassified)*",
                         taxon=="*unclassified_p__Firmicutes<br>c__Clostridia<br>o__Peptostreptococcales-Tissierellales<br>f__Family_XI<br>g_ Peptoniphilus*"~"*Peptoniphilus (unclassified)*",
                         taxon=="*unclassified_p__Firmicutes<br>c__Negativicutes<br>o__Veillonellales-Selenomonadales<br>f__Veillonellaceae<br>g_ Veillonella*"~"*Veillonella (unclassified)*",
                         taxon=="*unclassified_p__Fusobacteriota<br>c__Fusobacteriia<br>o__Fusobacteriales<br>f__Fusobacteriaceae<br>g_ Fusobacterium*"~"*Fusobacterium (unclassified)*",
                         taxon=="*uncultured_p__Firmicutes<br>c__Bacilli<br>o__Staphylococcales<br>f__Gemellaceae<br>g_ Gemella*"~"*Gemella (unclassified)*",
                         taxon=="*uncultured_p__Firmicutes<br>c__Clostridia<br>o__Peptostreptococcales-Tissierellales<br>f__Family_XI<br>g_ Finegoldia*"~"*Finegoldia (unclassified)*",
                         taxon=="*uncultured_p__Proteobacteria<br>c__Gammaproteobacteria<br>o__Burkholderiales<br>f__Neisseriaceae<br>g_ uncultured*"~"*Neisseriaceae (unclassified)*",
                         .default = as.character(taxon)))  %>%  # distinct(taxon) %>% pull
  
  # mutate(taxon= str_replace(taxon, "c__\\w<br>g__", "g__")) %>%
  # mutate(taxon =fct_reorder(taxon, mean_rel_abund,.desc = T)) %>% #arrange(desc(taxon))
  # mutate(taxon =fct_shift(taxon,n = 1)) %>%
  
  ggplot(aes(x=DiseaseStatusOracle_abbrev, y=mean_rel_abund, fill=taxon,
             color=NULL)) +
  geom_col() +
  facet_wrap(~Anatomical_Location, nrow = 3) +
  # facet_wrap(~Anatomical_Location+sex, nrow = 3) +
  labs(x=NULL,y=NULL
         # "Mean Relative Abundance (%)"
       ) +
  scale_fill_manual(name="Mean Relative Abundance (%)<br>of Species", #"Human Oral<br>Weighted<br>SILVA 138.1",
                    values= c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#1f78f4','#ffff99','#b15928',
                              '#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69','#ffed6f','#bf812d',
                              '#543005','#8c510a','#bc80bd','#fccde5','#d9d9d9','#ccebc5','#dfc27d','#f6e8c3','#f5f5f5',#'#c7eae5','#80cdc1','#35978f','#01665e','#003c30',
                              # '#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a',
                              "darkgray"))+
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0),n.breaks = 10)+
  theme_classic() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_markdown(face = "bold"),
        legend.text = element_markdown(face = "bold"),
        legend.justification = "top",
        # legend.key.size = unit(10,"pt"),
        strip.text.x.top = element_text(face = "bold", size = 9),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.y.left = element_text(face = "bold", size = 16))

species_plot


tailored_species_relabund <-
inner_join(species_rel_abund, species_pool, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other (< 3%)", taxon)) %>%
  group_by(#sex,
           Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
  summarise(mean_rel_abund = sum(mean_rel_abund), .groups = "drop") %>% #distinct(taxon) %>% pull
  mutate(taxon=case_when(taxon=="*unclassified_p__Bacteroidota<br>c__Bacteroidia<br>o__Bacteroidales<br>f__Prevotellaceae<br>g__Prevotella 7*"~"*Prevotella (unclassified)*",
                         taxon=="*unclassified_p__Firmicutes<br>c__Bacilli<br>o__Lactobacillales<br>f__Streptococcaceae<br>g_ Streptococcus*"~"*Streptococcus (unclassified)*",
                         taxon=="*unclassified_p__Firmicutes<br>c__Clostridia<br>o__Peptostreptococcales-Tissierellales<br>f__Family_XI<br>g_ Peptoniphilus*"~"*Peptoniphilus (unclassified)*",
                         taxon=="*unclassified_p__Firmicutes<br>c__Negativicutes<br>o__Veillonellales-Selenomonadales<br>f__Veillonellaceae<br>g_ Veillonella*"~"*Veillonella (unclassified)*",
                         taxon=="*unclassified_p__Fusobacteriota<br>c__Fusobacteriia<br>o__Fusobacteriales<br>f__Fusobacteriaceae<br>g_ Fusobacterium*"~"*Fusobacterium (unclassified)*",
                         taxon=="*uncultured_p__Firmicutes<br>c__Bacilli<br>o__Staphylococcales<br>f__Gemellaceae<br>g_ Gemella*"~"*Gemella (unclassified)*",
                         taxon=="*uncultured_p__Firmicutes<br>c__Clostridia<br>o__Peptostreptococcales-Tissierellales<br>f__Family_XI<br>g_ Finegoldia*"~"*Finegoldia (unclassified)*",
                         taxon=="*uncultured_p__Proteobacteria<br>c__Gammaproteobacteria<br>o__Burkholderiales<br>f__Neisseriaceae<br>g_ uncultured*"~"*Neisseriaceae (unclassified)*",
                         .default = as.character(taxon)))

tailored_species_relabund %>% 
  group_by(Anatomical_Location, taxon) %>%
  summarise(mean_rel_abund=mean(mean_rel_abund), .groups = "drop") %>%
  unique() %>%
  write_tsv("../06-final-barplots/oral-silva138_1/L7- taxa mean relative abundance by anatomical site.tsv")


tailored_species_relabund %>%
  group_by(Anatomical_Location, DiseaseStatusOracle_abbrev, taxon) %>%
  summarise(mean_rel_abund=mean(mean_rel_abund), .groups = "drop") %>%
  unique() %>%
  write_tsv("../06-final-barplots/oral-silva138_1/L7- taxa mean relative abundance of MS diagnosis in anatomical site.tsv")

tailored_phyla_relabund <-
  asv_rel_abund %>%
  filter(level=="phylum") %>%
  filter(excluded == "no") %>%
  group_by(#sex,
           Anatomical_Location,DiseaseStatusOracle_abbrev,sample_ID,taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
  summarise(mean_rel_abund = mean(rel_abund), .groups = "drop") %>%
  mutate(taxon =str_remove_all(taxon,"^\\w__")) %>%
  mutate(taxon = str_replace_all(taxon, "^(\\S*)$", "*\\1*")) 



tailored_phyla_relabund %>%
  group_by(Anatomical_Location, taxon) %>%
  summarise(mean_rel_abund=mean(mean_rel_abund),.groups = "drop") %>%
  unique() %>%
  write_tsv("../06-final-barplots/oral-silva138_1/L2- taxa mean relative abundance by anatomical site.tsv")


tailored_phyla_relabund %>% 
  group_by(Anatomical_Location, DiseaseStatusOracle_abbrev, taxon) %>%
  summarise(mean_rel_abund=mean(mean_rel_abund), .groups = "drop") %>%
  unique() %>%
  write_tsv("../06-final-barplots/oral-silva138_1/L2- taxa mean relative abundance of MS diagnosis in anatomical site.tsv")














# EXPORT
setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/qiime_v202309/07-taxa/06-final-barplots")
list.files()

ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1.pdf", 
       plot = phylum_plot , width = 5, height = 9)
ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1.svg", 
       plot = phylum_plot , width = 5, height = 9)

ggsave(filename =  "oral-silva138_1/L6-taxa-barchart-w-oral-silva138_1-other3.pdf", 
       plot = genus_plot , width = 5, height = 9)
ggsave(filename =  "oral-silva138_1/L6-taxa-barchart-w-oral-silva138_1.svg", 
       plot = genus_plot , width = 5, height = 9)

ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1-other3.pdf", 
       plot = species_plot , width = 8, height = 9)
ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1.svg", 
       plot = species_plot , width = 8, height = 12)

ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1-other3_v2.pdf", 
       plot = species_plot , width = 6, height = 9)
ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1_v2.svg", 
       plot = species_plot , width = 6, height = 12)
ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1_v2.pdf", 
       plot = phylum_plot , width = 5, height = 9)
ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1_v2.svg", 
       plot = phylum_plot , width = 5, height = 9)

# changed facet setting with sex cancelled in the plots
ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1_v2.pdf", 
       plot = phylum_plot , width = 4, height = 9)
ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1_v2.svg", 
       plot = phylum_plot , width = 4.6, height = 9)

ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1-other3_v3.pdf", 
       plot = species_plot , width = 4.6, height = 9)
ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1_v3.svg", 
       plot = species_plot , width = 4.6, height = 12)

# Export for eBioMedicine Revision v4

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/03b-EbioMedicine-Revise/")
dir.create("Adjusted_Figures/")
list.files()

ggsave(filename =  "02-Adjusted_Figures_Tables/L2-taxa-barchart-w-oral-silva138_1-v4.pdf", 
       plot = phylum_plot , width = 4, height = 9)
# ggsave(filename =  "oral-silva138_1/L2-taxa-barchart-w-oral-silva138_1_v2.svg", 
#        plot = phylum_plot , width = 4.6, height = 9)



ggsave(filename =  "02-Adjusted_Figures_Tables/L7-taxa-barchart-w-oral-silva138_1-other3_v4.pdf", 
       plot = species_plot , width = 4.6, height = 10.4)
# ggsave(filename =  "oral-silva138_1/L7-taxa-barchart-w-oral-silva138_1_v3.svg", 
#        plot = species_plot , width = 4.6, height = 12)


#############

phylum_plot_dodged <-
  asv_rel_abund %>%
  filter(level=="phylum") %>%
  filter(excluded == "no") %>%
  group_by(sex,Anatomical_Location,DiseaseStatusOracle_abbrev,sample_ID,taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups = "drop") %>%
  group_by(sex,Anatomical_Location,DiseaseStatusOracle_abbrev,taxon) %>%
  summarise(mean_rel_abund = mean(rel_abund), .groups = "drop") %>%
  mutate(taxon =str_remove_all(taxon,"^\\w__")) %>%
  mutate(taxon = str_replace_all(taxon, "^(\\S*)$", "*\\1*")) %>%
  
  # filter(sex =="Males") %>%
  mutate(taxon =fct_reorder(taxon, mean_rel_abund,.desc = T)) %>%
  ggplot(aes(x=DiseaseStatusOracle_abbrev, y= mean_rel_abund, fill= taxon)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~Anatomical_Location+taxon, nrow = 3) +
  labs(title = NULL,
       x=NULL,y="Mean Relative Abundance (%)") +
  scale_x_discrete()+
  scale_fill_manual(name="Human Oral<br>Weighted<br>SILVA 138.1",
                    values= c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',
                              '#8dd3c7','#bebada','#fb8072','#80b1d3',"darkgray",'#fdb462','#b3de69','#fccde5','#d9d9d9','#ccebc5','#ffed6f','#bc80bd',
                              # '#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30',
                              # '#e0e0e0','#bababa','#878787','#4d4d4d','#1a1a1a',
                              "darkgray"))+
  
  scale_y_continuous(limits = c(0,100),
                     expand = c(0,0),
                     n.breaks = 10)+
  theme_classic() +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.title = element_markdown(face = "bold"),
        legend.text = element_markdown(face = "bold"),
        legend.justification = "top",
        # legend.key.size = unit(10,"pt"),
        strip.text.x.top = element_markdown(face = "bold", size = 8),
        axis.text.y = element_text(face = "bold", size = 12, color = "black"),
        axis.text.x = element_text(face = "bold", size = 12, color = "black",
                                   angle = 45, hjust = 1, vjust = 1),
        axis.title.y.left = element_text(face = "bold", size = 16))

ggsave(filename =  "oral-silva138_1/L2-dodged-taxa-barchart-w-oral-silva138_1.pdf", 
       plot = phylum_plot_dodged ,
       width = 15, height = 9)
ggsave(filename =  "oral-silva138_1/L2-dodged-taxa-barchart-w-oral-silva138_1.svg", 
       # plot = phylum_plot , 
       width = 15, height = 3)

  
################ prevalence control, cut off value to remove/retain feature in samples: 10%



n_sample <- inner_join(metadata, asv_counts, by = c("sample_ID"="sample_ID")) %>%
  select(Anatomical_Location, DiseaseStatusOracle_abbrev, sample_ID, excluded) %>%
  filter(excluded == "no") %>%
  select(-excluded) %>%
  distinct() %>% 
  unite("Groups", c("Anatomical_Location", "DiseaseStatusOracle_abbrev"), remove = T) %>%
  group_by(Groups) %>%
  summarise(n_sample_per_group = n(), .groups = "drop")

ASV_prevalence_remove <-
  inner_join(metadata, asv_counts, by = c("sample_ID"="sample_ID")) %>%
  filter(str_detect(excluded, "no")) %>%
  select(-excluded) %>%
  select(sample_ID, Anatomical_Location, DiseaseStatusOracle_abbrev,feature_ID,count) %>%
  unite("Groups", c("Anatomical_Location", "DiseaseStatusOracle_abbrev"), remove = T) %>%
  group_by(Groups, sample_ID, feature_ID) %>%
  mutate(presence = if_else(count > 0, 1, 0)) %>% 
  ungroup() %>%
  group_by(Groups, feature_ID) %>%
  mutate(n_sample_w_feature_present_per_group = sum(presence)) %>%
  ungroup() %>%
  left_join(., n_sample, by = "Groups") %>%
  group_by(Groups, feature_ID) %>%
  summarise(prevalence = 100 * (n_sample_w_feature_present_per_group/n_sample_per_group), .groups = "drop") %>%
  group_by(feature_ID) %>%
  summarise(remove = max(prevalence) < 10, .groups = "drop")

maaslin_input_taxa <-
  inner_join(asv_rel_abund, ASV_prevalence_remove, by = "feature_ID") %>%
  filter(remove=="FALSE",
         excluded =="no") %>%
  select(feature_ID, sample_ID, rel_abund, level, taxon) %>%
  pivot_wider(., names_from = "level", values_from = "taxon",values_fill = "NULL") %>%
  unite("feature_name", c("kingdom","phylum","class","order","family", "genus", "species", "ASV"), sep = ";", remove = T) %>%
  select(-feature_ID) %>%
  pivot_wider(., names_from = feature_name, values_from = rel_abund, values_fill = 0) 

setwd("~/Dropbox (Partners HealthCare)/CL-04-NOMMS/qiime_v202309/09-MaAsLin")
write_tsv(x=maaslin_input_taxa,
          file = "input_ASV_rel_abund.tsv")

# library(Maaslin2)
# input_data = system.file("extdata", "HMP2_taxonomy.tsv", package="Maaslin2") # The abundance table file
# input_data
# df_input_data = read.table(file             = input_data,
#                            header           = TRUE,
#                            sep              = "\t", 
#                            row.names        = 1,
#                            stringsAsFactors = FALSE)













