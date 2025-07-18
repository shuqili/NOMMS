library(tidyverse)
library(ggtext)

getwd()
setwd("09-MaAsLin/03-output/")

nasal <-
  read_tsv("01-nostril/ASVlevel/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_location = "N")

gin <-
  read_tsv("02-gums/ASVlevel/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_location = "G")

oro <-
  read_tsv("03-throat/ASVlevel/tobacco_exposure + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  mutate(anatomical_location = "O")


# variables and maps
asv.abbrev <-
read_tsv("../04-blastn-significant-ASVs/ASV.abbreviation.list.v3.tsv") %>%
  select(number, feature)
p = 0.05
q = 0.05

l8.taxa.check.asv.abbrev <-
bind_rows(nasal, gin, oro) %>%
  filter(metadata == "tobacco_exposure") %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_location, levels = c("N", "G", "O"))) %>% 
  select(feature) %>% 
  mutate(feature = str_remove(feature, "d__.*g__"))  %>% 
  separate(feature, into = c("feature", "feature_id"), sep = ".ASV__") %>%
  left_join(., asv.abbrev, by = c("feature_id" = "feature")) 

asv.abbrev.added.jun12 <-
l8.taxa.check.asv.abbrev %>%
  filter(is.na(number)) %>% 
  select(number, feature_id) %>%
  rename(feature = feature_id) %>%
  arrange(feature) %>% # critical step for ASV abbrev match.
  mutate(number= seq(129,138,1)) %>%
  mutate(number = str_replace(number, "^(.*$)", "ASV_\\1")) 

asv.abbrev.added.jun12 %>% 
  bind_rows(asv.abbrev,.) %>%
  write_tsv("../04-blastn-significant-ASVs/ASV.abbreviation.list.v4.tsv")

getwd()
# Blast the added sig ASVs, Export as tsv, Upload the BLASTn server, Run and save results (hit table csv) for later data parsing and reformatting.
# Check the bottom of this document to see the detailed settings on BLASTn server
feature.list <-
  read_tsv("../../03-features/oral-silva138_1/sequences.tsv",col_names = F) %>%
  mutate(separate = if_else(str_detect(X1,">"),"feature","sequences")) %>%
  mutate(X1 = str_remove(X1, ">"), entry = rep(c(1:916), each = 2)) %>% # tibble row number 1832/2=916
  select(entry, separate, X1) %>%
  pivot_wider(., names_from = separate, values_from = X1) %>%
  select(-entry)
Jun12.blastn.pool <-
  left_join(asv.abbrev.added.jun12, feature.list, by = "feature") %>%
  select(number, sequences) %>%
  mutate(number = str_replace(number, "^(.*)$", ">\\1")) %>%
  mutate(entry = c(1:10)) %>% # use a new column as a anchor for pivot longer the table and match contents.
  select(entry, number, sequences) %>%
  pivot_longer(cols = number:sequences, names_to = "name", values_to = "value") %>%
  select(value)

Jun12.blastn.pool %>%
  write_tsv("../04-blastn-significant-ASVs/Jun12.blastn.pool.tsv",col_names = F)


# Plot l8

asv.new.list <-
  read_tsv("../04-blastn-significant-ASVs/ASV.abbreviation.list.v4.tsv")

l8.taxa <-
bind_rows(nasal, gin, oro) %>%
  filter(metadata == "tobacco_exposure") %>%
  filter(pval < p) %>%
  mutate(anatomical_site = factor(anatomical_location, levels = c("N", "G", "O"))) %>%
  separate(feature, into = c("feature", "feature_id"), sep = ".ASV__") %>% 
  left_join(., asv.new.list, by = c("feature_id" = "feature")) %>% #
  select(feature, number, everything()) %>%
  mutate(feature = str_remove(feature, "d__.*g__"))  %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***")) %>%
  mutate(feature = if_else(str_detect(feature, "s__un"), 
                           str_replace(feature, "^(.*).s__(.*)", "Unclassified \\1***"), 
                           str_replace(feature, "^(.*).s__(.*)_(.*)", "***\\2 \\3"))) %>%# print(., n= 30)
  mutate(feature = if_else(feature=="Unclassified ***d__Bacteria.p__Proteobacteria.c__Gammaproteobacteria.o__Enterobacterales.f__Pasteurellaceae.g_unclassified***",
                           "Unclassified ***Pasteurellaceae***",
                           feature))%>% 
  unite("feature", feature:number, sep = " ")
  
 
#
l8.taxa %>% # view to see max/min of coef
  select(feature, anatomical_site, everything()) %>%
  unite("feature.site", feature:anatomical_site, remove = F, sep = ":") %>%
  group_by(anatomical_site) %>%
  arrange(-coef, .by_group = T) %>%
  ungroup() %>%
  mutate(feature_pos = c(1:30)) %>%
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
                       limits=c(-4,4),
                       breaks = c(-4, 0, 4),
                       labels = c("-4","0","4")) +
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

ggsave("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/10-tobacco_exposure/tobacco_exposure-l8.pdf",
       width = 3.2, height = 7.2, dpi = 300)




#2. Upload the tsv to https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch#
# parameters setup:
# Set 1 - std. db. nr /OR/ Set 2 - 16S rRNA_typestrains
# check the box for ""Uncultured/environmental sample sequences"
# check optimization for "Highly similar sequences (megablast)"
# check the box for "Show results in a new window."
# Hit "BLAST"

# Hit "Download ALL" from the result page, go for the "Hit table (csv)"
# Save and rename to distinguish nr from 16SrRNA results