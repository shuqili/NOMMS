library(tidyverse)
library(ggtext)


getwd()
setwd( "~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin")
list.files()

# set constant variables
q = 0.05
p = 0.05
# set asv abbrev map
asv.match <-
read_tsv("04-blastn-significant-ASVs/ASV.abbreviation.list.v5.tsv")
# set input data
nasal.p <-
read_tsv("03-output/01-nostril/ASVlevel/Smoking_YN + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "smoking_YN") %>%
  filter(pval < p) %>%
  mutate(Anatomical_Location = "nasal")

gingival.p <-
read_tsv("03-output/02-gums/ASVlevel/Smoking_YN + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "smoking_YN") %>%
  filter(pval < p) %>%
  mutate(Anatomical_Location = "gingival")

oro.p <-
read_tsv("03-output/03-throat/ASVlevel/Smoking_YN + [age + sex + bmi ]/cplmfit-no-transformation-relab-asv-to-diseasestat/all_results.tsv") %>%
  filter(metadata == "smoking_YN") %>%
  filter(pval < p) %>%
  mutate(Anatomical_Location = "oropharyngeal")

# consolidate input data
l8.taxa <-
bind_rows(nasal.p, gingival.p, oro.p) %>%
  select(feature, Anatomical_Location, metadata, value, everything()) %>%
  rename(anatomical_site = Anatomical_Location) %>%
  mutate(anatomical_site = factor(anatomical_site, levels = c("nasal", "gingival", "oropharyngeal"))) %>% 
  separate(feature, into = c("feature_name", "ASV"), sep = "\\.ASV__", remove = F) %>% 
  left_join(., asv.match, by = c("ASV" = "feature"))%>% # print(., n=28)
  select(-feature_name, -ASV) %>%
  select(feature, number, anatomical_site, everything()) %>%
  mutate(feature = str_remove(feature, "d__.*g__"))  %>%
  mutate(feature = str_remove(feature, "\\.ASV__.*")) %>%
  mutate(feature = str_replace(feature, "(.*)", "***\\1***")) %>%
  mutate(feature = if_else(str_detect(feature, "s__un"), 
                           str_replace(feature, "^(.*).s__(.*)", "Unclassified \\1***"), 
                           str_replace(feature, "^(.*).s__(.*)_(.*)", "***\\2 \\3"))) %>% #print(n=28)
  mutate(feature = if_else(feature=="***.Eubacterium. brachy***",
                           "***Eubacterium brachy***",
                           feature)) %>% # print(n=28)
  unite("feature", feature:number, sep = " ")
  
l8.taxa %>%
    mutate(plabel= if_else(pval<0.05,"*","")) %>%
    mutate(qlabel=if_else(qval<0.05,"+","")) %>%
    unite("label", plabel:qlabel, sep = " ") %>%
    group_by(anatomical_site) %>%
    arrange(desc(coef),.by_group =T) %>%
    ungroup() %>%
    mutate(feature_pos = c(1:28)) %>%

  mutate(feature=fct_reorder(feature,feature_pos,.desc = F)) %>%
  complete(anatomical_site,feature) %>%
  ggplot(aes(x=anatomical_site, y= feature, fill=coef))+
  geom_tile(color="black", linewidth = 0.2) +
  geom_text(aes(label=label)) +
  scale_y_discrete(expand = c(0,0),position = "right") +
  scale_x_discrete(expand = c(0,0), position = "top",
                   breaks = c("nasal", "gingival","oropharyngeal"),
                   # labels = c("N", "G", "O"),
                   labels =c("Nasal", "Gingival","Oropharyngeal")) +
  scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with",
                       colours = c("#0433FF","white","#FF2600" ),
                       na.value = "transparent",
                       limits=c(-2.7,2.7),
                       breaks = c(-2.7, 0, 2.7),
                       labels = c("-2.7","0","2.7")) +
  theme_bw() +
  theme(plot.title = element_markdown(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_markdown(color = "black", face = "bold",
                                       angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_markdown(color = "black"),
        legend.title = element_markdown(face = "bold", hjust = 0),
        legend.text = element_text(face = "bold"),
        legend.position = "left",
        legend.justification = "top")  

getwd()
ggsave("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/05-smokingYN/smoking_YN.age_sex_BMI-l8.pdf",
       dpi = 300, width = 6, height = 7.6)





# # l6.taxa %>%  
# #   unite("feature.site", feature:anatomical_site, remove = F, sep = ":") %>%
# #   group_by(anatomical_site) %>%
# #   arrange(-coef, .by_group = T) %>%
# #   ungroup() %>%
# #   mutate(feature_pos = c(1:10)) %>%
# #   mutate(plabel= if_else(pval<0.05,"*","")) %>%
# #   mutate(qlabel=if_else(qval<0.05,"+","")) %>%
# #   unite("label", plabel:qlabel, sep = " ") %>% 
# #   complete(anatomical_site, feature) %>%
# #   mutate(feature = fct_reorder(feature, feature_pos, .desc = F)) %>%
# #   ggplot(aes(x=anatomical_site,y=feature, fill=coef)) +
# #   geom_tile(color="black") +
# #   geom_text(aes(label=label)) +
# #   scale_y_discrete(expand = c(0,0),position = "right") +
# #   scale_x_discrete(expand = c(0,0),position = "top",
# #                    breaks = c("nasal", "gingival","oropharyngeal"),
# #                    labels = c("N", "G", "O")) +
# #   scale_fill_gradientn(name = NULL, #"Bacteria<br>associated<br>with",
# #                        colours = c("#0433FF","white","#FF2600" ),
# #                        na.value = "transparent",
# #                        limits=c(-3,3),
# #                        breaks = c(-3, 0, 3),
# #                        labels = c("-3","0","3")) +
# #   theme_bw() +
# #   theme(axis.title = element_blank(),
# #         axis.text.x  = element_markdown(face = "bold", color = "black", 
# #                                         angle = 0, vjust = 0.5, hjust=0.5),
# #         axis.text.y = element_markdown(color = "black"),
# #         plot.title = element_markdown(face = "bold"),
# #         panel.grid = element_blank(),
# #         legend.position = "top",
# #         legend.justification = "left",
# #         legend.text = element_text(face = "bold", hjust = 0.5))
# # getwd()
# # ggsave("../../03b-EbioMedicine-Revise/00-github repository/03-differential analyses/06-visualization/05-smokingYN/smoking_YN-l6.pdf",
# #        width = 2.6, height = 4, dpi = 300)
# # 
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# smoking.yn.compilation.step1 <-
# bind_rows(nasal.p, gingival.p) %>%
#   bind_rows(., oro.p) %>% 
#   mutate(Anatomical_Location=factor(Anatomical_Location, 
#                                     labels = c("nasal","gingival","oropharyngeal"),
#                                     levels = c("nasal","gingival","oropharyngeal"))) %>%
#   mutate(plabel= if_else(pval<0.05,"*","")) %>%
#   mutate(qlabel=if_else(qval<0.05,"+","")) %>%
#   unite("label", plabel:qlabel, sep = " ") %>%
#   group_by(Anatomical_Location) %>%
#   arrange(desc(coef),.by_group =T) %>%
#   ungroup() %>%
#   mutate(feature_pos = c(1:28))
# 
# 
# 
# # asv.abbrev.list <-
# # asv.match %>% distinct(ASV.old) %>% pull
# 
#     
# existing.var.smoking.feature <-
#   smoking.yn.compilation.step1 %>% 
#   distinct(feature) %>%
#   mutate(smoking.feature = str_replace(feature,".*ASV__","")) %>%
#   mutate(match = if_else(smoking.feature %in% asv.match, 1, 0)) %>%
#   filter(match ==1) %>%
#     left_join(., asv.match, by = c("smoking.feature"="ASV.old")) %>%
#   select(-match)
# 
# 
# new.var.smoking.feature <-
#   smoking.yn.compilation.step1 %>% 
#   distinct(feature) %>%
#   mutate(smoking.feature = str_replace(feature,".*ASV__","")) %>%
#     mutate(match = if_else(smoking.feature %in% asv.abbrev.list, 1, 0)) %>%
#     filter(match==0) %>% # DESCREPENCY: assigned new asv number without arragne by feature_ID thus caused mis match between sub datasets.
#     mutate(ASV = "ASV",
#            n = c(93:106)) %>%
#     mutate(n=as.character(n)) %>%
#     mutate(ASV.abbrev = str_replace(n,"^(\\d+)$","ASV_\\1")) %>%
#   select(-match,-ASV, -n)
# 
# 
# var.smoking.feature <-
# bind_rows(existing.var.smoking.feature, new.var.smoking.feature) %>%
#   # mutate(feature = str_remove(feature, "ASV(.*)$")) %>%
#   select(feature, ASV.abbrev)
#   # tidyverse::unite("name", feature:ASV.abbrev,sep = "")
# 
# smoking.yn.compilation.step2 <-  
# smoking.yn.compilation.step1 %>%
#   left_join(., var.smoking.feature, by="feature") %>%
#   select(feature, ASV.abbrev, everything()) %>%
#   unite("feature", feature:ASV.abbrev,sep = ";",remove = T) %>%
#   mutate(feature = str_remove(feature, ".*g__")) %>% # distinct(feature) %>% pull
#   mutate(feature=case_when(feature == "Corynebacterium.s__Corynebacterium_kroppenstedtii.ASV__9f6c00ff9e1f455f42f211770d60dec0;ASV_93" ~ "***Corynebacterium kroppenstedtii*** ASV_93",
#                            feature == "Streptococcus.s__Streptococcus_parasanguinis.ASV__75e69ac35f8b5786cf860f3f0be74929;ASV_64" ~ "***Streptococcus parasanguinis*** ASV_64",
#                            feature == "Anaerococcus.s__Anaerococcus_octavius.ASV__3dde63a29d5649dd61c8ce85e2a82f64;ASV_94" ~ "***Anaerococcus octavius*** ASV_94",
#                            feature == "Anaerococcus.s__Anaerococcus_octavius.ASV__44cdc9c06c86362c29a54813aefd2f39;ASV_53" ~ "***Anaerococcus octavius*** ASV_53",
#                            feature == "Porphyromonas.s__unidentified.ASV__644049b2fa1363e0134dae0ab4a1ca50;ASV_95" ~ "***Porphyromonas*** ASV_95",
#                            feature == "Parvimonas.s__unclassified.ASV__6bf2414103f2ae11f0f872d03d711d31;ASV_71" ~ "***Parvimonas*** ASV_71",
#                            feature == "Haemophilus.s__unclassified.ASV__db3d96ce2e755dcb6db0a293e21cc111;ASV_96" ~ "***Haemophilus*** ASV_96",
#                            feature == "Prevotella.s__Prevotella_oris.ASV__fb8ddfb6886f3b08315ca8bbab50a788;ASV_97" ~ "***Prevotella oris*** ASV_97",
#                            feature == ".Eubacterium._brachy_group.s__.Eubacterium._brachy.ASV__de4edf45ab75f6490e1ada65b0e0306f;ASV_75" ~ "***Eubacterium brachy group*** ASV_75",
#                            feature == "Fusobacterium.s__uncultured_bacterium.ASV__4679cc695527888b23df9e041b1e7465;ASV_69" ~ "***Fusobacterium*** ASV_69",
#                            feature == "Fusobacterium.s__unclassified.ASV__c7fe8c5f65f8eaed352dacb33037e7de;ASV_98" ~ "***Fusobacterium*** ASV_98",
#                            feature == "Prevotella_7.s__Prevotella_denticola.ASV__83b33d4e62f2966de49db7cc6c47dfd1;ASV_68" ~ "***Prevotella denticola*** ASV_68",
#                            feature == "Rothia.s__unclassified.ASV__1ce0d1cfa2c68160f70e20616d2d4bc2;ASV_99" ~ "***Rothia*** ASV_99",
#                            feature == "Streptococcus.s__Streptococcus_parasanguinis.ASV__8aa917d28a15681c311a79b1a4954dbe;ASV_80" ~ "***Streptococcus parasanguinis*** ASV_80",
#                            feature == "Haemophilus.s__Haemophilus_parainfluenzae.ASV__5a11fc8118437f874b331f8a288af3d9;ASV_12" ~ "***Haemophilus parainfluenzae*** ASV_12",
#                            feature == "Neisseria.s__Neisseria_perflava.ASV__d5e4fe23604cf247e6ef369390ef0686;ASV_60" ~ "***Neisseria perflava*** ASV_60",
#                            feature == "Lautropia.s__uncultured_bacterium.ASV__f9a668f6c0e2112c51239d0b7c9b9280;ASV_52" ~ "***Lautropia*** ASV_52",
#                            feature == "Oribacterium.s__unidentified.ASV__2a28982df1e6b16f1f7c772d76746276;ASV_23" ~ "***Oribacterium*** ASV_23",
#                            feature == "Alloprevotella.s__uncultured_Bacteroidetes.ASV__28b09b766da0307573913dacd4645fc9;ASV_100" ~ "***Alloprevotella*** ASV_100",
#                            feature == "Treponema.s__Treponema_socranskii.ASV__c5bd025300b4ef0c652dd11cef56ac29;ASV_76" ~ "***Treponema socranskii*** ASV_76",
#                            feature == "Fusobacterium.s__unclassified.ASV__6d547ff67cfe5d9e34df30bdedd05db0;ASV_33" ~ "***Fusobacterium*** ASV_33",
#                            feature == "Actinomyces.s__Actinomyces_graevenitzii.ASV__e7b81c3cb2374da6f1de09a8c954cf7f;ASV_101" ~ "***Actinomyces graevenitzii*** ASV_101",
#                            feature == "Rothia.s__unclassified.ASV__af972921ecfdb7a8dcb9a561ab98a355;ASV_102" ~ "***Rothia*** ASV_102",
#                            feature == "Haemophilus.s__Haemophilus_parainfluenzae.ASV__25be8565bb6d4b0ebcfc117225491be2;ASV_103" ~ "***Haemophilus parainfluenzae*** ASV_103",
#                            feature == "Neisseria.s__uncultured_bacterium.ASV__123249ae7edc1a5664e2a52fb41d171c;ASV_65" ~ "***Neisseria*** ASV_65",
#                            feature == "Haemophilus.s__Haemophilus_parainfluenzae.ASV__e924271260411e885401de5355fc6814;ASV_104" ~ "***Haemophilus parainfluenzae*** ASV_104",
#                            feature == "Leptotrichia.s__Leptotrichia_sp..ASV__bc220f0642fc0da8bfaa703b10337cce;ASV_105" ~ "***Leptotrichia*** ASV_105",
#                            feature == "Prevotella_7.s__Prevotella_melaninogenica.ASV__8daab49003dfc975db2139efcdb6f57a;ASV_106" ~ "***Prevotella melaninogenica*** ASV_106"))
# 
# smoking.yn.compilation.step2 %>%
#   