library(tidyverse)

getwd()
setwd("~/Partners HealthCare Dropbox/Shuqi Li/CL-04-NOMMS/01-NOMMS_Final_qiime_v202309/09-MaAsLin")

#0. Preparing the Pool.
asv.list <- read_tsv("04-blastn-significant-ASVs/ASV.abbreviation.list.v3.tsv")

feature.list <-
  read_tsv("../03-features/oral-silva138_1/sequences.tsv",col_names = F) %>%
  mutate(separate = if_else(str_detect(X1,">"),"feature","sequences")) %>%
  mutate(X1 = str_remove(X1, ">"), entry = rep(c(1:916), each = 2)) %>% # tibble row number 1832/2=916
  select(entry, separate, X1) %>%
  pivot_wider(., names_from = separate, values_from = X1) %>%
  select(-entry)
  
batch.blastn.pool <-
  left_join(asv.list, feature.list, by = "feature") %>%
  select(number, sequences) %>%
  # batch.blastn.pool %>%
  mutate(number = str_replace(number, "^(.*)$", ">\\1")) %>%
  mutate(entry = c(1:128)) %>%
  select(entry, number, sequences) %>%
  pivot_longer(cols = number:sequences, names_to = "name", values_to = "value") %>%
  select(value)

#1. Export
batch.blastn.pool %>%
  write_tsv("batch-blastn-pool.tsv",col_names = F)
  

#2. Upload the tsv to https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch#
# parameters setup:
# Set 1 - std. db. nr /OR/ Set 2 - 16S rRNA_typestrains
# check the box for ""Uncultured/environmental sample sequences"
# check optimization for "Highly similar sequences (megablast)"
# check the box for "Show results in a new window."
# Hit "BLAST"

# Hit "Download ALL" from the result page, go for text file on the very top
# Save and rename to distinguish nr from 16SrRNA results

#3. Read-in results and manually put headers to match with webpage.
##1. typestrain results
res.typestrains <-
  read_tsv("04-blastn-significant-ASVs/YXAPE36M016-Alignment.txt", col_names = F) %>% 
  select(X1) 


pos.start <-
  res.typestrains %>%  
  mutate(entry = c(1:345625)) %>%
  filter(str_detect(X1,"ASV_\\d+")) %>%
  filter(entry != 2) %>%
  pull(entry)
  
# pos.end <-
# res.typestrains %>%  
#   mutate(entry = c(1:345625)) %>%
#   filter(str_detect(X1,"ASV_\\d+")) %>%
#   filter(entry != 2) %>%
#   mutate(entry = entry +104-1) %>%
#   pull(entry)


pos.start
getwd()
res.typestrains  %>%
  mutate(row.number = c(1:345625)) %>%
  mutate(pan.ding = if_else(row.number %in% pos.start, "yes", "no")) %>%
  slice(as.integer(outer(0:25, which(pan.ding == "yes"), `+`))) %>% # extract 26 lines when detected yes (row =0)
  mutate(pan.ding = if_else(str_detect(X1, "ASV_"), "yes", "no")) %>%
  slice(-as.integer(outer(1:3, which(pan.ding == "yes"), `+`))) %>% # remove 3 descriptive rows
  mutate(number = rep(c(1:128), each = 23)) %>% # 128 ASV was sent to blastn.
  mutate(number = str_replace(number,"^(\\d)$","0\\1")) %>% 
  mutate(number = str_replace(number,"^(\\d+)$","ASV_\\1")) %>% # assign ASV number to results
  filter(pan.ding != "yes" ) %>%# remove redundant header lines.
  select(-row.number) %>%
  write_tsv("test.tsv")
  
# Import to Excel by setting column width, check and save as tab-delimited txt.
# read back in R
typestrains.top10 <-
read_tsv("test.rRNA.txt", col_names = F) %>%
  slice(-1) %>%
  select(X1, X7, X8, X9, X10, X11) %>%
  rename(Description = X1, query.cover = X7, e.value = X8, per.ident = X9, acc.len = X10, accession = X11) %>%
  separate(accession, into =c("accession.version", "feature"), sep = "\tno\t") %>%
  mutate(db = "rRNA_typestrains") %>%
  group_by(feature) %>%
  mutate(hit.number = c(1:22)) %>% 
  ungroup() %>%
  select(db, feature, hit.number, everything()) %>%
  filter(hit.number < 11) 

##2. nr results

res.nr <-
  read_tsv("04-blastn-significant-ASVs/YXAMZ68R013-Alignment.txt", col_names = F) %>% 
  select(X1) 


pos.start <-
res.nr %>%  
  mutate(entry = c(1:487634 )) %>%
  filter(str_detect(X1,"ASV_\\d+")) %>%
  filter(entry != 2) %>%
  pull(entry)

pos.start

res.nr %>%
  mutate(row.number = c(1:487634)) %>%
  mutate(pan.ding = if_else(row.number %in% pos.start, "yes", "no")) %>%
  slice(as.integer(outer(0:25, which(pan.ding == "yes"), `+`))) %>% # extract 26 lines when detected yes (row =0)
  mutate(pan.ding = if_else(str_detect(X1, "ASV_"), "yes", "no")) %>%
  slice(-as.integer(outer(1:3, which(pan.ding == "yes"), `+`))) %>% # remove 3 descriptive rows
  mutate(number = rep(c(1:128), each = 23)) %>% # 128 ASV was sent to blastn.
  mutate(number = str_replace(number,"^(\\d)$","0\\1")) %>% 
  mutate(number = str_replace(number,"^(\\d+)$","ASV_\\1")) %>% # assign ASV number to results
  filter(pan.ding != "yes" ) %>%# remove redundant header lines.
  select(-row.number) %>%
  write_tsv("test.tsv")

# Import to Excel by setting column width, check and save as tab-delimited txt. watch closely for the separation.
# read back in R
core.nt.top10 <-
  read_tsv("test.core.nt.txt", col_names = F) %>%
  slice(-1) %>%
  select(X1, X7, X8, X9, X10, X11) %>%
  rename(Description = X1, query.cover = X7, e.value = X8, per.ident = X9, acc.len = X10, accession = X11) %>%
  separate(accession, into =c("accession.version", "feature"), sep = "\tno\t") %>%
  mutate(db = "core_nt") %>%
  group_by(feature) %>%
  mutate(hit.number = c(1:22)) %>% 
  ungroup() %>%
  select(db, feature, hit.number, everything()) %>%
  filter(hit.number < 11) 


# compile typestrain results with core_nt results

final.table <-
bind_rows(core.nt.top10, typestrains.top10) %>%
  select(feature, db, hit.number, everything()) %>%
  mutate(hit = str_c("hit # ", hit.number)) %>% 
  select(feature, db, hit, everything()) %>%
  unite("hit", db:hit, sep = " db ", remove = T) %>% 
  filter(hit.number <4) %>% 
  # write_tsv("prelim.tsv")
  select(-hit.number) %>%
  unite("all", Description:accession.version, sep = "^", remove = T) %>%
  pivot_wider(names_from = hit, values_from = all) %>% 
  separate_wider_delim(cols = -feature, delim = "^", names_sep = ":") %>% 
  rename_with(
    .fn = \(x){
    str_replace_all(string = x, pattern = fixed(":1"), replacement = ":Description") |>
      str_replace_all(pattern = fixed(":2"),replacement = ":query.cover") |>
      str_replace_all(pattern = fixed(":3"),replacement = ":e.value") |>
      str_replace_all(pattern = fixed(":4"),replacement = ":per.ident") |>
      str_replace_all(pattern = fixed(":5"),replacement = ":acc.len") |>
      str_replace_all(pattern = fixed(":6"),replacement = ":accession.version") 
    },
    .cols = -feature)
write_tsv(x = final.table, file = "final-table.tsv")

