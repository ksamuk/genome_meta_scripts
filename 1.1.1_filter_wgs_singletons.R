# filter out unique sites from WGS data

library(data.table)
library(dplyr)

#read in an format a cominbed stats file

out.dir <- "/home/sticklegroup/meta/SB_meta_data/fst_dxy/loose_filter/combined"
out.file.name <- list.files(out.dir,pattern="*.txt$", full.names = TRUE)[1]

stats.out <- fread(out.file.name)

print("wrangling stats file...")
stats.out <- stats.out %>%
  mutate(study_com = paste0(study,"_",comparison)) %>%
  mutate(lg_pos = paste0(lg,"_",pos))

# lg_pos's found in the grander dataset
print("finding sites in non-wgs data sets...")
all.sites <- stats.out %>%
  filter(!grepl("japan",study_com))%>%
  filter(!grepl("russian",study_com))%>%
  select(lg_pos) %>%
  unique %>%
  unlist

print("finding sites unique to japanese data...")
japan.sites <- stats.out %>%
  filter(grepl("japan",study_com))%>%
  select(lg_pos) %>%
  unique() %>%
  unlist

print("finding sites unique to russian data...")
russia.sites <- stats.out %>%
  filter(grepl("russian",study_com))%>%
  select(lg_pos) %>%
  unique() %>%
  unlist

print("identifying wgs vs. non-wgs overlap...")
japan.uniques <- japan.sites[!(japan.sites %in% all.sites)]
russia.uniques <- russia.sites[!(russia.sites %in% all.sites)]

print("filtering wgs singletons...")
stats.out <- stats.out[!(stats.out$lg_pos %in% japan.uniques), ]
stats.out <- stats.out[!(stats.out$lg_pos %in% russia.uniques), ]

print("writing to file...")
file.out <- gsub(".txt", ".filt.txt", out.file.name)
write.table(stats.out, file.out, row.names = FALSE)