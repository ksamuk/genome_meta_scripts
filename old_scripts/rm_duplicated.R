
# analysis files
library(dplyr)

out.folder <- "analysis_ready/clustering_tmp"
out.files <- list.files(out.folder)

out.files.exist <- out.files %>% gsub("\\d","",.) %>% gsub("-","",.) %>% gsub(".gz_.clustered.txt",".txt.gz",.)
stats.reformat <- list.files(stats.folder)
duplicated(out.files.exist) %>% which %>% out.files.exist[.]



<- list.files("analysis_ready/clustering_fst") 
