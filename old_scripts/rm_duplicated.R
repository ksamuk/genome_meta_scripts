
# analysis files
library(dplyr)

out.folder <- "analysis_ready/clustering_fst"
out.files <- list.files(out.folder)

out.files.exist <- out.files %>% gsub("\\d","",.) %>% gsub("-","",.) %>% gsub(".gz_.clustered.txt",".txt.gz",.)
dupes <- duplicated(out.files.exist) %>% which %>% out.files.exist[.]

out.files <- list.files(out.folder, full.names = TRUE)

for (i in dupes){
  dupe.list <- out.files[grep(i,out.files.exist)]
  print(paste0("removing ", dupe.list[2],"..."))
  file.remove(dupe.list[2])
}
