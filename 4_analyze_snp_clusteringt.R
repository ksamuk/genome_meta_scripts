### process clustering script output

library("dplyr")

cluster.folder <- "analysis_ready/clustering_tmp"
cluster.files <- list.files(cluster.folder, full.names = TRUE)

header.row <- names (read.table(file = cluster.files[1], header = TRUE, stringsAsFactors = FALSE))


read_files_chopped <- function(cluster.file){
  return(read.table(file = cluster.files[1], stringsAsFactors = FALSE, skip = 1))
}

file.chopped <- lapply(cluster.files, read_files_chopped)
cluster.df <- bind_rows(file.chopped)
names(cluster.df) <-  header.row 
