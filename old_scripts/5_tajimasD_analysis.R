# tajima's D analysis

# libraries

library(dplyr)
library(ggplot2)
library(magrittr)

#trust chrom to num function

chrom.to.num<-function(x){
  chrom.rom<-as.character(as.roman(c(1:21)))
  return(match(x,chrom.rom))
}

# get the taj d files

file.folder <- "ev_prep_scripts/tajimasD"
file.remove(file.path(getwd(),"ev_prep_scripts/tajimasD/tajimasD.combined.D"))
tajd.files <- list.files(file.folder, full.names = TRUE)

# bind the files together

tajima.list <- list()

for (i in 1:length(tajd.files)){
  
  # get study location and ecotype info
  file.current <- read.table(tajd.files[i], header = TRUE)
  study<- strsplit(tajd.files[i],split="/")[[1]][3] %>% strsplit(split="\\.") %>% .[[1]] %>% .[1]
  location <- strsplit(tajd.files[i],split="/")[[1]][3] %>% strsplit(split="\\.") %>% .[[1]] %>% .[2]
  location <- strsplit(location,split="_")[[1]][1]
  ecotype <- strsplit(tajd.files[i],split="/")[[1]][3] %>% strsplit(split="\\.") %>% unlist %>% .[2] %>% strsplit(split="_") %>% unlist
  
  if (length(ecotype) > 1){
    ecotype <- ecotype[2]
  }
  
  #rm scaffolds
  
  file.current <- file.current %>%
    filter(!grepl(CHROM, pattern = "scaff"))
  
  # chrom to num conversion
  
  file.current$CHROM <- chrom.to.num(gsub("group", "", file.current$CHROM))
  
  # add in window start pos
  
  file.current$pos1 <- file.current$BIN_START - 75000
  file.current$pos2 <- file.current$BIN_START - 1
  
  # format for binding
  tajima.list[[i]] <- data.frame(study = study, location = location, ecotype = ecotype, lg = file.current$CHROM, pos1 = file.current$pos1, pos2 = file.current$pos2, n.sites = file.current$N_SNPS , tajd = file.current$TajimaD)
}

tajd.out <- do.call("rbind", tajima.list)
tajd.out <- tajd.out[complete.cases(tajd.out),]

write.table(tajd.out, file = "ev_prep_scripts/tajimasD/tajimasD.combined.D")

