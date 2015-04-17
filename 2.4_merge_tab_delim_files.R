# merge tab delim vcf files

library(dplyr)
library(data.table)

# list of files 
file.list <- list.files(pattern="*\\.tab")

# list of chromosomes
chrom.list <- paste(rep("group",21), as.roman(1:21), sep="")

dir.create("merged_chromo", showWarnings = FALSE)

# function to subset by chromo
subsetByChromo <- function(chromo,file){
  lines.df <- system(paste("grep", "-w", chromo,file),intern=TRUE)
  lines.df <- strsplit(lines.df, split="\t")
  lines.df <- data.table(do.call("rbind", lines.df))
  setnames(lines.df,names(read.table(file, nrows=1, header=TRUE)))
  return(lines.df)
}

# init chromo loop

for (i in 1:length(chrom.list)){
  
  print(paste("Merging chromosome",i,"..."))
  print(paste("Reading file ",file.list[1],"...",sep=""))
  file.merged.chrom<-subsetByChromo(chrom.list[i],file.list[1])
  print("Done!")
  
  for (j in 2:length(file.list)){
    
    print(paste("Reading file ",file.list[j],"...",sep=""))
    file.current.chrom<-subsetByChromo(chrom.list[i],file.list[j])
    print("Done!")
    
    print(paste("Merging file ",file.list[j]," with stack...",sep=""))
    file.merged.chrom<-merge(file.merged.chrom,file.current.chrom,by=c("CHROM","POS"),all=TRUE,incomparables="NN")
    print("Done!")
    
    
  }
  print("Chromosome successfully merged! Sorting...")
  file.merged.chrom<-arrange(file.merged.chrom,CHROM,POS)
  
  print(paste("Sorting complete! Writing ",format(object.size(file.merged.chrom),units="Mb")," to file...",sep=""))
  date.stamp<-paste("_",format(Sys.time(),"%b-%d-%Y"),sep="")
  out.file.name<-paste(i,"_chrom_","_all_data",date.stamp,".tab",sep="")
  write.table(file.merged.chrom,file=file.path("merged_chromo",out.file.name),quote=FALSE,row.names=FALSE)
  rm(file.merged.chrom)
}


#### build concatenated file (optional)
library(data.table)
print("All chromosomes merged! Building master file...")
chromo.files <- list.files("merged_chromo")
date.stamp <- paste("_",format(Sys.time(), "%b-%d-%Y"), sep="")
merged.file.name <- paste("sb_meta_all", date.stamp, ".tab",sep="")

chrom.order <- as.numeric(sapply(strsplit(chromo.files,split="_"),function(x)x[1]))
chromo.files <- chromo.files[order(chrom.order)]

header.list<-list()
for (i in 1:length(chromo.files)){
  file.current <- fread(file.path("merged_chromo",chromo.files[i]),header=TRUE)
  header.list[[i]] <- names(file.current)
  header.tf <- i==1
  print(paste("Printing header?",header.tf))
  write.table(file.current,
              file=file.path("merged_chromo", merged.file.name), 
              col.names=header.tf, 
              append=!header.tf,
              row.names=FALSE,
              quote=FALSE)
}

if (length(unique(header.list))>1){
  print("WARNING: SOME HEADERS DO NOT MATCH!")
  print(unique(header.list))
}



