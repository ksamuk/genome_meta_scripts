#####filters and binds stats files into master file
####filters out invariant sites
####outputs to new dir "stats_filtered"
#####KS Feb 2015

#required libraries
require("dplyr")
require("data.table")

#roman to numeric function for lg naming 
chrom.to.num<-function(x){
  chrom.rom<-as.character(as.roman(c(1:21)))
  return(match(x,chrom.rom))
}

#the directory holding the stats files (one dir above "stats" to avoid pushing big data to github)
stats.dir<-paste(gsub("genome_meta_scripts","stats",getwd()))

#read in file names
filenames<-file.path(stats.dir,list.files(stats.dir,pattern="*sliding*"))

#headers for the bound file
#column.names<-sapply(strsplit(list.files(stats.dir,pattern="*.txt")),"_"), function(x) x[1:2])

#output file name and directory and with data stamp
dir.create(file.path(stats.dir, "stats_75k_filtered"), showWarnings = FALSE)
date.stamp<-paste("_",format(Sys.time(),"%b-%d-%Y"),sep="")
out.file.name<-file.path(stats.dir, "stats_75k_filtered",paste("stats_75k_master",date.stamp,".txt",sep=""))

#filters out invariant sites and reformats for rbinding
filter.fsts<-function(x){
  
  print(paste("Processing",x,"..."))
  matched.file<-fread(x,header=TRUE)
  
  ##filter invariant sites
  #matched.file<-matched.file[!is.na(matched.file$Fst),]
  #matched.file<-matched.file[!matched.file$Fst==Inf,]
  
  #study and comparison names
  study<-strsplit(gsub(paste(stats.dir,"/",sep=""),"",x),"_")[[1]][1]
  comparison<-strsplit(gsub(paste(stats.dir,"/",sep=""),"",x),"_")[[1]][2]
  
  print("Formatting for output...")
  print(study)
  print(comparison)
  #data format
  df.tmp<-data.table(lg=matched.file$Chr,
                     pos1=matched.file$StartPos,
                     pos2=matched.file$EndPos,
                     total.sites=matched.file$TotalSites,
                     variable.sites=matched.file$VariableSites,
                     dxy=matched.file$Dxy,
                     fst=matched.file$Fst,
                     hexp1=matched.file$Hexp1,
                     hexp2=matched.file$Hexp2)
  df.tmp$study<-study
  df.tmp$comparison<-comparison
  
  #remove scaffolds
  df.tmp<-df.tmp[grep("group",df.tmp$lg),]
  
  #rename lgs as numeric
  df.tmp$lg<-chrom.to.num(gsub("group","",df.tmp$lg))
  
  print("Writing to master file...")
  #append to file
  if (file.exists(out.file.name)){
    write.table(df.tmp,file=out.file.name,row.names=FALSE,append=TRUE,col.names=FALSE)
  }else{
    print("Master file does not exist, creating new file...")
    write.table(df.tmp,file=out.file.name,row.names=FALSE)
  }

}

#apply above function to filename list
lapply(filenames,filter.fsts)


  