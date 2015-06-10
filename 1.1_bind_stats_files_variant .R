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

#the directory holding the stats files
stats.dir<-"stats/snp_raw"

#read in file names
filenames<-file.path(stats.dir,list.files(stats.dir,pattern="*.gz"))

#headers for the bound file
#column.names<-sapply(strsplit(list.files(stats.dir,pattern="*.txt")),"_"), function(x) x[1:2])

#output file name and directory and with data stamp
dir.create(file.path("stats", "snp_filtered"), showWarnings = FALSE)
date.stamp<-paste("_",format(Sys.time(),"%Y-%m-%d"),sep="")
out.file.name<-file.path("stats/snp_filtered",paste("stats_master_variant",date.stamp,".txt",sep=""))

#filters out invariant sites and reformats for rbinding
filter.fsts<-function(x){
  
  print(paste("Processing",x,"..."))
  matched.file<-data.table(read.table(x,header=TRUE))
  
  ##filter invariant sites
  matched.file<-matched.file[!is.na(matched.file$Fst),]
  matched.file<-matched.file[!matched.file$Fst==Inf,]
  
  #study and comparison names
  study<-sapply(strsplit(x,"/"), function(x) sapply(strsplit(x[length(x)],"_"),function(x)x[1]))
  comparison<-sapply(strsplit(x,"/"), function(x) sapply(strsplit(x[length(x)],"_"),function(x)x[2]))
  
  print("Formatting for output...")
  #data format
  df.tmp<-data.table(lg=matched.file$CHROM,
                     pos=matched.file$POS,
                     n1=matched.file$N1,
                     n2=matched.file$N2,
                     ntotal=matched.file$NTotal,
                     dxy=matched.file$Dxy,
                     fst.num=matched.file$FstNum,
                     fst.den=matched.file$FstDenom,
                     fst=matched.file$Fst,
                     hexp1=matched.file$Hexp1,
                     hexp2=matched.file$Hexp2,
                     freqdif=matched.file$FreqDif,
                     study=study,
                     comparison=comparison)
  
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
lapply(filenames[1:length(filenames)],filter.fsts)
