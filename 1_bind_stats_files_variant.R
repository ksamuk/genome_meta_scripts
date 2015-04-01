##filters and binds stats files into master file
##filters out invariant sites, binds all dataframes together
##filters out sites with <50% representation across dataset

##run in matched file folder R CMD BATCH filter_fsts.R 
##outputs to new dir "stats_filtered"

library("dplyr")

#the directory holding the stats files
stats.dir<-"/Users/Kieran/Documents/Science/Projects/Ph.D./Genome Meta Analysis/stats"

filenames<-file.path(stats.dir,list.files(stats.dir,pattern="*.txt"))
column.names<-sapply(strsplit(filenames,"_"), function(x) x[1:2])

dir.create(file.path(stats.dir, "stats_filtered"), showWarnings = FALSE)

#filters out invariant sites and reformats for rbinding
filter.fsts<-function(x){
  matched.file<-read.table(x,header=TRUE)
  var.sites.filt1<-matched.file[!is.na(matched.file$Fst),]
  var.sites.filt2<-var.sites.filt1[!var.sites.filt1$Fst==Inf,]
  study<-sapply(strsplit(x,"_"), function(x) x[1])
  comparison<-sapply(strsplit(x,"_"), function(x) x[2])
  df.tmp<-data.frame(lg=var.sites.filt2$CHROM,
                     pos=var.sites.filt2$POS,
                     n1=var.sites.filt2$N1,
                     n2=var.sites.filt2$N2,
                     ntotal=var.sites.filt2$NTotal,
                     dxy=var.sites.filt2$Dxy,
                     fst.num=var.sites.filt2$FstNum,
                     fst.den=var.sites.filt2$FstDenom,
                     fst=var.sites.filt2$Fst,
                     hexp1=var.sites.filt2$Hexp1,
                     hexp2=var.sites.filt2$Hexp2,
                     freqdif=var.sites.filt2$FreqDif,
                     study=study,
                     comparison=comparison)
  print(x)
  return(df.tmp)
}

#apply above function to filename list
frame.list<-lapply(filenames,filter.fsts)

#bind into single df
df.out <- do.call("rbind", frame.list)

#remove scaffolds
df.out2<-df.out[-grep("scaffold*",df.out$lg),]

#rename lgs as numeric
df.out2$lg<-as.numeric(as.roman(sapply(strsplit(as.character(df.out2$lg),split="group"),function(x)x[2])))

#write to file with date stamp
date.stamp<-paste("_",format(Sys.time(),"%b-%d-%Y"),sep="")
file.name<-paste("stats_master",date.stamp,".txt",sep="")
write.table(df.out2,file=paste(file.path(stats.dir,"stats_filtered"),file.name,sep=""),row.names=FALSE)


  