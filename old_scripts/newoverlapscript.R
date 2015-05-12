####Match EVs to master stats file
####e.g. all_stats_all_pops type files
####Used for non-parallelism based analysis
####POWERED UP VERSION: IRanges instead of lapply -- about 100x faster!


library("IRanges")
library("GenomicRanges")
library("dplyr")

home.dir<-"E:/Genome Meta Analysis"
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/"
setwd(home.dir)
ev.dir<-file.path(home.dir,"evs/window")
ev.files<-list.files(ev.dir)
stats.file<-read.csv(file="all_stats_all_pops_feb9-2014.csv",header=TRUE)


##match evs to statsfile
##subsets each by chromosome to prevent mismatches
setwd(ev.dir)
ev.files<-list.files()
matched.all<-data.frame(stats.file)

for (i in 1:length(ev.files)){
  
  #read ev file
  ev<-read.table(file=ev.files[i],header=TRUE)
  
  #prep empty ev dataframe
  matched.evs<-data.frame()
  length(matched.evs)<-3
  names(matched.evs)<-c("lg","pos","ev")
  
  #sanitize file of NAs
  ev<-ev[complete.cases(ev),]
  ev<-arrange(ev,lg,pos1)
  
  for (j in 1:max(stats.file$lg)){
  
    #subset ev and stat by lg
    ev.chr<-subset(ev,ev$lg==j)
    stat.chr<-subset(stats.file,stats.file$lg==j)
    stat.chr<-stat.chr[,1:2]
    
    #build IRanges objects for overlap finding
    ev.range<-IRanges(start=ev.chr$pos1,end=ev.chr$pos2)
    stat.range<-IRanges(start=stat.chr$pos,width=1)
    overlap<-data.frame(as.matrix(findOverlaps(stat.range,ev.range)))
    
    #remove duplicates (sometimes a single query hits two e.g. exons, this doesn't matter for our purposes)
    overlap<-overlap[unique(overlap$queryHits),]
    
    #build list of matched pos (from stat file) and ev values (from ev file)
    matched.evs.chr<-data.frame(pos=stat.chr$pos[overlap$queryHits],ev=ev.chr[,length(ev.chr)][overlap$subjectHits])
    
    #join matches for this chromosome to master file (i.e. get the lg column values from the stat file)
    matched.evs.chr<-left_join(stat.chr,matched.evs.chr,by="pos")
    
    #rbind to master file
    matched.evs<-rbind(matched.evs,matched.evs.chr,row.names=FALSE)
    print(paste("matching",ev.files[i],"to stats file","lg",j,"complete"))
  }
  
  names(matched.evs)[3]<-sapply(strsplit(ev.files[i],split=".txt"),function(x)x[1])
  matched.evs.chr<-left_join(matched.all,matched.evs,by=c("pos","lg"))
}

setwd(home.dir)
write.table(stats.file.evs,file="all_stats_all_pops_evs.txt",row.names=FALSE)




