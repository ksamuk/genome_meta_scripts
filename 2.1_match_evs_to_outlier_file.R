####Match EVs to master stats file
####e.g. all_stats_all_pops type files
####Used for non-parallelism based analysis
####POWERED UP VERSION: IRanges instead of lapply -- about 100x faster!

#washy washy
rm(list=ls())

library("IRanges")
library("GenomicRanges")
library("dplyr")

#home.dir<-"E:/Genome Meta Analysis"
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/"
setwd(home.dir)
ev.dir<-file.path(home.dir,"evs/window")
ev.files<-list.files(ev.dir)
stats.file<-read.csv(file="outliers_mar-1-2015.csv",header=TRUE)
stats.file<-arrange(stats.file,lg,pos)


##match evs to statsfile
##subsets each by chromosome to prevent mismatches

#find the ev files
setwd(ev.dir)
ev.files<-list.files()

#initialize matched data variable
matched.all<-data.frame(stats.file)
start.time<-Sys.time()

###looping through all the ev files
for (i in 1:length(ev.files)){
  
  #read ev file
  ev<-read.table(file=ev.files[1],header=TRUE)
  
  #prep empty ev dataframe
  matched.evs<-data.frame()
  length(matched.evs)<-4
  names(matched.evs)<-c("lg","pos1","pos2","ev")
  
  #sanitize file of NAs and duped rows
  ev<-ev[complete.cases(ev),]
  ev<-arrange(ev,lg,pos1)
  ev<-unique(ev)
  
  #VERY IMPORTANT NOTE: this collapses evs down if any of them overlap (i.e. makes a minimal set of ev ranges)
  #this was added to deal with overlapping alternatively spliced exons 
  #the concensus ranges then all get the same value (1 in this case)
  #this could really mess stuff up if your data aren't formatted correctly
  
  #this was wi
  ev.range.all<-IRanges(start=ev$pos1,end=ev$pos2)
  
  
  ###loop through lgs, matching evs as we god
  for (j in 1:max(stats.file$lg)){
  for (j in 1:1){
    #subset ev and stat by lg
    ev.chr<-subset(ev,ev$lg==j)
    stat.chr<-subset(stats.file,stats.file$lg==j)
    stat.chr<-stat.chr[,1:3]
    
    #build IRanges objects for overlap finding
    ev.range<-IRanges(start=ev.chr$pos1,end=ev.chr$pos2)
    ev.range<-reduce(ev.range)
    stat.range<-IRanges(start=stat.chr$pos1,end=stat.chr$pos2)
    
    #find ovelaps amd build an "overlap df"
    overlap<-findOverlaps(stat.range,ev.range,select="all")
    overlap.df<-data.frame(stat.chr[queryHits(overlap),],
                           ev.start=start(ev.range[subjectHits(overlap)]),
                           ev.end=end(ev.range[subjectHits(overlap)]),
                           width=width(ev.range[subjectHits(overlap)]),
                           width.percent=width(ev.range[subjectHits(overlap)])/75000)
    overlap.df<-unique(overlap.df)
    
    #calculate the proportion of overlap
    prop.overlap<-overlap.df%>%
      group_by(pos1,pos2)%>%
      summarise(sum(width.percent))
    
    prop.overlap<-overlap.df%>%
      group_by(pos1,pos2)%>%
      summarise(sum(width)/75000)
    
    
    #build list of matched pos (from stat file) and ev values (from ev file)
    matched.evs.chr<-data.frame(lg=stat.chr$lg,pos1=stat.chr$pos1,pos2=stat.chr$pos2,ev=ev.chr[,length(ev.chr)][overlap])
        
    #rbind to master file
    matched.evs<-rbind(matched.evs,matched.evs.chr)
    print(paste("matching",ev.files[i],"to stats file","lg",j,"complete"))
  }
  ###end lg loop
  
  #attach real name of ev and cbind to stats file
  names(matched.evs)[3]<-sapply(strsplit(ev.files[i],split=".txt"),function(x)x[1])
  matched.all<-cbind(matched.all,matched.evs[,3])
  names(matched.all)[length(matched.all)]<-names(matched.evs)[3]
  
}
###end ev loop

#output report and write file
end.time<-Sys.time()
print(paste("complete!","matching took",(end.time-start.time),"seconds"))
print("preview of output file:")
print(head(matched.all))
setwd(home.dir)
write.table(matched.all,file="all_stats_all_pops_evs.txt",row.names=FALSE)




