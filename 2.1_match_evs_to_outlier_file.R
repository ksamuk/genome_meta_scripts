####Match EVs to master stats file
####e.g. all_stats_all_pops type files
####Used for non-parallelism based analysis
####POWERED UP VERSION: IRanges instead of lapply -- about 100x faster!

#washy washy
rm(list=ls())

library("IRanges")
library("dplyr")

#home dir setup
home.dir<-"E:/Genome Meta Analysis"
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/"
setwd(home.dir)

#ev dir setup
ev.dir<-file.path(home.dir,"evs/window")
ev.files<-list.files(ev.dir)

#the "stats" file, specifcally an "outlier" file
stats.file<-read.csv(file="outliers_mar-1-2015.csv",header=TRUE)
#stats.file<-arrange(stats.file,lg,pos1)

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
  ev<-read.table(file=ev.files[i],header=TRUE)
  
  #prep empty ev dataframe
  matched.evs<-data.frame()
  length(matched.evs)<-4
  names(matched.evs)<-c("lg","pos1","pos2","ev")
  
  #sanitize file of NAs and duped rows
  ev<-ev[complete.cases(ev),]
  ev<-arrange(ev,lg,pos1)
  ev<-unique(ev)
  
  ###loop through lgs, matching evs as we god
  for (j in 1:max(stats.file$lg)){
  #for (j in 1:2){
    #subset ev and stat by lg
    ev.chr<-subset(ev,ev$lg==j)
    stat.chr<-subset(stats.file,stats.file$lg==j)
    stat.chr<-stat.chr[,1:3]
    
    #build IRanges objects for overlap finding
    ev.range<-IRanges(start=ev.chr$pos1,end=ev.chr$pos2)
    #minimal set. note that this means overlapping evs end up returning the same value, rather than say, their average
    #(this shouldn't matter for most types of evs)
    ev.range<-reduce(ev.range)
    stat.range<-IRanges(start=stat.chr$pos1,end=stat.chr$pos2)
    
    #find ovelaps amd build an "overlap df"
    overlap<-findOverlaps(stat.range,ev.range,select="all")
    overlap.df<-data.frame(lg=stat.chr[queryHits(overlap),]$lg,
                           pos1=stat.chr[queryHits(overlap),]$pos1,
                           pos2=stat.chr[queryHits(overlap),]$pos2,
                           ev.start=start(ev.range[subjectHits(overlap)]),
                           ev.end=end(ev.range[subjectHits(overlap)]),
                           width=width(ev.range[subjectHits(overlap)]))
    
    overlap.df$ev<-ev.chr[subjectHits(overlap),4]
    overlap.df<-unique(overlap.df)
    
    #calculate the proportion of total overlap
    prop.overlap<-overlap.df%>%
      group_by(pos1)%>%
      summarise(total.overlap=sum(width))
    
    #add it back into the overlap dataframe
    prop.overlap<-data.frame(prop.overlap)
    overlap.df$total.overlap<-prop.overlap$total.overlap[match(overlap.df$pos1,prop.overlap$pos1)]
    
    #the weighted contribution this to the total (deals with multiple ev hits to a single window)
    overlap.df$ev.prop<-overlap.df$ev*(overlap.df$width/overlap.df$total.overlap)
    
    #calculate the weighted evs for each window
    weighted.evs<-overlap.df%>%
      group_by(pos1)%>%
      summarise(ev=sum(ev.prop))
    weighted.evs<-data.frame(weighted.evs)
    #add back into df
    overlap.df$ev.avg<- weighted.evs$ev[match(overlap.df$pos1, weighted.evs$pos1)]
    
    #build output file
    matched.evs.chr<-data.frame(lg=overlap.df$lg,pos1=overlap.df$pos1,pos2=overlap.df$pos2,ev=overlap.df$ev.avg)
    matched.evs.chr<-unique(matched.evs.chr)
    
    #rbind to master file
    matched.evs<-rbind(matched.evs,matched.evs.chr)
    print(paste("matching",ev.files[i],"to stats file","lg",j,"complete"))
  }
  ###end lg loop
  
  #attach real name of ev and cbind to stats file
  names(matched.evs)[4]<-sapply(strsplit(ev.files[i],split=".txt"),function(x)x[1])
  matched.all<-cbind(matched.all,matched.evs[,4])
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




