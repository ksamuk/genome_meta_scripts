####Match EVs to master stats file
####e.g. all_stats_all_pops type files
####Used for non-parallelism based analysis
####POWERED UP VERSION: IRanges instead of lapply -- about 100x faster!

#washy washy
rm(list=ls())

library("IRanges")
library("dplyr")
install_github("hadley/readr")


home.dir<-gsub("/genome_meta_scripts","",getwd())

ev.dir<-file.path(getwd(),"evs")

stats.file<-read.csv(file=file.path(home.dir,"all_stats_all_pops_feb9-2014.csv"),header=TRUE)
stats.file<-arrange(stats.file,lg,pos)


##match evs to statsfile
##subsets each by chromosome to prevent mismatches

#find the ev files
ev.files<-list.files(ev.dir)

#initialize matched data variable
matched.all<-data.frame(stats.file)
start.time<-Sys.time()

###looping through all the ev files
for (i in 1:length(ev.files)){
  
  #read ev file
  ev<-read.table(file=file.path(ev.dir,ev.files[i]),header=TRUE)
  
  #prep empty ev dataframe
  matched.evs<-data.frame()
  length(matched.evs)<-3
  names(matched.evs)<-c("lg","pos","ev")
  
  #sanitize file of NAs
  ev<-ev[complete.cases(ev),]
  ev<-arrange(ev,lg,pos1)
  
  ###loop through lgs, matching evs as we god
  for (j in 1:max(stats.file$lg)){

    #subset ev and stat by lg
    ev.chr<-subset(ev,ev$lg==j)
    stat.chr<-subset(stats.file,stats.file$lg==j)
    stat.chr<-stat.chr[,1:2]
    
    #build IRanges objects for overlap finding
    ev.range<-IRanges(start=ev.chr$pos1,end=ev.chr$pos2)
    stat.range<-IRanges(start=stat.chr$pos,width=1)
    
    #find ovelaps, returns list of length query with NAs for no hits
    overlap<-findOverlaps(stat.range,ev.range,select="first")
    
    #build list of matched pos (from stat file) and ev values (from ev file)
    matched.evs.chr<-data.frame(lg=stat.chr$lg,pos=stat.chr$pos,ev=ev.chr[,length(ev.chr)][overlap])
        
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

date.stamp<-paste("_",format(Sys.time(),"%b-%d-%Y"),sep="")
file.name<-paste("snp_outliers_analysis",date.stamp,".gz",sep="")
write.table(matched.all,file=gzfile(file.path(getwd(),"analysis_ready",file.name)),row.names=FALSE)




