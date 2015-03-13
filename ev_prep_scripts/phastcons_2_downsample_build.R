####Sliding window for downsampling high-res EV data
####KS 02/2015
####Specific goal is to downsample Gacu PHASTCONs data from UCSC genome
####input file format:- csv: lg,pos1,pos2,phastcons (or any ev)

#great R library or greatest R library?
library(data.table)
library(dplyr)

#WINDOWS
home.dir<-"E:/Genome Meta Analysis/evs/nonwindow/conservation"
out.dir<-"E:/Genome Meta Analysis/evs/window"

#MAC OS
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/nonwindow/conservation"
out.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/window"
setwd(home.dir)

#find all the csv files in the home dir
file.list<-grep(".csv",value=TRUE,list.files())

#define window size(in base pairs), ev name and outdir 

window.size<-75000
ev.name<-"phastcons"

#calc windowed stats for each coordinate run
start<-Sys.time()
ev.chr.window.list<-list()
for (i in 1: length(file.list)){
  
  #set up the windows for the ev file
  ev.chr<-fread(file.list[i])
  window.num<-as.integer(ev.chr$pos1/window.size)+1

  ev.chr<-data.table(lg=ev.chr$lg,
                     pos1=ev.chr$pos1,
                     pos2=ev.chr$pos2,
                     window.num,
                     ev=ev.chr[,get(ev.name)])
  
  #summarize windows by window number
  ev.chr.window<-ev.chr%>%
    group_by(window.num)%>%
    summarise(mean(ev))
  
  #coerce to df
  ev.chr.window<-data.frame(ev.chr.window)
  
  #calc window boundaries
  window.pos1<-(ev.chr.window$window.num-1)*75000
  window.pos1[window.pos1==0]<-1
  window.pos2<-((ev.chr.window$window.num)*75000)-1
  
  #final ev.chr
  ev.chr.window<-data.frame(lg=i,pos1=window.pos1,pos2=window.pos2,ev=ev.chr.window[,2])
  names(ev.chr.window)[4]<-ev.name
  ev.chr.window.list[[i]]<-ev.chr.window
}
print(Sys.time()-start)

ev.out<-do.call("rbind",ev.chr.window.list)

write.table(ev.out,file=file.path(out.dir,paste(ev.name,".txt",sep="")),row.names=FALSE)

