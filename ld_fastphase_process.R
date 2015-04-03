#fastphase formatting script

library(ggplot2)

closeAllConnections()
rm(list=ls())

########TRIM A FASTPHASE INPUT FILE FOR EXPERIMENTATION PURPOSES
fast.file<-file(file.path(getwd(),"ev_prep_scripts/fastphase/test.in"))
file.lines<-readLines(fast.file)

trim <- function (x){
  if(nchar(x)>100){
    return(substr(x,1,1000))
  }else{
    return(x)
  }
}

trimed.file<-unlist(lapply(file.lines,trim))
fileConn<-file(file.path(getwd(),"ev_prep_scripts/fastphase/test_trim.in"))
writeLines(trimed.file,fileConn)
close(fileConn)
########END

########PARSE A FASTPHASE HAPGUESS FILE

#read in file
hap.file<-file(file.path(getwd(),"ev_prep_scripts/fastphase/fastphase_hapguess_switch.out"))
hap.lines<-readLines(hap.file)
closeAllConnections()

#assumes id lines are indicated by a '#'
id.lines<-grep("#",hap.lines)
hap.df<-data.frame()
for (i in 1:length(id.lines)){
  id<-gsub("# ","",hap.lines[id.lines[i]])
  hap1<-strsplit(hap.lines[id.lines[i]+1]," ")[[1]]
  hap2<-strsplit(hap.lines[id.lines[i]+2]," ")[[1]]
  pos<-c(1:length(hap1))
  hap.df.tmp<-data.frame(id=rep(id,2),hap=c(rep(1,length(hap1)),rep(2,length(hap1))),pos=c(pos,pos),state=c(hap1,hap2))
  hap.df<-rbind(hap.df,hap.df.tmp)
}

ggplot(data=hap.df,aes(xmin=pos,xmax=pos+1,ymin=0,ymax=0.5,color=state))+geom_rect()+facet_wrap(~id)

