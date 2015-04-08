#fastphase formatting script

library(dplyr)
library(ggplot2)
library(RColorBrewer)

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

#convert read lines into dataframe
hap.df<-data.frame()
for (i in 1:length(id.lines)){
  id<-gsub("# ","",hap.lines[id.lines[i]])
  hap1<-strsplit(hap.lines[id.lines[i]+1]," ")[[1]]
  hap2<-strsplit(hap.lines[id.lines[i]+2]," ")[[1]]
  pos<-c(1:length(hap1))
  hap.df.tmp<-data.frame(id=rep(id,2),hap=c(rep(1,length(hap1)),rep(2,length(hap1))),pos=c(pos,pos),state=c(hap1,hap2))
  hap.df<-rbind(hap.df,hap.df.tmp)
}

#match in physical positions
loci.list<-file("ev_prep_scripts/fastphase/marines_pac.groupI.locilist.txt")
loci.list<-as.numeric(readLines(loci.list))

loci.sub.df<-data.frame(index=c(1:max(hap.df$pos)),pos=loci.list[1:max(hap.df$pos)])

hap.df$pos<-loci.sub.df$pos[match(hap.df$pos,loci.sub.df$index)]

#blue/red palatte
cbPalette<-c("#CCCCCC","#FF3300","#0000FF")

#blue/red palatte (no missing)
cbPalette<-c("#FF3300","#0000FF")

#blue/orange palatte
cbPalette<-c("#CCCCCC","#FF9900","#3399FF")

#plot chromosomes with ggplot

hap.df.fp<-hap.df

hap.df.fp%>%
  filter(pos<10000)%>%
  ggplot(.,aes(xmin=pos,xmax=pos+50,ymin=((as.numeric(id)/4)+(hap/10)),ymax=((as.numeric(id)/4)+(hap/10)+0.07),fill=state))+
  geom_rect()+
  geom_text(aes(x=min(pos),y=0.30+as.numeric(id)/4,label=id,hjust=0))+
  scale_fill_manual(values=cbPalette)+
  theme_classic()

