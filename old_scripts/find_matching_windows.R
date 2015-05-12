####formatting parallel outlier file 
####adds in genomic coordinates
####KS feb 21 2015

library("dplyr")
library("reshape")
library("magrittr")

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/"
setwd(home.dir)

#read in outlier file
#wide format

outliers<-read.csv(file="parallel_outliers_feb-23.csv",header=TRUE,na.strings="NA")

#convert windows back to physical coordinates and join with EV file
#needs the original stats_ev file (all stats all pops evs)

window.size<-75000

#convert the outlier window groups back to genomic coordinates
outliers$lg<-sapply(strsplit(as.character(outliers$lg_bp_group),split="_"),function(x)x[1])
outliers$bp_group<-sapply(strsplit(as.character(outliers$lg_bp_group),split="_"),function(x)x[2])

pos1<-((as.numeric(outliers$bp_group)-1)*window.size)
pos1[pos1<1]<-1
pos2<-((as.numeric(outliers$bp_group))*window.size)-1

outliers$pos1<-pos1
outliers$pos2<-pos2

#write the output file

comparison.columns<-grep("fst",names(outliers))

####ORIGINAL WIDE FORMAT
outliers.out<-data.frame(lg=outliers$lg,
                         pos1=outliers$pos1,
                         pos2=outliers$pos2,
                         outliers[,comparison.columns],
                         outliers[,(max(comparison.columns)+1):(length(outliers)-4)])
####

####LONG FORMAT

#trim wide file of prop calculations
comparison.columns<-grep("fst",names(outliers.out))
outliers.long<-data.frame(lg=outliers.out$lg,
                          pos1=outliers.out$pos1,
                          pos2=outliers.out$pos2,
                          outliers.out[,comparison.columns])
#convert to long
comparison.columns<-grep("fst",names(outliers.long),value=TRUE)
outliers.long<-reshape(outliers.long,varying=comparison.columns,timevar="comparison",times=comparison.columns,v.name="outlier",direction="long")

#fix row names and remove redundant "ID" variable
row.names(outliers.long)<-NULL
outliers.long<-outliers.long[,-6]

#fix comparison names
comparison.names<-tolower(sapply(strsplit(outliers.long$comparison,"_"),function(x)x[2]))
outliers.long$comparison<-comparison.names

#add ecotype data
add.ecotype<-function(x){
  if (x=="boot"|x=="joes"|x=="constance"|x=="geneva"|x=="misty"|x=="roberts"){
    ecotype<-"lake.stream"
  }
  if (x=="catch1"|x=="catch2"|x=="fer"|x=="hohenlohe"){
    ecotype<-"marine.fresh"
  }
  if (x=="lq"|x=="pri"|x=="pax"){
    ecotype<-"ben.lim"
  }
  return(ecotype)
}

outliers.long$ecotype<-sapply(outliers.long$comparison,add.ecotype)
####

