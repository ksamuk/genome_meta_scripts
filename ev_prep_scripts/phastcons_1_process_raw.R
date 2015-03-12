##ucsc phast cons scores conversion
##downloaded phastcons data via ftp and unzipped in "conservation" directory
##.data files are actually wiggle (.wig) files

#washy washy
rm(list=ls())

#home dir
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/Conservation"
setwd(home.dir)

#libraries
library(rtracklayer)
library(biomaRt)

#list of wig files in the home dir
files.wig<-list.files()[grep(".data",list.files())]
#trimmed names for writing to file
wig.names<-(unlist(strsplit(files.wig,".data")))

for (i in 1:length(files.wig)){
  wig.tmp<-import.wig(con=files.wig[i])
  df <- data.frame(lg=seqnames(wig.tmp),
                   starts=start(wig.tmp)-1,
                   ends=end(wig.tmp),
                   phastcons=score(wig.tmp)
  )
  write.csv(df,file=paste(wig.names[i],".csv",sep=""),row.names=FALSE)
}


df <- data.frame(lg=seqnames(wig.test),
                 starts=start(wig.test)-1,
                 ends=end(wig.test),
                 phastcons=score(wig.test)
                 )