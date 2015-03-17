#fix for recomb rate window estimates
#forces windows to be non-overlapping

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis"
ev.dir<-file.path(home.dir,"evs/window")

setwd(ev.dir)

recomb<-read.table(file="recomb_rate.txt",header=TRUE)

recomb$pos1<-recomb$pos1+1
recomb$pos2<-recomb$pos2-1

recomb<-recomb[complete.cases(recomb),]
head(recomb[,-c(4:5)]

write.table(recomb,file="recomb_rate.txt",row.names=FALSE)
