#filter ds and dn values

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis"
ev.dir<-file.path(home.dir,"evs/window")

setwd(ev.dir)

ds.dat<-read.table(file="ds.txt",header=TRUE)
dn.dat<-read.table(file="dn.txt",header=TRUE)

#only keep rows in both dataframe where ds >= 2 (probably alignment errors)
ds.out<-ds.dat[ds.dat$ds<=2,]
dn.out<-dn.dat[ds.dat$ds<=2,]

write.table(ds.out,file="ds.txt",row.names=FALSE)
write.table(dn.out,file="dn.txt",row.names=FALSE)

