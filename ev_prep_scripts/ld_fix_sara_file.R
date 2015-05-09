# fix LD files from sara

ev.dir<-file.path(getwd(),"evs")
ld.files<-list.files(ev.dir,"*LD*", full.names=TRUE)

ev.dir<-file.path(getwd(),"evs")

ld.file<-read.csv(ld.files[1])
ld.file$pos1<-ld.file$pos1-1
ld.file$pos2<-ld.file$pos2-1
ld.file <- with(ld.file,data.frame(lg=lg,pos1=pos1,pos2=pos2,ld_r2_atl=ld_r2))
write.table(ld.file,file=file.path(ev.dir,"ld_atl.txt"))

ld.file<-read.csv(ld.files[2])
ld.file$pos2<-ld.file$pos1+10000
ld.file$pos1<-ld.file$pos1-1
ld.file$pos2<-ld.file$pos2-1
ld.file$pos2<-ld.file$pos2-1
ld.file <- with(ld.file,data.frame(lg=lg,pos1=pos1,pos2=pos2,ld_r2_pac=ld_r2))
write.table(ld.file,file=file.path(ev.dir,"ld_pac.txt"))

