####prep pi data for ev processing
####converts pi files from GO to KS ev format
####KS Mar 07-2015

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/raw"
out.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/window"

#file names of raw pi files from GO
pi.file.names<-c("marines.GATK.atl.pi.10000.txt",
                 "marines.GATK.pac.pi.10000.txt",
                 "marines.GATK.total.filtered.atl.pi.75000.sw.txt",
                 "marines.GATK.total.filtered.pac.pi.75000.sw.txt")

#file names of output files
#ORDER MUST MATCH ABOVE OR YOU FAIL
pi.file.names.out<-c("pi_atl_10k.txt",
                     "pi_pac_10k.txt",
                     "pi_atl_75k.txt",
                     "pi_pac_75k.txt")

#process file files in ev files
for (i in 1:length(pi.file.names)){
  
  setwd(home.dir)
  pi.file<-read.table(file=pi.file.names[i],header=TRUE,fill=TRUE)
  
  lg<-as.numeric(as.roman(sapply(strsplit(as.character(pi.file$Chr),"group"),function(x)x[2])))
  pos1<-pi.file$StartPos
  pos2<-pi.file$EndPos
  pi<-pi.file$Hexp
  pi.out<-data.frame(lg,pos1,pos2,pi)
  
  setwd(out.dir)
  write.table(pi.out,file=pi.file.names.out[i],row.names=FALSE)
  
}