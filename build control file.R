####build paml control files for a list of paml alignments
####there's probably a way to not do it like this...but I can't figure it out

library("dplyr")

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev prep scripts/paml_analysis"
setwd(home.dir)

#empty control file
#made from default with some mods
empty.ctl<-readLines(file("codeml_empty.ctl"))


##the header to be made and added to a empty control file
#seqfile = alignments/filename.best.phy
#treefile = alignments/filename.best.dnd
#outfile = codeml_results/filename.ml

#read in file list as vector
file.list<-scan(file="paml_file_list.txt",what="character")
file.list<-sapply(strsplit(file.list,split=".best.phy"),function(x)x[1])

for (i 1:length(file.list)){
  line1<-paste("seqfile = alignments/",file.list[1],".best.phy",sep="")
  line2<-paste("treefile = alignments/",file.list[1],".best.dnd",sep="")
  line3<-paste("outfile = codeml_results/",file.list[1],".ml",sep="")
  full.ctl<-c(line1,line2,line3,empty.ctl)
  cat(full.ctl,fill=1,file="test.ctl.txt")
}
