####create data table from (pre-processed) codeml output
####outputs gene ids and ds values for gacu ONLY
####KS mar-3-2015

library(ape)

#home dir set up
#home.dir<-"E:/Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"

setwd(home.dir)

#contained pre-processed paml output (the dn and ds trees)
paml.output.dir<-""

#dat file list
file.list<-list.files(paml.output.dir)

#set up empty output dataframe
gene.id<-vector(mode="character",length=length(file.list))
ds<-vector(mode="numeric",length=length(file.list))
gacu.ds<-data.frame(gene.id,ds)

for (i in 1:length(file.list)){
  
  file.con<-file(file.list[i])
  file.lines<-readLines(file.con)
  ds.tree<-read.tree(text=file.lines[2])
  
  gene.name<-sapply(strsplit(file.list[i],split=".cml.txt"),function(x)x[1])
  gacu.ds$gene.id[i]<-gene.name
  gacu.ds$ds[i]<-ds.tree$edge.length[which.edge(ds.tree, gene.name)]
  if (i%%1000==0){
    print(i)
  }
}


