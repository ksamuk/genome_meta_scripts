####create data table from (pre-processed) codeml output
####outputs gene ids and ds values for gacu ONLY
####KS mar-3-2015

library(ape)

#home dir set up
#home.dir<-"E:/Genome Meta Analysis/ev_prep_scripts/paml_analysis"
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/review/analysis/gma/ev_prep_scripts/paml_analysis"

#contained pre-processed paml output (the dn and ds trees)
paml.output.dir<-file.path(home.dir,"output")
setwd(paml.output.dir)

#dat file list
file.list<-list.files()

#set up empty output dataframe
gene.id<-vector(mode="character",length=length(file.list))
ds<-vector(mode="numeric",length=length(file.list))
dn<-vector(mode="numeric",length=length(file.list))
num.species<-vector(mode="numeric",length=length(file.list))
gacu.ds<-data.frame(gene.id,ds,dn,num.species,stringsAsFactors=FALSE)

for (i in 1:length(file.list)){
  
  #read in file
  file.con<-file(file.list[i])
  file.lines<-readLines(file.con)
  closeAllConnections()
  
  if (!is.na(file.lines[2])){  
    
    #make the tree(s)
    num.species<-substr(gsub(" ","",file.lines[1]),1,1)
    ds.tree<-read.tree(text=file.lines[3])
    dn.tree<-read.tree(text=file.lines[5])
    
    #get gene name and find it in the ds tree
    gene.name<-sapply(strsplit(file.list[i],split=".cml.txt"),function(x)x[1])
    
    #pull out ds value for gacu
    gacu.ds$gene.id[i]<-gene.name
    
    if(is.null(ds.tree$edge.length[which.edge(ds.tree,gene.name)])){
      gacu.ds$ds[i]<-NA
    }else{
      gacu.ds$ds[i]<-ds.tree$edge.length[which.edge(ds.tree,gene.name)]
      
    }
    
    if(is.null(dn.tree$edge.length[which.edge(dn.tree,gene.name)])){
      gacu.ds$dn[i]<-NA
    }else{
      gacu.ds$dn[i]<-dn.tree$edge.length[which.edge(dn.tree,gene.name)]
    }
    
    gacu.ds$num.species[i]<-num.species
    
    #tracker
    if (i%%1000==0){
      cat(i,"...",sep="")
    }
    
  }
  if (is.na(file.lines[2])){
    gacu.ds$gene.id[i]<-sapply(strsplit(file.list[i],split=".cml.txt"),function(x)x[1])
    gacu.ds$ds[i]<-NA
    gacu.ds$dn[i]<-NA
    gacu.ds$num.species<-NA
  }
  
}

gacu.ds<-gacu.ds[complete.cases(gacu.ds),]

#output
setwd(home.dir)
write.table(gacu.ds,file="ds_estimates_gacu.txt")

