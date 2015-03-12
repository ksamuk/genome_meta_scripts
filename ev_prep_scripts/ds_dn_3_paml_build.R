####create data table from (pre-processed) codeml output
####outputs gene ids and ds values for gacu ONLY
####KS mar-3-2015

########HEAD

library(ape)
library(biomaRt)
library(dplyr)
library(magrittr)

#home dir set up
home.dir<-"E:/Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/review/analysis/gma/ev_prep_scripts/paml_analysis"

#contains pre-processed paml output (the dn and ds tree lines+first file of paml file)
#NB: Files MUST be named exactly as the sequence name in the the PAML file
#e.g. ENSGACP00000014140.cml.txt, and the target branch name in PAML file is "ENSGACP00000014140"
#script will FAIL if any of the files are missing the matching sequence
#I had to manually remove files with missing sequences
paml.output.dir<-file.path(home.dir,"output")

#the output director (window evs)
out.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/evs/window"
out.dir<-"E:/Genome Meta Analysis/evs/window"

########END HEAD

########PROCESS PAML OUTPUT FILES

setwd(paml.output.dir)

#dat file list
file.list<-list.files()

#set up empty output dataframe
gene.id<-vector(mode="character",length=length(file.list))
ds<-vector(mode="numeric",length=length(file.list))
dn<-vector(mode="numeric",length=length(file.list))
gacu.ds<-data.frame(gene.id,ds,dn,stringsAsFactors=FALSE)

for (i in 1:length(file.list)){
  
  #read in file
  file.con<-file(file.list[i])
  file.lines<-readLines(file.con)
  closeAllConnections()
  
  if ((!is.na(file.lines[2]))){  
    
    #make the tree(s)
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
    
    #tracker
    if (i%%1000==0){
      cat(i,"...",sep="")
    }
    
  }
  if (is.na(file.lines[2])){
    gacu.ds$gene.id[i]<-sapply(strsplit(file.list[i],split=".cml.txt"),function(x)x[1])
    gacu.ds$ds[i]<-NA
    gacu.ds$dn[i]<-NA
  }
  
}

gacu.ds<-gacu.ds[complete.cases(gacu.ds),]


########END PROCESS PAML OUTPUT FILES

#####match gene.ids to genomic coordinates

#initialize gacu ensembl
ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")

#the attributes of interest
attributes.feat<-c("ensembl_gene_id",
                  "ensembl_peptide_id",
                   "start_position",
                   "end_position",
                   "chromosome_name")

#query ensembl
coords<-getBM(attributes=attributes.feat,values=gacu.ds$gene.id,filters=c("ensembl_peptide_id"),mart = ensembl)
coords<-coords[!duplicated(coords$ensembl_peptide_id),]
########CLEAN FOR OUTPUT

#the output dataframe
gacu.out<-data.frame(lg=coords$chromosome_name,
                    pos1=coords$start_position,
                    pos2=coords$end_position,
                    ds=gacu.ds$ds,
                    dn=gacu.ds$dn)

#remove scaffolds
gacu.out.2<-gacu.out[grep("scaff*",gacu.out$lg,invert=TRUE),]
#convert lg to numeric (!)
gacu.out.2$lg<-as.numeric(as.roman(sapply(strsplit(as.character(gacu.out.2$lg),"group"),function(x)x[2])))
#arrange
gacu.out.2<-arrange(gacu.out.2,lg)

gacu.out.ds<-gacu.out.2[,1:4]
gacu.out.dn<-gacu.out.2[,c(1:3,5)]

#output
setwd(out.dir)
write.table(gacu.out.ds,file="ds.txt",row.names=FALSE)
write.table(gacu.out.dn,file="dn.txt",row.names=FALSE)

#########END CLEAN FOR OUTPUT
