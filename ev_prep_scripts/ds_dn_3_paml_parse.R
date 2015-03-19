#parse paml file to data frame
#new approach (uses raw paml files to avoid shell script shenanigans)

home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/review/analysis/gma/ev_prep_scripts/paml_analysis"

#contains pre-processed paml output (the dn and ds tree lines+first file of paml file)
paml.output.dir<-file.path(home.dir,"alignments_all")

setwd(paml.output.dir)

#dat file list
file.list<-list.files()

gene.id<-vector(mode="character",length=length(file.list))
ds<-vector(mode="numeric",length=length(file.list))
dn<-vector(mode="numeric",length=length(file.list))
gacu.dnds<-data.frame(gene.id,ds,dn,stringsAsFactors=FALSE)

for (i in 1:length(file.list)){
  
  #read in file
  file.con<-file(file.list[i])
  file.lines<-readLines(file.con)
  closeAllConnections()
  
  #find the "ds tree" line
  ds.tree.line<-grep("dS tree",file.lines)
  
  #find gene name from the file name
  gene.name<-sapply(strsplit(file.list[i],split=".cml"),function(x)x[1])
  gacu.dnds$gene.id[i]<-gene.name
  #if no ds tree or target gene, skip file
  #uncomment for commentary on the quality of your data files
  
  if(length(ds.tree.line)==0){
    print(paste(file.list[i],"is missing dS tree."))
    gacu.dnds$ds[i]<-NA
    gacu.dnds$dn[i]<-NA
  }else if(length(grep(paste(gene.name,":",sep=""),file.lines))==0){
    print(paste(file.list[i],"has a dS tree, but is missing the target gene."))
    gacu.dnds$ds[i]<-NA
    gacu.dnds$dn[i]<-NA
  }else{
    
    #make the tree(s)
    ds.tree<-read.tree(text=file.lines[ds.tree.line+1])
    dn.tree<-read.tree(text=file.lines[ds.tree.line+3])
    
    
    
    #if there is no dn or ds value, assign NA, otherwise grab value
    if(is.null(ds.tree$edge.length[which.edge(ds.tree,gene.name)])){
      gacu.dnds$ds[i]<-NA
    }else{
      gacu.dnds$ds[i]<-ds.tree$edge.length[which.edge(ds.tree,gene.name)]
    }
    
    if(is.null(dn.tree$edge.length[which.edge(dn.tree,gene.name)])){
      gacu.dnds$dn[i]<-NA
    }else{
      gacu.dnds$dn[i]<-dn.tree$edge.length[which.edge(dn.tree,gene.name)]
    }
    
    #progress bar
    if (i%%1000==0){
      cat(i,"...",sep="")
    }
  }
}

#remove lines containing NAs
gacu.dnds<-gacu.dnds[complete.cases(gacu.dnds),]

