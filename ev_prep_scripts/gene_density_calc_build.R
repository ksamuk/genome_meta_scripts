####Gene density calulator
####Outputs two evs: "gene_density" and "gene_coverage"
####gene_density = number of genes in the window
####    _coverage = proportion of window that is "gene"


library(rtracklayer)
library(biomaRt)
library(dplyr)
library(data.table)

########BIOMART: GENE LOCATIONS

#set the biomart
gacu.ensembl <- useMart("ensembl",dataset="gaculeatus_gene_ensembl")

#gene attributes from the "feature" page
attributes.feat<-c("start_position",
                    "end_position",
                    "chromosome_name")


gene.list.raw<-getBM(attributes=attributes.feat,values="*",mart = gacu.ensembl)

#remove scaffolds
gene.list<-gene.list.raw[grep("group",gene.list.raw$chromosome_name),]

#chromo to numeric
gene.list$chromosome_name<-as.numeric(as.roman(gsub("group","",gene.list$chromosome_name)))
gene.list<-arrange(gene.list,chromosome_name,start_position)

#standard formatting
gene.list<-data.frame(lg=gene.list$chromosome_name,pos1=gene.list$start_position,pos2=gene.list$end_position)

########END BIOMART: GENE LOCATIONS

########SLIDING WINDOWS OF GENE COUNT

window.size<-75000

#kind of a grody implementation

gene.counts.chr<-list()
for (i in 1:max(gene.list$lg)){
  gene.list.chr<-subset(gene.list,gene.list$lg==i)
  
  #convert nuc position to window index
  window.nums<-as.integer(gene.list.chr$pos1/window.size)+1
  
  #tally up genes in windows
  window.names<-as.numeric(names(table(window.nums)))
  window.counts<-as.numeric(table(window.nums))
  
  #reformatting 
  gene.counts.chr[[i]]<-data.frame(win.num=window.names,gene.count=window.counts)
  
  #calc window boundaries
  window.pos1<-(gene.counts.chr[[i]]$win.num-1)*window.size
  window.pos1[window.pos1==0]<-1
  window.pos2<-((gene.counts.chr[[i]]$win.num)*window.size)-1
  
  #the output for this chromo
  lg<-rep(i,length(window.pos2))
  gene.counts.chr[[i]]<-data.frame(lg=lg,pos1=window.pos1,pos2=window.pos2,gene.count=window.counts)
}

#bind
gene.counts.out<-do.call("rbind",gene.counts.chr)

#write to file
write.table(gene.counts.out,file=file.path(getwd(),"evs","gene_count.txt"),row.names=FALSE)

########END SLIDING WINDOWS OF GENE COUNT
