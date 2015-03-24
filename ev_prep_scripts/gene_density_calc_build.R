####Gene density calulator
####Outputs two evs: "gene_density" and "gene_coverage"
####gene_density = number of genes in the window
####    _coverage = proportion of window that is "gene"


library(rtracklayer)
library(biomaRt)
library(dplyr)
library(data.table)

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

#combine duplicate rows in feature dataset

DT <- data.table(gene.list, key="ensembl_gene_id")
DT2<-DT[, list(chromosome_name=chromosome_name,
               start_position=start_position, 
               end_position=end_position, 
               gene_biotype=gene_biotype,
               go_ids=list(go_id)),
        by=ensembl_gene_id]

go.df<-data.frame(ensembl_gene_id=DT$ensembl_gene_id,go_id=DT$go_id)

#collapse duplicated rows
master.list.feat.collapsed<-data.frame(DT2)
master.list.feat.collapsed<-master.list.feat.collapsed[!duplicated(master.list.feat.collapsed$ensembl_gene_id),]