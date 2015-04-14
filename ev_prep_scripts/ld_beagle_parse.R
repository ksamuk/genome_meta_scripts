######## process a phased VCF file (e.g. BEAGLE output) to a dataframe
######## 1. plot with ggplot

library("ggplot2")
library("dplyr")
library("data.table")

#the phased vcf file
vcf.file <- file("ev_prep_scripts/beagle/IBD_groupIV_beagle.vcf.gz")
vcf <- readLines(vcf.file)

#finds the header so can skip first rows of file
header.line <- grep("#CHROM", vcf)

#vcf to data.frame
vcf.df.raw <- read.table(file=vcf.file,header=TRUE,sep="\t",skip=(header.line-1),comment.char=" ",stringsAsFactors=FALSE)

####filtering and formatting VCF

#1. remove non-biallelic/indels
vcf.df <- vcf.df.raw
vcf.df <- vcf.df[nchar(vcf.df$REF)<=1, ]
vcf.df <- vcf.df[nchar(vcf.df$ALT)<=1, ]

#2. crunch down unused columns
vcf.df <- vcf.df[,c(1,2,4,5,10:length(vcf.df))]

#3. extract phased genotypes and format for ggplot

extractGenotype<- function(x){
  split <- unlist(strsplit(x, ""))
  return(list(split[1], split[3]))
}

df.list<-list()
for (i in 5:length(vcf.df)){
  ind <- vcf.df[,i]
  ind.id <- names(vcf.df)[i]
  geno.ind <- unlist(lapply(as.list(vcf.df[,i]), extractGenotype))
  hap1 <- geno.ind[seq(1, length(geno.ind), 2)]
  hap2 <- geno.ind[seq(2, length(geno.ind), 2)]
  df.len <- length(hap1)
  df.list[[i-4]] <- data.frame(id=rep(ind.id,df.len), 
                            hap=c(rep(1,df.len),rep(2,df.len)), 
                            pos=vcf.df$POS,
                            state=c(hap1,hap2))
}
hap.df<-do.call("rbind",df.list)

########PLOTS

hap.df.beagle<-hap.df

#blue/red palatte
cbPalette<-c("#FF3300","#0000FF")

#blue/orange palatte
cbPalette<-c("#FF9900","#3399FF")

#plot chromosomes with ggplot

#find invariant sites

hap.df.beagle%>%
  filter(pos<100000)%>%
ggplot(.,aes(xmin=pos,xmax=pos+50,ymin=((as.numeric(id)/4)+(hap/10)),ymax=((as.numeric(id)/4)+(hap/10)+0.07),fill=state))+
  geom_rect()+
  geom_text(aes(x=min(pos),y=0.30+as.numeric(id)/4,label=id,hjust=0))+
  scale_fill_manual(values=cbPalette)+
  theme_classic()

tmp<-hap.df%>%
  group_by(pos)%>%
  summarise(num.state=sum(as.numeric(as.character(state))))
