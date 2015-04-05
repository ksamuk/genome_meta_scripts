######## fastPHASE test script
######## 1. generates sets of haplotypes that vary in: length, physical position, LD, missing data, variant/invariant %, num individuals
######## 2. converts to genotypes and runs through fastPHASE via the console                     
######## 3. graphically assesses performance of fastPHASE (plots real vs. inferred haplotypes)

#### libraries
library ("ggplot2")
library ("dplyr")

#### function definitions

generateHaplotype <- function(individual, size, sites, missing){
  
  num.sites <- round(sites*size)
  num.missing <- round(missing*num.sites)
  num.present.sites <- num.sites-num.missing
  
  #chromatid 1
  loci.1 <- sample(size,num.sites)
  loci.state.1 <- rep(1,round(num.present.sites/2))
  loci.state.2 <- rep(2,round(num.present.sites/2))
  loci.missing <- rep("?",num.missing)
  state.1<-c(loci.state.1,loci.state.2,loci.missing)
  
  #chromatid 2
  loci.2 <- sample(size,num.sites)
  loci.state.1 <- rep(1,round(num.present.sites/2))
  loci.state.2 <- rep(2,round(num.present.sites/2))
  loci.missing <- rep("?",num.missing)
  state.2<-c(loci.state.1,loci.state.2,loci.missing)
  
  #print(c(length(loci.state.1),length(loci.state.2),length(loci.missing)))
  hap.out <- data.frame(id=rep(i,num.sites),chrom=c(rep(1,num.sites),rep(2,num.sites)),pos=c(loci.1,loci.2),state=c(state.1,state.2))
  hap.out <- arrange(hap.out,chrom,pos)
  
  return(hap.out)
}

#### body

## haplotype parameters

# number of samples (individuals)
hap.samples <- 10

# size of the region (base pairs)
hap.size <- 10000000

# percentage of sequenced sites
hap.sites <- 0.0001

# percent missing data
hap.missing <- 0.10

# TBD: missing data variability
# hap.missing.var <- 0.0

# TBD: state '1' frequency
# hap.freq<-0.5

# TBD: spacing of sites
# hap.spacing <- 1/10000

# TBD: variability of spacing
#hap.spacing.var <- 0

# TBD: percent invariant sites
#hap.invariant <- 0.10

##generate haplotypes

hap.list<-list()
length(hap.list)<-hap.samples
for (i in 1:hap.samples){
  hap.list[[i]]<-generateHaplotype(i,hap.size,hap.sites,hap.missing)
}

hap.df <- do.call("rbind",hap.list)

##plot generated haplotypes

#blue/orange palatte
cbPalette<-c("#CCCCCC","#FF9900","#3399FF")

#plot chromosomes with ggplot
ggplot(data=hap.df,aes(xmin=pos,xmax=pos+0.25,ymin=((as.numeric(id)/4)+(chrom/10)),ymax=((as.numeric(id)/4)+(chrom/10)+0.07),color=state))+
  geom_rect()+
  scale_color_manual(values=cbPalette)+
  theme_classic()









