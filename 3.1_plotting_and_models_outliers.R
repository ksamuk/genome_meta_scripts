#statistical analysis of dvs vs evs

library("dplyr")
library("ggplot2")
library("Hmisc")
library("nlme")
library("lme4")
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/genome_meta_scripts"
#setwd(home.dir)

outlier.folder<-"analysis_ready"
outlier.file<-file.path(outlier.folder,"outliers_analysis_Mar-24-2015.txt")

#read in outlier data
outlier.dat<-read.table(file=outlier.file,header=TRUE,na.strings=c("NA","<NA>"))

####FILTERING EVS

#1. Recombination distances >25cM
outlier.dat$recomb_rate[outlier.dat$recomb_rate>=25]<-NA

#2. Gene true/false
outlier.dat$in.a.gene<-as.numeric(!is.na(outlier.dat$gene_id))

#3. KS
outlier.dat$ks[outlier.dat$ks>=1]<-NA

#4. Outlier
#strips whitespace (...why is there whitespace??) and converts to logical
outlier.dat$outlier<-as.logical(gsub("[[:space:]]", "", as.character(outlier.dat$outlier)))

#5. gene_count
outlier.dat$gene_count[outlier.dat$gene_count>=20]<-NA


#### END FILTERING EVS

####VISUALIZING EVS

###EV NAMES
# [1] "dn"          "ds"          "exon"        "gene_id"     "phastcons"   "pi_atl_10k"  "pi_atl_75k"  "pi_pac_10k" 
# [9] "pi_pac_75k"  "recomb_rate" "utr3"        "utr5"        "in.a.gene"  

#ks vs ds (no relationship)
ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point(alpha=1,size=1)

ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point(alpha=1,size=2)+geom_smooth(method="lm")

#recombination
ggplot(data=outlier.dat,aes(x=log(ds+1),y=log(phastcons+1)))+geom_point()#+facet_wrap(~lg)

#genes
ggplot(data=outlier.dat,aes(x=pos1,y=gene_density))+geom_smooth()+facet_wrap(~lg)

####END VISUALIZING EVS


####VISUALIZING EVS VS. FST
ggplot(data=outlier.dat,aes(x=pos1,as,y=as.numeric(outlier),color=ecotype))+stat_smooth(n=20)+facet_wrap(~lg)

ggplot(data=all.data.out,aes(x=pos,y=fst,color=in.a.gene))+geom_point()+facet_wrap(~lg)

ggplot(data=all.data.out,aes(y=pi_pac_marine,x=as.factor(in.a.gene)))+geom_boxplot()

ggplot(data=all.data.out,aes(x=recomb_rates))+geom_histogram()+facet_wrap(~in.a.gene)

ggplot(data=all.data.out,aes(x=as.factor(in.a.gene),y=recomb_rates))+geom_boxplot()

summary(lm(log(pi_pac_marine+1)~in.a.gene,data=all.data.out))


####THE CATEGORICAL APPROACH
#for each snp, calculate parallelism within and between all ecotypes
#collapse this into a not parallel, parallel, super parallel measure
#NB: MAKE SURE OUTLIER IS STRIPPED OF WHITESPACE AND IS LOGICAL (in filter section above)

#loop through each ecotype and summarize num.obs and num.outlier (recapitulates original file from Diana)
para.out<-list()
for (i in 1:length(levels(outlier.dat$ecotype))){
  para.out[[i]]<-subset(outlier.dat,outlier.dat$ecotype==levels(outlier.dat$ecotype)[i])%>%
    group_by(lg,pos1)%>%
    summarise(num.obs=sum(!is.na(outlier)),num.outlier=sum(as.numeric(outlier),na.rm=TRUE))%>%
    ungroup()
  para.out[[i]]<-data.frame(para.out[[i]])
  prop<-as.numeric(para.out[[i]][,4])/as.numeric(para.out[[i]][,3])
  prop[is.nan(prop)]<-NA
  para.out[[i]]$prop<-prop
  
  names(para.out[[i]])[3]<-paste0(levels(outlier.dat$ecotype)[i],"_",names(para.out[[i]])[3])
  names(para.out[[i]])[4]<-paste0(levels(outlier.dat$ecotype)[i],"_",names(para.out[[i]])[4])
  names(para.out[[i]])[5]<-paste0(levels(outlier.dat$ecotype)[i],"_",names(para.out[[i]])[5])
}



#squish together ecotype counts
#not scaled, do later
para.df<-left_join(para.out[[1]],para.out[[2]])
para.df<-left_join(para.df,para.out[[3]])

#join wide data back with match evs
outlier.big<-data.frame(inner_join(para.df,outlier.dat))
#rm dupes
outlier.big<-outlier.big[duplicated(outlier.big[,1:2]),]

#summarizing across ecotypes
outlier.big$num.obs.all<-rowSums(outlier.big[,grep(".+num.obs",names(outlier.big))])
outlier.big$num.outlier.all<-rowSums(outlier.big[,grep(".+num.outlier",names(outlier.big))])
outlier.big$prop.outlier.all<-outlier.big$num.outlier.all/outlier.big$num.obs.all
outlier.big$prop.outlier.scale<-rowMeans(outlier.big[,grep(".+prop",names(outlier.big))]) #makes any NA = NA
outlier.big$outlier.any<-outlier.big$prop.outlier.scale!=0 #as above, missing in any ecotype = NA
outlier.big$outlier.super<-outlier.big$prop.outlier.scale>=.33 #as above, missing in any ecotype = NA
outlier.big$outlier.ultra<-outlier.big$prop.outlier.scale>=.60 #as above, missing in any ecotype = NA

#distribution of outliers across the genome
ggplot(data=outlier.big,aes(x=pos1,y=gene_count))+geom_point()+facet_wrap(~lg)

#ks vs. scaled
ggplot(data=outlier.big,aes(x=ks,y=prop.outlier.scale))+geom_point()#+facet_wrap(~lg)

ggplot(data=outlier.big,aes(x=ds))+geom_histogram()#+facet_wrap(~lg)

#recomb rate in outlier vs. not outlier
ggplot(data=outlier.big,aes(y=recomb_rate,x=outlier.any))+geom_boxplot()#+facet_wrap(~lg)

#only outlier.any==TRUE (nz = no zeroes)
outlier.big.nz<-subset(outlier.big,outlier.big$outlier.any==TRUE)
outlier.big.nz$super<-outlier.big.nz$prop.outlier.scale>=.33 #as above, missing in any ecotype = NA
outlier.big.nz$ultra<-outlier.big.nz$prop.outlier.scale>=.60 #as above, missing in any ecotype = NA

outlier.big.nz$para.fact<-cut2(outlier.big.nz$prop.outlier.scale, g=5)      # quantile groups

outlier.big.nz%>%<-subset(outlier.big,outlier.big$outlier.any==TRUE)

ggplot(data=outlier.big.nz,aes(y=prop.outlier.scale,x=recomb_rate))+
  geom_point(size=1)+
  stat_smooth(n=20)

#cut into 5 equal sized groups?
outlier.big.nz$para.fact<-cut2(outlier.big.nz$prop.outlier.scale, g=20) 

#recomb
ggplot(data=outlier.big.nz,aes(y=recomb_rate,x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

#pi
ggplot(data=outlier.big.nz,aes(y=pi_pac_75k,x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

#ds
ggplot(data=outlier.big.nz,aes(y=ds,x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

#ks 9sp
ggplot(data=outlier.big.nz,aes(y=ks,x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

#ks 9sp
ggplot(data=outlier.big.nz,aes(y=as.numeric(in.a.gene),x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))



anova(lm(prop.para~recomb_rate,data=outlier.big))
  
####LINEAR MODELS
mod1<-lme(fst~pi_pac_marine+recomb_rates,random=~lg|study,data=all.data,na.action="na.omit")

###EV NAMES
# [1] "dn"          "ds"          "exon"        "gene_id"     "phastcons"   "pi_atl_10k"  "pi_atl_75k"  "pi_pac_10k" 
# [9] "pi_pac_75k"  "recomb_rate" "utr3"        "utr5"        "in.a.gene"  

mod2<-glm(outlier~comparison+recomb_rate+pi_pac_10k+phastcons+in.a.gene,data=outlier.dat,na.action="na.omit",family="binomial")

mod1<-lme(log(fst+1)~in.a.gene,random=~1|pop,data=all.data.out,na.action="na.omit")

ggplot(data=all.data,aes(x=pos,y=fst,color=pop))+geom_smooth()+facet_wrap(~lg)
ggplot(data=all.data.filt,aes(y=log(pi_pac_marine),x=outlier))+geom_boxplot()






