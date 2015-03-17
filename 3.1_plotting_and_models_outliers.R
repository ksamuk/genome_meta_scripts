#statistical analysis of dvs vs evs

library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("nlme")
library("lme4")
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis"
setwd(home.dir)

#read in outlier data
outlier.dat<-read.table(file="outliers_analysis_Mar-17-2015.txt",header=TRUE,na.strings=c("NA","<NA>"))

####FILTERING EVS

#1. Recombination distances >25cM
outlier.dat$recomb_rate[outlier.dat$recomb_rate>=25]<-NA

#2. Gene true/false
outlier.dat$in.a.gene<-as.numeric(!is.na(outlier.dat$gene_id))

#3. KS
outlier.dat$ks[outlier.dat$ks>=1]<-NA

#### END FILTERING EVS

####VISUALIZING EVS

###EV NAMES
# [1] "dn"          "ds"          "exon"        "gene_id"     "phastcons"   "pi_atl_10k"  "pi_atl_75k"  "pi_pac_10k" 
# [9] "pi_pac_75k"  "recomb_rate" "utr3"        "utr5"        "in.a.gene"  

#pacific marine pi
ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point(alpha=1,size=2)+ylim(c(0,0.02))+geom_smooth()

ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point(alpha=1,size=2)

#recombination
ggplot(data=outlier.dat,aes(x=log(ds),y=log(phastcons)))+geom_point()#+facet_wrap(~lg)

#genes
ggplot(data=all.data.out,aes(x=pos1,y=in.a.gene))+geom_smooth()+facet_wrap(~lg)

####END VISUALIZING EVS


####VISUALIZING EVS VS. FST
ggplot(data=all.data.out,aes(x=pos,y=fst,color=as.factor(in.a.gene)))+geom_smooth()+facet_wrap(~lg)

ggplot(data=all.data.out,aes(x=pos,y=fst,color=in.a.gene))+geom_point()+facet_wrap(~lg)

ggplot(data=all.data.out,aes(y=pi_pac_marine,x=as.factor(in.a.gene)))+geom_boxplot()

ggplot(data=all.data.out,aes(x=recomb_rates))+geom_histogram()+facet_wrap(~in.a.gene)

ggplot(data=all.data.out,aes(x=as.factor(in.a.gene),y=recomb_rates))+geom_boxplot()

summary(lm(log(pi_pac_marine+1)~in.a.gene,data=all.data.out))


####LINEAR MODELS
mod1<-lme(fst~pi_pac_marine+recomb_rates,random=~lg|study,data=all.data,na.action="na.omit")
mod2<-glmer(outlier~pi_pac_marine+recomb_rates+(1|study),data=all.data.out,na.action="na.omit",family=binomial)

mod1<-lme(log(fst+1)~in.a.gene,random=~1|pop,data=all.data.out,na.action="na.omit")



ggplot(data=all.data,aes(x=pos,y=fst,color=pop))+geom_smooth()+facet_wrap(~lg)
ggplot(data=all.data.filt,aes(y=log(pi_pac_marine),x=outlier))+geom_boxplot()






