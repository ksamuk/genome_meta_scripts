#statistical analysis of dvs vs evs

library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("nlme")
library("lme4")
home.dir<-"E:/Genome Meta Analysis"
setwd(home.dir)

all.data<-read.table(file="all_stats_all_pops_evs.txt",header=TRUE,na.strings=c("NA","<NA>"))

#set negative FSTs to NA
#this is sketchy, but for now i'm doing it
all.data[all.data$fst<0,]$fst<-NA
all.data.sub<-subset(all.data,all.data$recomb_rates<50)
hist(all.data$fst)

####magic janky outlier detection with dplyr
is.outlier<-function(x){
  return(x>quantile(x,na.rm=TRUE,probs=0.95)[1])
}

all.data.filt<-all.data.sub%>%
  group_by(pop)%>%
  mutate(outlier=is.outlier(fst))

all.data.out<-data.frame(ungroup(all.data.filt))

####end magic janky outlier detection with dplyr

####FILTERING EVS

#1. Recombination distances >25cM
all.data.out$recomb_rates[all.data.out$recomb_rates>=25]<-NA

#2. Gene true/false
all.data.out$in.a.gene<-as.numeric(!is.na(all.data.out$gene_id))


#### END FILTERING EVS

####VISUALIZING EVS

#pacific marine pi
ggplot(data=all.data,aes(x=pos,y=pi_pac_marine))+geom_smooth()+facet_wrap(~lg)

#recombination
ggplot(data=all.data.sub,aes(x=pos,y=recomb_rates))+stat_smooth(span=1)+facet_wrap(~lg)

#genes
ggplot(data=all.data.out,aes(x=pos,y=in.a.gene))+geom_smooth()+facet_wrap(~lg)

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






