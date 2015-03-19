#statistical analysis of dvs vs evs

library("dplyr")
library("reshape2")
library("ggplot2")
library("nlme")
library("lme4")

#home.dir<-"E:/Genome Meta Analysis"
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/"
setwd(home.dir)

all.data<-read.table(file="snp_analysis_Mar-19-2015.txt",header=TRUE,na.strings=c("NA","<NA>"))

#set negative FSTs to NA
#this is sketchy, but for now i'm doing it
all.data[all.data$fst<0,]$fst<-NA
hist(all.data$fst)

####magic janky outlier detection with dplyr
is.outlier<-function(x){
  return(x>quantile(x,na.rm=TRUE,probs=0.95)[1])
}

all.data.filt<-all.data%>%
  group_by(pop)%>%
  mutate(outlier=is.outlier(fst))

all.data.out<-data.frame(ungroup(all.data.filt))

####end magic janky outlier detection with dplyr

####FILTERING EVS

#1. Recombination distances >25cM
all.data.out$recomb_rate[all.data.out$recomb_rate>=25]<-NA

#2. Gene true/false
all.data.out$in.a.gene<-as.numeric(!is.na(all.data.out$gene_id))

#3. KS
all.data.out$ks[all.data.out$ks>=2]<-NA

#### END FILTERING EVS

####VISUALIZING EVS

#pacific marine vs . atl marine pi
ggplot(data=all.data.out,aes(x=pi_pac_10k,y=pi_atl_10k))+geom_point()+facet_wrap(~lg)

#ds vs. pi?
ggplot(data=all.data.out,aes(x=pi_pac_10k,y=ds))+geom_smooth()+facet_wrap(~lg)

#ds vs. ks
ggplot(data=all.data.out,aes(x=ks,y=ds))+geom_point()#+facet_wrap(~lg)

#ds vs. pos
ggplot(data=all.data.out,aes(x=pos,y=ds))+geom_point()+geom_smooth()+facet_wrap(~lg)

#recombination across genome
ggplot(data=all.data.out,aes(x=pos,y=recomb_rate))+geom_point()+facet_wrap(~lg)

#genes
ggplot(data=all.data.out,aes(x=pos,y=in.a.gene))+geom_poi()+facet_wrap(~lg)

####END VISUALIZING EVS


####VISUALIZING EVS VS. FST
ggplot(data=all.data.out,aes(x=pos,y=fst,color=as.factor(in.a.gene)))+geom_smooth()+facet_wrap(~lg)

ggplot(data=all.data.out,aes(x=pos,y=fst,color=in.a.gene))+geom_point()+facet_wrap(~lg)

ggplot(data=all.data.out,aes(y=ds,x=as.factor(outlier)))+geom_boxplot()

ggplot(data=all.data.out,aes(x=recomb_rates))+geom_histogram()+facet_wrap(~in.a.gene)

ggplot(data=all.data.out,aes(x=as.factor(in.a.gene),y=recomb_rates))+geom_boxplot()

summary(lm(log(pi_pac_marine+1)~in.a.gene,data=all.data.out))


####LINEAR MODELS
mod1<-lme(fst~in.a.gene+ds+pi_pac_10k+recomb_rate,random=~1|study,data=all.data.out,na.action="na.omit")
anova(mod1)
mod2<-glmer(outlier~pi_pac_marine+recomb_rates+(1|study),data=all.data.out,na.action="na.omit",family=binomial)

mod1<-lme(log(fst+1)~in.a.gene,random=~1|pop,data=all.data.out,na.action="na.omit")



ggplot(data=all.data,aes(x=pos,y=fst,color=pop))+geom_smooth()+facet_wrap(~lg)
ggplot(data=all.data.filt,aes(y=log(pi_pac_marine),x=outlier))+geom_boxplot()






