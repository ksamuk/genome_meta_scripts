#statistical analysis of dvs vs evs

library("dplyr")
library("reshape2")
library("ggplot2")
library("nlme")
library("lme4")

all.data<-read.table(file=file.path("analysis_ready","stats_75k_Apr-11-2015.txt"),header=TRUE,na.strings=c("NA","<NA>"))

#set negative FSTs to NA
#this is sketchy, but for now i'm doing it
all.data[!is.na(all.data$fst) & all.data$fst<0, ]$fst <- NA
hist(all.data$fst)
hist(all.data$dxy)
####magic janky outlier detection with dplyr
is.outlier<-function(x){
  return(x>quantile(x,na.rm=TRUE,probs=0.95)[1])
}

all.data.filt<-all.data%>%
  group_by(comparison)%>%
  mutate(fst.outlier=is.outlier(fst))%>%
  mutate(dxy.outlier=is.outlier(dxy))

outlier.dat<-data.frame(ungroup(all.data.filt))

####end magic janky outlier detection with dplyr

####FILTERING OUTLIER FILE

#1. Recombination distances >25cM
outlier.dat$recomb_rate[outlier.dat$recomb_rate>=25]<-NA

#2. Gene true/false
outlier.dat$in.a.gene<-as.numeric(!is.na(outlier.dat$gene_id))

#3. KS
outlier.dat$ks[outlier.dat$ks>=1]<-NA

#3. dS
outlier.dat$ds[outlier.dat$ds>=1]<-NA

#4. Outlier
#strips whitespace (...why is there whitespace??) and converts to logical
outlier.dat$outlier<-as.logical(gsub("[[:space:]]", "", as.character(outlier.dat$outlier)))

#5. gene_count
all.data.out$gene_count[outlier.dat$gene_count>=20]<-NA

#### END FILTERING

####VISUALIZING EVS

#pacific marine vs . atl marine pi
ggplot(data=outlier.dat,aes(x=pi_pac_10k,y=pi_atl_10k))+geom_point()+facet_wrap(~lg)

#ds vs. pi?
ggplot(data=outlier.dat,aes(x=pi_pac_10k,y=ds))+geom_smooth()+facet_wrap(~lg)

#ds vs. ks FIXED HAHAHAHA
ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point()+geom_smooth(method="lm")

#ds vs. pos
ggplot(data=outlier.dat,aes(x=pos,y=ds))+geom_point()+geom_point()+facet_wrap(~lg)

ggplot(data=outlier.datt,aes(x=pos,y=ds))+geom_point()+stat_smooth(n=5)+facet_wrap(~lg)

#recombination across genome
ggplot(data=outlier.dat,aes(x=pos1,y=recomb_rate))+geom_point()+facet_wrap(~lg)

#genes
ggplot(data=outlier.dat,aes(x=pos,y=in.a.gene))+geom_poi()+facet_wrap(~lg)

####END VISUALIZING EVS


####VISUALIZING EVS VS. FST


ggplot(data=outlier.dat,aes(x=pos,y=fst,color=as.factor(in.a.gene)))+geom_smooth()+facet_wrap(~lg)

ggplot(data=outlier.dat,aes(x=recomb_rate,y=dxy))+geom_point()+geom_smooth()

ggplot(data=all.data.out,aes(y=ds,x=as.factor(outlier)))+geom_boxplot()

ggplot(data=all.data.out,aes(x=recomb_rates))+geom_histogram()+facet_wrap(~in.a.gene)

ggplot(data=all.data.out,aes(x=as.factor(in.a.gene),y=recomb_rates))+geom_boxplot()

summary(lm(log(pi_pac_marine+1)~in.a.gene,data=all.data.out))


####LINEAR MODELS

mod1<-lm(dxy~recomb_rate,data=outlier.dat,na.action="na.omit")
anova(mod1)
mod2<-glmer(outlier~pi_pac_marine+recomb_rates+(1|study),data=all.data.out,na.action="na.omit",family=binomial)

mod1<-lme(log(fst+1)~in.a.gene,random=~1|pop,data=all.data.out,na.action="na.omit")



ggplot(data=all.data,aes(x=pos,y=fst,color=pop))+geom_smooth()+facet_wrap(~lg)
ggplot(data=all.data.filt,aes(y=log(pi_pac_marine),x=outlier))+geom_boxplot()






