#statistical analysis of dvs vs evs


library("dplyr")
library("reshape2")
library("ggplot2")
library("nlme")
library("lme4")
home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis"
setwd(home.dir)

#read in outlier data
outlier.dat<-read.table(file="outliers_analysis_Mar-20-2015.txt",header=TRUE,na.strings=c("NA","<NA>"))

####FILTERING EVS

#1. Recombination distances >25cM
outlier.dat$recomb_rate[outlier.dat$recomb_rate>=25]<-NA

#2. Gene true/false
outlier.dat$in.a.gene<-as.numeric(!is.na(outlier.dat$gene_id))

#3. KS
outlier.dat$ks[outlier.dat$ks>=1]<-NA

#outlier to logical
#strips whitespace (...why is there whitespace??) and converts to logical

outlier.dat$outlier<-as.logical(gsub("[[:space:]]", "", as.character(outlier.dat$outlier)))

#### END FILTERING EVS

####VISUALIZING EVS

###EV NAMES
# [1] "dn"          "ds"          "exon"        "gene_id"     "phastcons"   "pi_atl_10k"  "pi_atl_75k"  "pi_pac_10k" 
# [9] "pi_pac_75k"  "recomb_rate" "utr3"        "utr5"        "in.a.gene"  

#pacific marine pi
ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point(alpha=1,size=1)

ggplot(data=outlier.dat,aes(x=ks,y=ds))+geom_point(alpha=1,size=2)+geom_smooth(method="lm")

#recombination
ggplot(data=outlier.dat,aes(x=log(ds),y=log(phastcons)))+geom_point()#+facet_wrap(~lg)

#genes
ggplot(data=all.data.out,aes(x=pos1,y=in.a.gene))+geom_smooth()+facet_wrap(~lg)

####END VISUALIZING EVS


####VISUALIZING EVS VS. FST
ggplot(data=all.data.out,aes(x=pos,as,y=as.numeric(outlier),color=study))+stat_smooth(n=10)+facet_wrap(~lg)

ggplot(data=all.data.out,aes(x=pos,y=fst,color=in.a.gene))+geom_point()+facet_wrap(~lg)

ggplot(data=all.data.out,aes(y=pi_pac_marine,x=as.factor(in.a.gene)))+geom_boxplot()

ggplot(data=all.data.out,aes(x=recomb_rates))+geom_histogram()+facet_wrap(~in.a.gene)

ggplot(data=all.data.out,aes(x=as.factor(in.a.gene),y=recomb_rates))+geom_boxplot()

summary(lm(log(pi_pac_marine+1)~in.a.gene,data=all.data.out))


####THE CATEGORICAL APPROACH
#for each snp, calculate parallelism within and between all ecotypes
#collapse this into a not parallel, parallel, super parallel measure
#do stuff?

para.eco<-outlier.dat%>%
  group_by(ecotype,lg,pos1)%>%
  summarise(para.ecotype=mean(outlier,na.rm=TRUE),count=n())%>%
  ungroup()

ggplot(data=para.eco,aes(x=pos1,y=prop.para,color=ecotype))+stat_smooth(n=10)+facet_wrap(~lg)

#using a loop, wow such shame

para.out<-list()
for (i in 1:length(levels(outlier.dat$ecotype))){
  para.out[[i]]<-subset(outlier.dat,outlier.dat$ecotype==levels(outlier.dat$ecotype)[i])%>%
    group_by(lg,pos1)%>%
    summarise(num.obs=sum(!is.na(outlier)),num.outlier=sum(as.numeric(outlier),na.rm=TRUE))%>%
    ungroup()
  
  names(para.out[[i]])[3]<-paste0(levels(outlier.dat$ecotype)[i],"_",names(para.out[[i]])[3])
  names(para.out[[i]])[4]<-paste0(levels(outlier.dat$ecotype)[i],"_",names(para.out[[i]])[4])
}

para.df<-left_join(para.out[[1]],para.out[[2]])
para.df<-left_join(para.df,para.out[[3]])

#join wide data back with match evs
outlier.big<-data.frame(inner_join(para.df,outlier.dat))
#rm dupes
outlier.big<-outlier.big[duplicated(outlier.big[,1:2]),]

outlier.big$prop.para
  prop.1<-(outlier.big[,4]/outlier.big[,3])
  prop.1<-!is.nan(prop.1)*prop.1
  prop.2<-(outlier.big[,6]/outlier.big[,5])
  prop.3<-(outlier.big[,8]/outlier.big[,7])
outlier.big$prop.para<-outlier.big$prop.para/12

ggplot(data=outlier.big,aes(x=pos1,y=prop.para))+stat_smooth(n=10)+facet_wrap(~lg)

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






