#statistical analysis of dvs vs evs

library("dplyr")
library("ggplot2")
library("Hmisc")
library("nlme")
library("lme4")
library("car")
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/genome_meta_scripts"
#setwd(home.dir)

outlier.folder<-"analysis_ready"
outlier.file<-file.path(outlier.folder,"outliers_analysis_Mar-25-2015.txt")

#read in outlier data
outlier.dat<-read.table(file=outlier.file,header=TRUE,na.strings=c("NA","<NA>"))

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
outlier.dat$gene_count[outlier.dat$gene_count>=20]<-NA

#5. pi
outlier.dat$pi_pac_10k[outlier.dat$pi_pac_10k>=0.015]<-NA

#### END FILTERING

####REFORMATTING FOR OUTLIER COUNTS (i.e. WIDE FORMAT)
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

####END REFORMATTING FOR OUTLIER COUNTS (i.e. WIDE FORMAT)

#########################LINEAR MODELS###############################

###EV NAMES
# [1] "dn"          "ds"          "exon"        "gene_id"     "phastcons"   "pi_atl_10k"  "pi_atl_75k"  "pi_pac_10k" 
# [9] "pi_pac_75k"  "recomb_rate" "utr3"        "utr5"        "in.a.gene"  

###NEW DATAFRAME FOR ELIMINATING PSEUDOREPLICATION (due to long -> wide reformat above, each window shows up many times)

outlier.big.red<-outlier.big%>%
  select(-comparison,-ecotype,-outlier,-pos2)%>%
  distinct()

outlier.big.red$para.fact<-as.factor(cut(outlier.big.red$prop.outlier.scale, breaks=c(seq(0,1,by=0.1))))

#para.cat: number of ecotypes with a paralell outlier at that locus in at least one comparison
outlier.big.red$para.cat<-as.numeric(outlier.big.red[,4]>0)+as.numeric(outlier.big.red[,7]>0)+as.numeric(outlier.big.red[,10]>0)

#prop.outlier.all: proportion of comparisons (any ecotype) with an outlier

#prop.outlier.scale: average proportion of parallel outliers within ecotypes (1=outlier in all comparisons, 0= in none)

#para.

####prop.outlier.scale model
mod1<-outlier.big.red%>%
  filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  with(.,glm(prop.outlier.scale~gene_density*gene_count*pi_pac_10k*recomb_rate*ds,na.action="na.omit",family=quasibinomial))

#type 1 ANOVA
anova(mod1)

#type 2 ANOVA
Anova(mod1,type=2)

####outlier.any model
mod2<-outlier.big.red%>%
  filter(lg!=19)%>%
  with(.,glm(outlier.any~gene_density*gene_count*pi_pac_10k*recomb_rate*ds,na.action="na.omit",family=quasibinomial))

#type 1 ANOVA
anova(mod2)

#type 2 ANOVA
Anova(mod2,type=2)
aov(mod2)



##################PARALLELISM PLOTS######################

#####SCALE PROP PARA VS ALL
outlier.big.red%>%
  #filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
    ggplot(data=.)+
    stat_smooth(aes(x=pos1,y=prop.outlier.scale*20,color="1_parallelism"),n=30,size=1,level=0)+
    stat_smooth(aes(x=pos1,y=recomb_rate,color="2_recombination"),n=30,size=1,level=0)+
    stat_smooth(aes(x=pos1,y=pi_pac_10k*750,color="3_pi"),n=30,size=1,level=0)+
    #stat_smooth(aes(x=pos1,y=phastcons*30,color="phastcons"),n=30,size=1,level=0)+
    stat_smooth(aes(x=pos1,y=gene_count,color="4_gene_count"),n=30,size=1,level=0)+
    stat_smooth(aes(x=pos1,y=ds*10,color="5_ds"),n=30,size=1,level=0)+
    theme(text=element_text(size=16))+
    xlab("Position (bp)")+
    ylab("Genomic variable (various units)")+
    facet_wrap(~lg)
#####

#####PROB OUTLIER
outlier.big.red%>%
  #filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.)+
  stat_smooth(aes(x=pos1,y=as.numeric(outlier.any)*10,color="1_prob_outlier"),n=30,size=1,level=0)+
  stat_smooth(aes(x=pos1,y=recomb_rate,color="2_recombination"),n=30,size=1,level=0)+
  stat_smooth(aes(x=pos1,y=pi_pac_10k*750,color="3_pi"),n=30,size=1,level=0)+
  #stat_smooth(aes(x=pos1,y=phastcons*30,color="phastcons"),n=30,size=1,level=0)+
  stat_smooth(aes(x=pos1,y=gene_count,color="4_gene_count"),n=30,size=1,level=0)+
  stat_smooth(aes(x=pos1,y=ds*10,color="5_ds"),n=30,size=1,level=0)+
  theme(text=element_text(size=16))+
  xlab("Position (bp)")+
  ylab("Genomic variable (various units)")+
  facet_wrap(~lg)
#####

#####PROB OUTLIER
outlier.big.red%>%
  #filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.)+
  stat_smooth(aes(x=pos1,y=as.numeric(outlier.any),color="1_prob_outlier"),n=30,size=1,level=0)+
  #geom_point(aes(x=pos1,y=as.numeric(outlier.any)))+
  theme(text=element_text(size=16))+
  xlab("Position (bp)")+
  ylab("Prob. Outlier")+
  ylim(c(0,1))+
  facet_wrap(~lg)
#####

outlier.big.red%>%
  #filter(prop.outlier.scale>0)%>%
  #filter(lg!=19)%>%
  ggplot(data=.)+
  geom_point(aes(x=recomb_rate,y=gene_count,color="gene_count"))+
  geom_point(aes(x=pos1,y=recomb_rate,color="recombination"))
  #stat_smooth(aes(x=pos1,y=pi_pac_10k*750,color="pi_pac_10k"),n=30,size=1,level=0)+
  #stat_smooth(aes(x=pos1,y=phastcons*30,color="phastcons"),n=30,size=1,level=0)+
  #stat_smooth(aes(x=pos1,y=gene_count,color="gene_count"),n=30,size=1,level=0)+
  #stat_smooth(aes(x=pos1,y=ds*10,color="ds"),n=30,size=1,level=0)+

  facet_wrap(~lg)

#########################VISUALIZING PARALLELISM (WORK IN PROGRESS)###############################

####PROP OUTLIER VS. RECOMB
outlier.big.red%>%
  filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.,aes(x=recomb_rate,y=prop.outlier.scale))+
    geom_point(size=3,alpha=0.5)+
    #stat_smooth(method ="lm",formula=y~log(x),size=1)+
    stat_smooth(method ="loess",n=5)+
    #geom_smooth(formula=y~poly(x,3))+
    theme(text=element_text(size=16))+
    xlab("Recombination rate (cM/Mb)")+
    ylab("Scaled prob. outlier")

####PROP OUTLIER VS. PI
outlier.big.red%>%
  filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.,aes(x=pi_pac_10k,y=prop.outlier.scale))+
  geom_point(size=3,alpha=0.5)+
  geom_smooth(method="lm",formula=y~poly(x,3))+
  #stat_smooth(n=8,size=1)+
  theme(text=element_text(size=16))+
  xlab("Marine nucleotide diversity (pi)")+
  ylab("Scaled prob. outlier")

####PROP OUTLIER VS. PI
outlier.big.red%>%
  #filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.,aes(x=recomb_rate,y=gene_count))+
  geom_hex()+
  #geom_smooth(method="lm")+
  #stat_smooth(n=8,size=1)+
  theme(text=element_text(size=16))
  #xlab("Marine nucleotide diversity (pi)")+
  #ylab("Scaled prob. outlier")

ggplot(data=outlier.big.red,aes(y=log(recomb_rate),x=as.factor(para.cat)))+
  geom_boxplot()+
  theme(text = element_text(size=16))+
  xlab("Number of ecotypes with outlier")+
  ylab("Recombination rate (cM/Mb)")+
  facet_wrap(~lg)

ggplot(data=outlier.big.red,aes(y=log(recomb_rate),x=prop.outlier.all))+
  geom_point()+
  geom_smooth(method="lm")+
  theme(text = element_text(size=16))+
  xlab("Number of ecotypes with outlier")+
  ylab("Recombination rate (cM/Mb)")+
  facet_wrap(~lg)

####model
mod1<-outlier.big.red%>%
  filter(lg!=19)%>%
  with(.,lm(log(recomb_rate)~as.factor(para.cat),na.action="na.omit"))

summary(mod1)

#type 1 ANOVA
anova(mod1)

ggplot(data=outlier.big.red,aes(y=pi_pac_10k,x=as.factor(para.cat)))+
  geom_boxplot()+
  theme(text = element_text(size=16))+
  xlab("Number of ecotypes with outlier")+
  ylab("Marine nucleotide diversity (pi)")

ggplot(data=outlier.big.red,aes(y=ds,x=as.factor(para.cat)))+
  geom_boxplot()+
  theme(text = element_text(size=16))+
  xlab("Number of ecotypes with outlier")+
  ylab("Mutation rate (dS)")

ggplot(data=outlier.big.red,aes(y=gene_count,x=as.factor(para.cat)))+
  geom_boxplot()+
  theme(text = element_text(size=16))+
  xlab("Number of ecotypes with outlier")+
  ylab("Gene density (number of genes)")

ggplot(data=outlier.big.red,aes(y=gene_density,x=as.factor(para.cat)))+
  geom_boxplot()+
  theme(text = element_text(size=16))+
  xlab("Number of ecotypes with outlier")+
  ylab("Gene coverage (% window genic)")

####PROP OUTLIER VS. PI
outlier.big.red%>%
  filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.,aes(x=recomb_rate,y=gene_count))+
  geom_hex()+
  #geom_smooth(method="lm")+
  #stat_smooth(n=8,size=1)+
  theme(text=element_text(size=16))
#xlab("Marine nucleotide diversity (pi)")+
#ylab("Scaled prob. outlier")



outlier.big.red%>%
  filter(prop.outlier.scale>0)%>%
  filter(lg!=19)%>%
  ggplot(data=.,aes(x=recomb_rate,y=gene_count))+
  geom_hex()+
  #geom_smooth(method="lm")+
  #stat_smooth(n=8,size=1)+
  theme(text=element_text(size=16))

#########################PROP ON PROP PLOT(WORK IN PROGRESS)###############################



#########################VISUALIZING EVS (WORK IN PROGRESS)###############################

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

#distribution of outliers across the genome
ggplot(data=outlier.big,aes(x=pos1,y=gene_count))+geom_point()+facet_wrap(~lg)

#ks vs. scaled
outlier.big%>%
  filter(prop.outlier.scale>0)%>%
  ggplot(data=.)+
  geom_point(aes(x=pos1,y=prop.outlier.scale*20,color="red"))+
  geom_point(aes(x=pos1,y=recomb_rate))+
  facet_wrap(~lg)

outlier.big%>%
  filter(prop.outlier.scale>0)%>%
  ggplot(data=.)+
  stat_smooth(aes(x=pos1,y=prop.outlier.scale*20,color="a_parallelism"),n=20)+
  stat_smooth(aes(x=pos1,y=recomb_rate,color="recombination_rate"),n=20)+
  stat_smooth(aes(x=pos1,y=pi_pac_10k*500,color="pi_pac_10k"),n=20)+
  facet_wrap(~lg)

outlier.big%>%
  #filter(prop.outlier.scale>0)%>%
  ggplot(data=.,aes(x=recomb_rate))+geom_histogram()#+facet_wrap(~lg)

#recomb rate in outlier vs. not outlier
ggplot(data=outlier.big,aes(y=recomb_rate,x=outlier.any))+geom_boxplot()#+facet_wrap(~lg)

#cut into equal sized groups?
outlier.big$para.fact<-as.factor(cut(outlier.big$prop.outlier.scale, breaks=c(seq(0,1,by=0.1))))

#recomb
ggplot(data=outlier.big,aes(y=recomb_rate,x=para.fact))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

#pi
ggplot(data=outlier.big,aes(y=pi_pac_75k,x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

########OUTLIER BAR PLOTS
outlier.big.red%>%
  filter(!is.na(para.fact))%>%
  ggplot(data=.,aes(x=para.fact))+
  geom_bar()+
  theme(text=element_text(size=16))+
  xlab("# comparisons sharing outlier window")+
  ylab("# of outlier windows")

outlier.big.red%>%
  filter(!is.na(para.cat))%>%
  ggplot(data=.,aes(x=as.factor(para.cat)))+
  geom_bar()+
  theme(text=element_text(size=16))+
  xlab("# ecotypes sharing outlier window")+
  ylab("# of outlier windows")


#ks 9sp
ggplot(data=outlier.big.nz,aes(y=ks,x=para.fact))+
  geom_violin()+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))

#ks 9sp
ggplot(data=outlier.big,aes(y=ds,x=prop.outlier.scale))+
  geom_point()+
  geom_smooth(method="lm")+
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.25))



anova(lm(prop.para~recomb_rate,data=outlier.big))

