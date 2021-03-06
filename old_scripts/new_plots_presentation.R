############################################################################

# statistical analysis of dvs vs evs

rm(list=ls())

library("dplyr")
library("reshape2")
library("ggplot2")
library("nlme")
library("lme4")
library("car")
library("visreg")

# data
all.data <- read.table(file=file.path("analysis_ready","stats_75k_2015-06-14.gz"),header=TRUE,na.strings=c("NA","<NA>"))

# set negative FSTs to NA
# this is sketchy, but for now i'm doing it
all.data[!is.na(all.data$fst) & all.data$fst<0, ]$fst <- NA
hist(all.data$fst)
hist(all.data$dxy)

#### magic janky outlier detection with dplyr
is.outlier<-function(x){
  return(x>quantile(x,na.rm=TRUE,probs=0.95)[1])
}

#percentile rank of individual scores (not included)
perc.rank <- function(x) {
  trunc(rank(x)/length(x))
}

all.data.filt<-all.data%>%
  mutate(study_com=paste0(study,comparison)) %>%
  filter(!is.na(fst)) %>%
  #filter(comparison == "paxton")%>%
  group_by(study_com)%>%
  mutate(fst.outlier = is.outlier(fst))%>%
  mutate(dxy.outlier = is.outlier(dxy))%>%
  mutate(both.outlier = dxy.outlier == TRUE & fst.outlier == TRUE)

outlier.dat<-data.frame(ungroup(all.data.filt))

#inject japan as no gene flow

outlier.dat$geography[outlier.dat$study == "japan"] <- "parapatric.d"
outlier.dat$study_com[outlier.dat$study_com == "japangroup1"] <- "parapatric.djapangroup1"

#### end magic janky outlier detection with dplyr

#### FILTERING EVS

# Recombination distances >25cM
outlier.dat$recomb_rate[outlier.dat$recomb_rate>=25]<-NA

# Gene true/false
outlier.dat$in.a.gene<-as.numeric(!is.na(outlier.dat$gene_id))

# KS
outlier.dat$ks[outlier.dat$ks>=1]<-NA

# dS
outlier.dat$ds[outlier.dat$ds>=1]<-NA

# dN
outlier.dat$dn[outlier.dat$dn>=0.4]<-NA

# gene_count
outlier.dat$gene_count[outlier.dat$gene_count>=20]<-NA

# marine pi
outlier.dat$pi_pac_10k[outlier.dat$pi_pac_10k>=0.02]<-NA

############################################################################

# q & d ecotype generator

make.ecotype <- function (x){
  
  study.com.new <- "artificial"
  
  if (grepl("benlim",x)){
    study.com.new <- paste0("1_benlim_",gsub("benlim","",x))
  }
  
  if (grepl("catchen",x)|grepl("fer",x)|grepl("russia",x)|grepl("hohenlohe",x)){
    study.com.new <- paste0("2_marfsh_",x)
  }
  
  if (grepl("roesti",x)){
    study.com.new <- paste0("3_lkstrm_",x)
  }
  
  if (grepl("whtcmn",x)){
    study.com.new <- paste0("4_whtcmn_",x)
  }
  
  if (grepl("japan",x)){
    study.com.new <- paste0("5_japmar_",x)
  }
  
  return(study.com.new)

}

outlier.dat<- outlier.dat %>% 
  rowwise%>%
  mutate(ecotype_study = make.ecotype(study_com))

#### Plots of genomic variables

# RECOMBINATION
outlier.dat %>% 
  filter(lg == 2) %>%
  #ggplot(aes(x = pos1, y = recomb_rate))+
  ggplot(aes(x = pos1, y = recomb_rate))+
  geom_point(size=2)+
  geom_smooth(size=1.5)+
  #geom_line()+
  theme_classic()+
  theme(panel.border = element_rect(fill = 0, colour=1),
        text=element_text(size=18))
  #facet_wrap(~lg)


# MUTATION RATE: ds
outlier.dat %>% 
  ggplot(aes(x = pos1, y = ds))+
  geom_point(size=2)+
  #geom_line()+
  geom_smooth(size=1.5)+
  theme_classic()+
  theme(panel.border = element_rect(fill = 0, colour=1),
        text=element_text(size=18))+
  facet_wrap(~lg)

# GENE DENSITY
outlier.dat %>% 
  ggplot(aes(x = pos1, y = gene_count))+
  geom_point(size=2)+
  #geom_line()+
  geom_smooth(size=1.5)+
  #scale_colour_gradient(limits=c(0, 16), low="blue", high="red", space="Lab")+
  theme_classic()+
  theme(panel.border = element_rect(fill = 0, colour=1),
        text=element_text(size=18))+
  facet_wrap(~lg)

# MARINE PI
outlier.dat %>% 
  ggplot(aes(x = pos1, y = pi_pac_10k))+
  geom_point(size=2)+
  #geom_line()+
  geom_smooth(size=1.5)+
  #scale_colour_gradient(limits=c(0, 0.015), low="blue", high="red", space="Lab")+
  theme_classic()+
  theme(panel.border = element_rect(fill = 0, colour=1),
        text=element_text(size=18))+
  facet_wrap(~lg)

  
#### PLOTS OF DIVERGENCE

# FST HISTOGRAM
outlier.dat %>%
  ggplot(aes(y = fst, x = study_com, fill = geography))+
  geom_boxplot()+
  theme_classic()+
  facet_wrap(~geography, scale = "free_x")+
  theme(text = element_text(size=28))

outlier.dat %>%
  ggplot(aes(x = fst, fill = study_com))+
  geom_histogram()+
  theme_classic()+
  theme(legend.position = "none")+
  facet_wrap(~geography, scale = "free")+
  theme(text = element_text(size=10))

outlier.dat %>%
  with(unique(study_com))
  group_by(geography,study_com) %>%
  tally()
  
  
  ggplot(aes(x = pos1, y = as.numeric(fst.outlier)))+
  geom_point(size = 3)+
  geom_smooth(size = 2, se = FALSE)+
  theme_classic()+
  theme(text = element_text(size=28))+
  
  facet_wrap(~lg)

# FST OUTLIER GENOME

    facet_grid(ecotype_study~lg,space="free_x")

# DXY OUTLIER  GENOME
outlier.dat %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
    ggplot(aes(x = pos1, y = as.numeric(dxy.outlier), color=ecotype_study))+
    geom_smooth(size = 1, se=FALSE)+
    theme_classic()+
    facet_grid(ecotype_study~lg,space="free_x")
    facet_wrap(~lg)

# DXY GENOME
outlier.dat %>%
    filter(!grepl("allo", study_com)) %>%
    filter(!grepl("para", study_com)) %>%
      ggplot(aes(x = pos1, y = dxy, color=ecotype_study))+
      geom_smooth(size = 1, se=FALSE)+
      theme_classic()+
      facet_grid(ecotype_study~lg,space="free_x")

  
# FST GENOME
  outlier.dat %>%
    filter(!grepl("allo", study_com)) %>%
    filter(!grepl("para", study_com)) %>%
    ggplot(aes(x = pos1, y = fst, color=ecotype_study))+
    geom_smooth(size = 1, se=FALSE)+
    theme_classic()+
    facet_grid(ecotype_study~lg,space="free_x")

  
# FST BOX PLOTS  
outlier.dat %>%
    filter(!is.na(fst)) %>%
    filter(!is.na(dxy)) %>%
    filter(!grepl("allo", study_com)) %>%
    filter(!grepl("para", study_com)) %>%
    mutate(ecotype = substr(ecotype_study,1,1) ) %>%
    ggplot(aes(x = reorder(study_com,fst,FUN=function(x)median(x)*-1), y = fst, fill= ecotype))+
    geom_boxplot()+
    #theme_classic()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    facet_grid(~ecotype,space = "free",scale="free_x")

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  ggplot(aes(x = fst, fill= ecotype))+
  geom_histogram()
  #theme_classic()+

# DXY BOX PLOTS  

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  ggplot(aes(x = reorder(study_com,dxy,FUN=function(x)median(x)*-1), y = log(dxy), fill= ecotype))+
  geom_boxplot()+
  #theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(~ecotype,space = "free",scale="free_x")

### PLOTS OF DIVERGENCE VS EVS

##### RECOMBINATION

# REcomb vs. FST per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = fst))+
  geom_point()+
  geom_smooth(method="loess",size=1.5)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))

# REcomb vs. FST outlier per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier)))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,5))

# REcomb vs. DXY per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,dxy,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = dxy))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))

# REcomb vs. DXY.outlier per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,dxy,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = as.numeric(dxy.outlier)))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  #mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = fst.outlier, y = recomb_rate, fill=study.reorder))+
  geom_boxplot()

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  #mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,dxy,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = dxy.outlier, y = recomb_rate, fill=study.reorder))+
  geom_boxplot()



##### RECOMBINATION

##### MARINE PI

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = pi_pac_10k, y = fst))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)

# MARINE PI vs. FST outlier per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = pi_pac_10k, y = as.numeric(fst.outlier)))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)

# MARINE PI vs. DXY per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,dxy,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = pi_pac_10k, y = dxy))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))

# MARINE PI vs. DXY.outlier per comparisons

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,dxy,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = pi_pac_10k, y = as.numeric(dxy.outlier)))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))

#### GENE COUNT / KS : NO PATTERN

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = gene_count, y = fst))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  filter(!grepl("allo", study_com)) %>%
  filter(!grepl("para", study_com)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = ds, y = fst))+
  geom_point()+
  geom_smooth(method="loess",size=2)+
  facet_wrap(~study.reorder)


#### GEOGRAPHY PLOTS

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  #filter(geography == "allopatric.d") %>%
  #filter(geography == "allopatric.s") %>%
  #filter(geography == "parapatric.d") %>%
  #filter(geography == "allopatric.d") %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier)))+
  #geom_point(aes(color = geography))+
  geom_smooth(method="loess",size=2, se=FALSE)+
  facet_grid(~study_com)+
  coord_cartesian(xlim=c(0,10),ylim=c(0,0.20))

outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = fst))+
  geom_point(aes(color = geography))+
  geom_smooth(method="loess",size=2, se=FALSE)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))



outlier.dat %>%
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  mutate(ecotype = substr(ecotype_study,1,1) ) %>%
  mutate(study.reorder = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier)))+
  geom_point(aes(color = geography))+
  geom_smooth(method="loess",size=2, se=FALSE)+
  facet_wrap(~study.reorder)+
  coord_cartesian(xlim=c(0,10))

outlier.dat %>% 
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  #filter(geography == "allopatric.d") %>%
  #filter(geography == "allopatric.s") %>%
  filter(geography == "parapatric.d") %>%
  #filter(geography == "allopatric.d") %>%
  mutate(study_com = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = fst.outlier, y = log(recomb_rate), fill = study_com))+
  #geom_point()+
  geom_boxplot()
  #facet_grid(~geography)

outlier.dat %>% 
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  #filter(geography == "allopatric.d") %>%
  #filter(geography == "allopatric.s") %>%
  #filter(geography == "parapatric.d") %>%
  #filter(geography == "allopatric.d") %>%
  #mutate(study_com = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = fst.outlier, y = log(recomb_rate), fill = study_com))+
  #geom_point()+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position = "none")+
  facet_grid(~geography)
  

outlier.dat %>% 
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  #filter(geography == "allopatric.d") %>%
  #filter(geography == "allopatric.s") %>%
  #filter(geography == "parapatric.d") %>%
  #filter(geography == "allopatric.d") %>%
  #mutate(study_com = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = dxy.outlier, y = log(recomb_rate), fill = geography))+
  #geom_point()+
  theme_bw()+
  geom_boxplot()

  #facet_grid(~geography)

outlier.dat %>% 
  filter(!is.na(fst)) %>%
  filter(!is.na(dxy)) %>%
  mutate(study_com = reorder(study_com,fst,FUN=function(x)mean(x)*-1)) %>%
  ggplot(aes(x = fst.outlier, y = log(recomb_rate), fill = study_com))+
  #geom_point()+
  geom_boxplot()+
  facet_grid(~geography)

outlier.dat %>% 
  ggplot(aes(x = pos1, y = dxy, color = geography))+
  geom_smooth(size=2, se=FALSE)+
  facet_wrap(~lg)

outlier.dat %>% 
  filter(!is.na(fst)) %>%
  ggplot(aes(x = fst.outlier, y = as.numeric(fst.outlier), color = geography))+
  geom_smooth(se=FALSE)+
  facet_wrap(~lg)


outlier.dat %>% 
  ggplot(aes(x = pos1, y = as.numeric(dxy.outlier), color = geography))+
  geom_smooth(size=2, se=FALSE)+
  facet_wrap(~lg)

outlier.dat %>% 
  ggplot(aes(x = pos1, y = as.numeric(fst.outlier), color = geography))+
  geom_smooth(size=2, se=FALSE)+
  facet_wrap(~lg)

mod1<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glm(dxy.outlier~recomb_rate*
               gene.flow*
               div.selection,
             na.action="na.omit",
             family=binomial))


summary(mod1)
Anova(mod1)
anova(mod1)

(-1)/((outlier.dat$hs)^2) %>% hist

mod.full <- mod2
mod.full<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(fst.outlier~log(hs)+
                recomb_rate*
               gene.flow*
               div.selection+
               (1|study_com),
             na.action="na.omit",
             family=binomial))

mod.gene_flow<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(fst.outlier~hs+
                 recomb_rate*
                 gene.flow+
                 (1|study_com),
               na.action="na.omit",
               family=binomial))

mod.selection<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(fst.outlier~hs+
                 recomb_rate*
                 div.selection+
                 (1|study_com),
               na.action="na.omit",
               family=binomial))

mod.recomb<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(fst.outlier~hs+
                 recomb_rate+
                 (1|study_com),
               na.action="na.omit",
               family=binomial))

anova(mod.recomb,mod.selection)
anova(mod.recomb,mod.gene_flow)
anova(mod.recomb,mod.selection, mod.full)

######### DXY

mod.full<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(dxy.outlier~recomb_rate*
                 gene.flow*
                 div.selection+
                 (1|study_com),
               na.action="na.omit",
               family=binomial))


mod.selection<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(dxy.outlier~recomb_rate*
                 div.selection+
                 (1|study_com),
               na.action="na.omit",
               family=binomial))

mod.recomb<-outlier.dat%>%
  filter(lg!=19)%>%
  mutate(gene.flow=as.factor(!grepl("allopatric",study)))%>%
  mutate(div.selection=as.factor(!grepl(".s",study)))%>%
  with(.,glmer(dxy.outlier~recomb_rate+
                 (1|study_com),
               na.action="na.omit",
               family=binomial))

anova(mod.recomb,mod.selection, mod.full)


anova(mod2)

