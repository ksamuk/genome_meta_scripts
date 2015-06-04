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
all.data<-read.table(file=file.path("analysis_ready","stats_75k_2015-06-03.gz"),header=TRUE,na.strings=c("NA","<NA>"))

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


#extract list of outliers
outlier.list <- outlier.dat %>%
  filter(geography=="parapatric.d") %>%
  filter(fst.outlier==TRUE | dxy.outlier == TRUE ) %>%
  select(lg,pos1,pos2,study,fst.outlier,dxy.outlier)



write.table(outlier.list, "sb_meta_outliers.txt")

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


#### Plots of genomic variables

# RECOMBINATION
outlier.dat %>% 
  ggplot(aes(x = pos1, y = recomb_rate, color=recomb_rate))+
  geom_point(size=3)+
  #geom_line()+
  scale_colour_gradient(limits=c(0, 25), low="blue", high="red", space="Lab")+
  theme_classic()+
  theme(panel.border = element_rect(fill = 0, colour=1))+
  facet_wrap(~lg)

# MUTATION RATE
outlier.dat %>% 
  ggplot(aes(x = pos1, y = ks, color=ks))+
  geom_point(size=3)+
  #geom_line()+
  scale_colour_gradient(limits=c(0, 1), low="blue", high="red", space="Lab")+
  theme_classic()+
  theme(panel.border = element_rect(fill = 0, colour=1))+
  facet_wrap(~lg)
  


outlier.dat %>% 
  ggplot(aes(x = pos1, y = fst, color = geography))+
  geom_smooth(size=2, se=FALSE)+
  facet_wrap(~lg)

outlier.dat %>% 
  ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier), color = geography))+
  #geom_point()+
  geom_smooth(size=2, se=FALSE)
  #facet_wrap(~lg)

outlier.dat %>% 
  ggplot(aes(x = fst.outlier, y = recomb_rate, color = geography))+
  #geom_point()+
  geom_boxplot()
#facet_wrap(~lg)

outlier.dat %>% 
  ggplot(aes(x = pos1, y = dxy, color = geography))+
  geom_smooth(size=2, se=FALSE)+
  facet_wrap(~lg)

outlier.dat %>% 
  filter(!is.na(fst)) %>%
  ggplot(aes(x = pos1, y = as.numeric(fst.outlier), color = study_com))+
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


