
-----------------
#read in new file from Kieran w/ invariable sites removed
testall_data2 <- read.table("latest_fst_allpop.txt", header = T)


#create a population variable rather than study and comparison
testall_data2$pop <-paste0(testall_data2$study, testall_data2$comparison)

#add column indicating ecotype 

testall_data2$ecotype <- "streamlake" 
testall_data2$ecotype[testall_data2$pop =="catchengroup1"] <- "marinefresh"
testall_data2$ecotype[testall_data2$pop =="catchengroup2"] <- "marinefresh"
testall_data2$ecotype[testall_data2$pop =="fer"] <- "marinefresh"
testall_data2$ecotype[testall_data2$pop =="hohenlohe"] <- "marinefresh"
testall_data2$ecotype[testall_data2$pop =="benlimlq"] <- "benthiclimnetic"
testall_data2$ecotype[testall_data2$pop =="benlimpaxton"] <- "benthiclimnetic"
testall_data2$ecotype[testall_data2$pop =="benlimpriest"] <- "benthiclimnetic"

#write out this revised file 

write.csv(testall_data2, file="revisedalldata.csv")
testall_data2 <- read.csv("revisedalldata.csv")
-----
#create a vector with population names:

population <- unique(testall_data2$pop)

#add a counting variable 
testall_data2$cnt <- 1


#Identify threshold for top 5% of fst values 
z <- tapply(testall_data2$fst, INDEX=testall_data2$pop, quantile, probs=0.95)

#Generates separate dataframes for each population
for (i in 1:length(population)){i_temp <- testall_data2[testall_data2 $pop ==population[i],]; assign(paste("ddat",i,sep="_"),i_temp)}
#now have thirteen ddat data frames 

win=75000 # set window size. This is what we chose as the optimal 

#still maybe want to write a loop for this 
ddat_1$bp_group <- as.integer(ddat_1$pos/win)+1 #LQ 
ddat_1$lg_bp_group <- sprintf("%s_%s",ddat_1$lg,ddat_1$bp_group)
ddat_2$bp_group <- as.integer(ddat_2$pos/win)+1 #Pax
ddat_2$lg_bp_group <- sprintf("%s_%s",ddat_2$lg,ddat_2$bp_group)
ddat_3$bp_group <- as.integer(ddat_3$pos/win)+1 #Pri 
ddat_3$lg_bp_group <- sprintf("%s_%s",ddat_3$lg,ddat_3$bp_group)
ddat_4$bp_group <- as.integer(ddat_4$pos/win)+1 #Catch1
ddat_4$lg_bp_group <- sprintf("%s_%s",ddat_4$lg,ddat_4$bp_group)
ddat_5$bp_group <- as.integer(ddat_5$pos/win)+1 #Catch2
ddat_5$lg_bp_group <- sprintf("%s_%s",ddat_5$lg,ddat_5$bp_group)
ddat_6$bp_group <- as.integer(ddat_6$pos/win)+1 #Fer
ddat_6$lg_bp_group <- sprintf("%s_%s",ddat_6$lg,ddat_6$bp_group)
ddat_7$bp_group <- as.integer(ddat_7$pos/win)+1 #hohenlo
ddat_7$lg_bp_group <- sprintf("%s_%s",ddat_7$lg,ddat_7$bp_group)
ddat_8$bp_group <- as.integer(ddat_8$pos/win)+1 #boot
ddat_8$lg_bp_group <- sprintf("%s_%s",ddat_8$lg,ddat_8$bp_group)
ddat_9$bp_group <- as.integer(ddat_9$pos/win)+1 #constance
ddat_9$lg_bp_group <- sprintf("%s_%s",ddat_9$lg,ddat_9$bp_group)
ddat_10$bp_group <- as.integer(ddat_10$pos/win)+1 #geneva
ddat_10$lg_bp_group <- sprintf("%s_%s",ddat_10$lg,ddat_10$bp_group)
ddat_11$bp_group <- as.integer(ddat_11$pos/win)+1 #joes
ddat_11$lg_bp_group <- sprintf("%s_%s",ddat_11$lg,ddat_11$bp_group)
ddat_12$bp_group <- as.integer(ddat_12$pos/win)+1 #misty
ddat_12$lg_bp_group <- sprintf("%s_%s",ddat_12$lg,ddat_12$bp_group)
ddat_13$bp_group <- as.integer(ddat_13$pos/win)+1 #roberts
ddat_13$lg_bp_group <- sprintf("%s_%s",ddat_13$lg,ddat_13$bp_group)

library(plyr)
library(parallel)

#note below says prop_3perc but this is actually 5perc variable just not recoded

island1<- ddply(ddat_1,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[1])/sum(cnt))
island1 <- island1[with(island1, order(chrom, bp_group)), ]
island1 <- subset(island1, island1$nsnps>2)

#This bit was written for checking window sizes 
#island2 <- subset(island1, island1$nsnps>2)
#x <- nrow(island1)
#y <- nrow(island2)
#y/x			
			
permutation1 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[1]))/nsnps)
}
phi1 <- as.data.frame(do.call(rbind, mclapply(island1$nsnps, permutation1, data=ddat_1, mc.cores=4)))
island_perm_fst1 <- data.frame(island1,phi1)
save(island_perm_fst1,file="island_N_perm_fst1.R")
#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst1)) {
  perm<-island_perm_fst1[r,7:1006]
  h<-island_perm_fst1[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst1<-data.frame(island1,pval)
save(island_perm_pval_fst1,file="island_N_perm_pval_fst1.R")

#P value detection .05 threshold 
island_perm_pval_fst1$padjust<-p.adjust(island_perm_pval_fst1$pval,method="fdr")
island_perm_pval_fst1$fst_0.05<-island_perm_pval_fst1$padjust<=0.05 
island_perm_pval_fst1$win<-(island_perm_pval_fst1$bp_group*75000)-37500 

#check for number of sig windows 

del=subset(island_perm_pval_fst1,island_perm_pval_fst1$fst_0.05=="TRUE")
nrow(del)
288/3838

#-------------- Further code for additional pops this needs to be edited still 

island2<- ddply(ddat_2,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[2])/sum(cnt))
island2 <- island2[with(island2, order(chrom, bp_group)), ]
island2 <- subset(island2, island2$nsnps>2)
			
			
permutation2 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[2]))/nsnps)
}			
phi2 <- as.data.frame(do.call(rbind, mclapply(island2$nsnps, permutation2, data=ddat_2, mc.cores=4)))
island_perm_fst2 <- data.frame(island2,phi2)
save(island_perm_fst2,file="island_N_perm_fst2.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst2)) {
  perm<-island_perm_fst2[r,7:1006]
  h<-island_perm_fst2[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst2<-data.frame(island2,pval)
save(island_perm_pval_fst2,file="island_N_perm_pval_fst2.R")

#P value detection .05 threshold 
island_perm_pval_fst2$padjust<-p.adjust(island_perm_pval_fst2$pval,method="fdr")
island_perm_pval_fst2$fst_0.05<-island_perm_pval_fst2$padjust<=0.05 
island_perm_pval_fst2$win<-(island_perm_pval_fst2$bp_group*75000)-37500 

del=subset(island_perm_pval_fst2,island_perm_pval_fst2$fst_0.05=="TRUE")
nrow(del)
390/4615 


#----
island3<- ddply(ddat_3,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[3])/sum(cnt))
island3 <- island3[with(island3, order(chrom, bp_group)), ]
island3 <- subset(island3, island3$nsnps>2)

permutation3 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[3]))/nsnps)
}			
phi3 <- as.data.frame(do.call(rbind, mclapply(island3$nsnps, permutation3, data=ddat_3, mc.cores=4)))
island_perm_fst3 <- data.frame(island3,phi3)
save(island_perm_fst3,file="island_N_perm_fst3.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst3)) {
  perm<-island_perm_fst3[r,7:1006]
  h<-island_perm_fst3[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst3<-data.frame(island3,pval)
save(island_perm_pval_fst3,file="island_N_perm_pval_fst3.R")

#P value detection .05 threshold 
island_perm_pval_fst3$padjust<-p.adjust(island_perm_pval_fst3$pval,method="fdr")
island_perm_pval_fst3$fst_0.05<-island_perm_pval_fst3$padjust<=0.05 
island_perm_pval_fst3$win<-(island_perm_pval_fst3$bp_group*75000)-37500 


del=subset(island_perm_pval_fst3,island_perm_pval_fst3$fst_0.05=="TRUE")
nrow(del)

414/4653

#-----------
island4<- ddply(ddat_4,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[4])/sum(cnt))
island4 <- island4[with(island4, order(chrom, bp_group)), ]
island4 <- subset(island4, island4$nsnps>2)
permutation4 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[4]))/nsnps)
}		
phi4 <- as.data.frame(do.call(rbind, mclapply(island4$nsnps, permutation4, data=ddat_4, mc.cores=4)))	
island_perm_fst4 <- data.frame(island4,phi4)
save(island_perm_fst4,file="island_N_perm_fst4.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst4)) {
  perm<-island_perm_fst4[r,7:1006]
  h<-island_perm_fst4[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst4<-data.frame(island4,pval)
save(island_perm_pval_fst4,file="island_N_perm_pval_fst4.R")

#P value detection .05 threshold 
island_perm_pval_fst4$padjust<-p.adjust(island_perm_pval_fst4$pval,method="fdr")
island_perm_pval_fst4$fst_0.05<-island_perm_pval_fst4$padjust<=0.05 
island_perm_pval_fst4$win<-(island_perm_pval_fst4$bp_group*75000)-37500 

del=subset(island_perm_pval_fst4,island_perm_pval_fst4$fst_0.05=="TRUE")
nrow(del)
285/4713

#------
island5<- ddply(ddat_5,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[5])/sum(cnt))
island5 <- island5[with(island5, order(chrom, bp_group)), ]
island5 <- subset(island5, island5$nsnps>2)

permutation5 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[5]))/nsnps)
}	
phi5 <- as.data.frame(do.call(rbind, mclapply(island5$nsnps, permutation5, data=ddat_5, mc.cores=4)))		
island_perm_fst5 <- data.frame(island5,phi5)
save(island_perm_fst5,file="island_N_perm_fst5.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst5)) {
  perm<-island_perm_fst5[r,7:1006]
  h<-island_perm_fst5[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst5<-data.frame(island5,pval)
save(island_perm_pval_fst5,file="island_N_perm_pval_fst5.R")

#P value detection .05 threshold 
island_perm_pval_fst5$padjust<-p.adjust(island_perm_pval_fst5$pval,method="fdr")
island_perm_pval_fst5$fst_0.05<-island_perm_pval_fst5$padjust<=0.05 
island_perm_pval_fst5$win<-(island_perm_pval_fst5$bp_group*75000)-37500 

del=subset(island_perm_pval_fst5,island_perm_pval_fst5$fst_0.05=="TRUE")
nrow(del)
287/4674

#---------
	
island6<- ddply(ddat_6,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[6])/sum(cnt))
island6 <- island6[with(island6, order(chrom, bp_group)), ]
island6 <- subset(island6, island6$nsnps>2)

permutation6 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[6]))/nsnps)
}	
phi6 <- as.data.frame(do.call(rbind, mclapply(island6$nsnps, permutation6, data=ddat_6, mc.cores=4)))	
island_perm_fst6 <- data.frame(island6,phi6)
save(island_perm_fst6,file="island_N_perm_fst6.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst6)) {
  perm<-island_perm_fst6[r,7:1006]
  h<-island_perm_fst6[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst6<-data.frame(island6,pval)
save(island_perm_pval_fst6,file="island_N_perm_pval_fst6.R")


#P value detection .05 threshold 
island_perm_pval_fst6$padjust<-p.adjust(island_perm_pval_fst6$pval,method="fdr")
island_perm_pval_fst6$fst_0.05<-island_perm_pval_fst6$padjust<=0.05 
island_perm_pval_fst6$win<-(island_perm_pval_fst6$bp_group*75000)-37500 

del=subset(island_perm_pval_fst6,island_perm_pval_fst6$fst_0.05=="TRUE")
nrow(del)

98/3838

#-------------
		
island7<- ddply(ddat_7,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[7])/sum(cnt))
island7 <- island7[with(island7, order(chrom, bp_group)), ]
island7 <- subset(island7, island7$nsnps>2)
			
permutation7 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[7]))/nsnps)
}
phi7 <- as.data.frame(do.call(rbind, mclapply(island7$nsnps, permutation7, data=ddat_7, mc.cores=4)))
island_perm_fst7 <- data.frame(island7,phi7)
save(island_perm_fst7,file="island_N_perm_fst7.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst7)) {
  perm<-island_perm_fst7[r,7:1006]
  h<-island_perm_fst7[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst7<-data.frame(island7,pval)
save(island_perm_pval_fst7,file="island_N_perm_pval_fst7.R")

#P value detection .05 threshold 
island_perm_pval_fst7$padjust<-p.adjust(island_perm_pval_fst7$pval,method="fdr")
island_perm_pval_fst7$fst_0.05<-island_perm_pval_fst7$padjust<=0.05 
island_perm_pval_fst7$win<-(island_perm_pval_fst7$bp_group*75000)-37500 

#check for number of sig windows 

del=subset(island_perm_pval_fst7,island_perm_pval_fst7$fst_0.05=="TRUE")
nrow(del)
196/3213

#--------------------
island8<- ddply(ddat_8,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[8])/sum(cnt))
island8 <- island8[with(island8, order(chrom, bp_group)), ]
island8 <- subset(island8, island8$nsnps>2)
permutation8 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[8]))/nsnps)
}
phi8 <- as.data.frame(do.call(rbind, mclapply(island8$nsnps, permutation8, data=ddat_8, mc.cores=4)))		
island_perm_fst8 <- data.frame(island8,phi8)
save(island_perm_fst8,file="island_N_perm_fst8.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst8)) {
  perm<-island_perm_fst8[r,7:1006]
  h<-island_perm_fst8[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst8<-data.frame(island8,pval)
save(island_perm_pval_fst8,file="island_N_perm_pval_fst8.R")

#P value detection .05 threshold 
island_perm_pval_fst8$padjust<-p.adjust(island_perm_pval_fst8$pval,method="fdr")
island_perm_pval_fst8$fst_0.05<-island_perm_pval_fst8$padjust<=0.05 
island_perm_pval_fst8$win<-(island_perm_pval_fst8$bp_group*75000)-37500 

del=subset(island_perm_pval_fst8,island_perm_pval_fst8$fst_0.05=="TRUE")
nrow(del)
143/2867

#-----------	
island9<- ddply(ddat_9,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[9])/sum(cnt))
island9 <- island9[with(island9, order(chrom, bp_group)), ]
island9 <- subset(island9, island9$nsnps>2)
permutation9 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[9]))/nsnps)
}
phi9 <- as.data.frame(do.call(rbind, mclapply(island9$nsnps, permutation9, data=ddat_9, mc.cores=4)))
island_perm_fst9 <- data.frame(island9,phi9)
save(island_perm_fst9,file="island_N_perm_fst9.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst9)) {
  perm<-island_perm_fst9[r,7:1006]
  h<-island_perm_fst9[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst9<-data.frame(island9,pval)
save(island_perm_pval_fst9,file="island_N_perm_pval_fst9.R")

#P value detection .05 threshold 
island_perm_pval_fst9$padjust<-p.adjust(island_perm_pval_fst9$pval,method="fdr")
island_perm_pval_fst9$fst_0.05<-island_perm_pval_fst9$padjust<=0.05 
island_perm_pval_fst9$win<-(island_perm_pval_fst9$bp_group*75000)-37500 

del=subset(island_perm_pval_fst9,island_perm_pval_fst9$fst_0.05=="TRUE")
nrow(del)
87/3374

#-----			
island10<- ddply(ddat_10,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[10])/sum(cnt))
island10 <- island10[with(island10, order(chrom, bp_group)), ]
island10 <- subset(island10, island10$nsnps>2)
		
permutation10 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[10]))/nsnps)
}	
phi10 <- as.data.frame(do.call(rbind, mclapply(island10$nsnps, permutation10, data=ddat_10, mc.cores=4)))	
island_perm_fst10 <- data.frame(island10,phi10)
save(island_perm_fst10,file="island_N_perm_fst10.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst10)) {
  perm<-island_perm_fst10[r,7:1006]
  h<-island_perm_fst10[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst10<-data.frame(island10,pval)
save(island_perm_pval_fst10,file="island_N_perm_pval_fst10.R")

#P value detection .05 threshold 
island_perm_pval_fst10$padjust<-p.adjust(island_perm_pval_fst10$pval,method="fdr")
island_perm_pval_fst10$fst_0.05<-island_perm_pval_fst10$padjust<=0.05 
island_perm_pval_fst10$win<-(island_perm_pval_fst10$bp_group*75000)-37500 

del=subset(island_perm_pval_fst10,island_perm_pval_fst10$fst_0.05=="TRUE")
nrow(del)
66/1954

#-------
island11<- ddply(ddat_11,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[11])/sum(cnt))
island11 <- island11[with(island11, order(chrom, bp_group)), ]
island11 <- subset(island11, island11$nsnps>2)
		
permutation11 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[11]))/nsnps)
}	
phi11 <- as.data.frame(do.call(rbind, mclapply(island11$nsnps, permutation11, data=ddat_11, mc.cores=4)))	
island_perm_fst11 <- data.frame(island11,phi11)
save(island_perm_fst11,file="island_N_perm_fst11.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst11)) {
  perm<-island_perm_fst11[r,7:1006]
  h<-island_perm_fst11[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst11<-data.frame(island11,pval)
save(island_perm_pval_fst11,file="island_N_perm_pval_fst11.R")

#P value detection .05 threshold 
island_perm_pval_fst11$padjust<-p.adjust(island_perm_pval_fst11$pval,method="fdr")
island_perm_pval_fst11$fst_0.05<-island_perm_pval_fst11$padjust<=0.05 
island_perm_pval_fst11$win<-(island_perm_pval_fst11$bp_group*75000)-37500 

del=subset(island_perm_pval_fst11,island_perm_pval_fst11$fst_0.05=="TRUE")
nrow(del)
129/2797

#---
island12<- ddply(ddat_12,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[12])/sum(cnt))
island12 <- island12[with(island12, order(chrom, bp_group)), ]
island12 <- subset(island12, island12$nsnps>2)
		
permutation12 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[12]))/nsnps)
}	
phi12 <- as.data.frame(do.call(rbind, mclapply(island12$nsnps, permutation12, data=ddat_12, mc.cores=4)))	
island_perm_fst12 <- data.frame(island12,phi12)
save(island_perm_fst12,file="island_N_perm_fst12.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst12)) {
  perm<-island_perm_fst12[r,7:1006]
  h<-island_perm_fst12[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst12<-data.frame(island12,pval)
save(island_perm_pval_fst12,file="island_N_perm_pval_fst12.R")


#P value detection .05 threshold 
island_perm_pval_fst12$padjust<-p.adjust(island_perm_pval_fst12$pval,method="fdr")
island_perm_pval_fst12$fst_0.05<-island_perm_pval_fst12$padjust<=0.05 
island_perm_pval_fst12$win<-(island_perm_pval_fst12$bp_group*75000)-37500 

del=subset(island_perm_pval_fst12,island_perm_pval_fst12$fst_0.05=="TRUE")
nrow(del)
80/3095

#---
island13<- ddply(ddat_13,~lg_bp_group,summarise,
			chrom=unique(lg),
			bp_group=unique(bp_group),
			nsnps=length(fst),
			mean=mean(fst),
			prop_3perc=sum(fst>= z[13])/sum(cnt))
island13 <- island13[with(island13, order(chrom, bp_group)), ]
island13 <- subset(island13, island13$nsnps>2)
		
permutation13 <- function(data,nsnps){
  perm <- 1:1000
  sapply(perm,function (x) sum(data[sample(1:nrow(data), nsnps,replace=FALSE),]$fst>= (z[13]))/nsnps)
}	
phi13 <- as.data.frame(do.call(rbind, mclapply(island13$nsnps, permutation13, data=ddat_13, mc.cores=4)))	
island_perm_fst13 <- data.frame(island13,phi13)
save(island_perm_fst13,file="island_N_perm_fst13.R")

#assign pvalues
pval_fun<-function(perm,h){length(perm[perm>h])/1000}
pval<-data.frame()
for (r in 1:nrow(island_perm_fst13)) {
  perm<-island_perm_fst13[r,7:1006]
  h<-island_perm_fst13[r,6]
  pvalh<-pval_fun(perm,h)
  pval<-rbind(pval,pvalh)
}  ###need to double check that col 6 (prop_3perc is correct) input here
colnames(pval)<-c("pval")
island_perm_pval_fst13<-data.frame(island13,pval)
save(island_perm_pval_fst13,file="island_N_perm_pval_fst13.R")

#P value detection .05 threshold 
island_perm_pval_fst13$padjust<-p.adjust(island_perm_pval_fst13$pval,method="fdr")
island_perm_pval_fst13$fst_0.05<-island_perm_pval_fst13$padjust<=0.05 
island_perm_pval_fst13$win<-(island_perm_pval_fst13$bp_group*75000)-37500 

del=subset(island_perm_pval_fst13,island_perm_pval_fst13$fst_0.05=="TRUE")
nrow(del)
118/3474
	
***********	
#we will now merge all of these data frames together using the bp_group colum as a reference
#combined <- merge(island_perm_pval_fst2, island_perm_pval_fst1, by= "lg_bp_group") #this is just written for two pops need to change for more 

combined <- read.csv("combining.csv") #this is just a dateframe containing all the lg_bp_group across all datasets and thus used to tie them together 

combined$fst1_0.05 <- island_perm_pval_fst1$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst1$lg_bp_group)]
combined$fst2_0.05 <- island_perm_pval_fst2$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst2$lg_bp_group)]
combined$fst3_0.05 <- island_perm_pval_fst3$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst3$lg_bp_group)]
combined$fst4_0.05 <- island_perm_pval_fst4$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst4$lg_bp_group)]
combined$fst5_0.05 <- island_perm_pval_fst5$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst5$lg_bp_group)]
combined$fst6_0.05 <- island_perm_pval_fst6$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst6$lg_bp_group)]
combined$fst7_0.05 <- island_perm_pval_fst7$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst7$lg_bp_group)]
combined$fst8_0.05 <- island_perm_pval_fst8$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst8$lg_bp_group)]
combined$fst9_0.05 <- island_perm_pval_fst9$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst9$lg_bp_group)]
combined$fst10_0.05 <- island_perm_pval_fst10$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst10$lg_bp_group)]
combined$fst11_0.05 <- island_perm_pval_fst11$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst11$lg_bp_group)]
combined$fst12_0.05 <- island_perm_pval_fst12$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst12$lg_bp_group)]
combined$fst13_0.05 <- island_perm_pval_fst13$fst_0.05[match(combined $lg_bp_group, island_perm_pval_fst13$lg_bp_group)]

#Total comparison 
combined$overalltrue <- sapply(levels,function(x)rowSums(combined=="TRUE",na.rm=TRUE)) 
combined$overallfalse <- sapply(levels,function(x)rowSums(combined[2:14]=="FALSE",na.rm=TRUE)) 


#B-L comparison 
combined$bltrue <- sapply(levels,function(x)rowSums(combined[2:4]=="TRUE",na.rm=TRUE)) 
combined$blfalse <- sapply(levels,function(x)rowSums(combined[2:4]=="FALSE",na.rm=TRUE)) 

#M-L comparison
combined$mltrue <- sapply(levels,function(x)rowSums(combined[5:8]=="TRUE",na.rm=TRUE)) 
combined$mlfalse <- sapply(levels,function(x)rowSums(combined[5:8]=="FALSE",na.rm=TRUE)) 

#S-L comparison
combined$sltrue <- sapply(levels,function(x)rowSums(combined[9:14]=="TRUE",na.rm=TRUE)) 
combined$slfalse <- sapply(levels,function(x)rowSums(combined[9:14]=="FALSE",na.rm=TRUE)) 


write.csv(combined, file="parallel.csv")


#Testing parallelism
#create a column to deterime number of parallel windows
#combined$parallel <- 0
#combined$parallelprop <- 0
#combined$ecotypeprop <- 0 

2:4 benthic-limnetic
5:8 marine-fresh
9:14 stream-lake 


#count the number of outliers between populations (paralism) need to revise for multiple pops
combined$parallel[combined$fst_0.05.x =="TRUE" | combined$fst_0.05.y =="TRUE"] <- 1 
combined$parallel[combined$fst_0.05.x =="TRUE" & combined$fst_0.05.y =="TRUE"] <- 2

#Significance testing of parallelism ***This needs to be redone for more pops****

#create a column to deterime number of parallel windows
combined$parallel <- 0
#count the number of outliers between populations (paralism)
combined$parallel[combined$fst_0.05.x =="TRUE" | combined$fst_0.05.y =="TRUE"] <- 1 
combined$parallel[combined$fst_0.05.x =="TRUE" & combined$fst_0.05.y =="TRUE"] <- 2

########
# pvalue generating code 

pval_fun<-function(perm,h){length(perm[perm>=h])/perms} #replace this 

perm=apply(permutation_output,2,FUN=function(x) length(which(x=='TRUE')))
h=sum(combined$parallel==2)

pval<-pval_fun(perm,h)

		
#will want to check for parallelism within ecotypes and between 
#we can subset out ecotypes for separate tests then run whole dataframe together. 

##??? continuous measure of parallelism???

test <- colnames(combined[,(2:14)])

scp rennison@zoology.ubc.ca:/Newdisk/SciBorg/array0/rennison/devil.pdf ./


scp rennison@zoology.ubc.ca:/Newdisk/SciBorg/array0/rennison/parallel.csv ./