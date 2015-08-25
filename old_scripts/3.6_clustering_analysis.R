# calculate nearest neighbour distances

library(data.table)
library(dplyr)
library(parallel)
library(ggplot2)

# read in snp file
snp.file <- list.files(file.path("stats/snp_filtered"),pattern=".txt$", full.names = TRUE)
snp.file <- fread(snp.file)

# remove weird outliers
snp.file[!is.na(snp.file$fst) & snp.file$fst<0, ]$fst <- NA

# find outliers

is.outlier<-function(x){
  return(x >= quantile(x,na.rm=TRUE,probs=0.95)[1])
}

snp.file<-snp.file%>%
  group_by(study_com)%>%
  mutate(fst.outlier = is.outlier(fst)) %>% 
  ungroup

############### NEW APPROACH: PERMUTATION FTW

sample_sites <- function (gen.pos, num.outliers) {
  
  site.sample <- gen.pos %>% sample_n(num.outliers) %>% arrange(gen.pos) %>% select(gen.pos) %>% with(as.numeric(gen.pos)) %>% diff %>% mean(na.rm = TRUE)
  return(site.sample)
  
}

permute_distances_snp_list <- function(x, num.samples){
  
  mean.distances.lg.all <- list()  
  
  for (j in unique(x$lg)){
    
    dist.sub.lg <- x %>%
      filter(lg == j) %>%
      filter(!is.na(gen.pos)) %>%
      arrange(gen.pos)
    
    num.outliers <- sum(as.numeric(dist.sub.lg$fst.outlier), na.rm = TRUE)
    
    if (num.outliers > 1){
      
      mean.distances.lg <- replicate(num.samples, sample_sites(dist.sub.lg, num.outliers))
      mean.distances.lg.df <- data.frame(mean.distance = mean.distances.lg)
      mean.distances.lg.df$lg <- j
      mean.distances.lg.all[[j]] <- mean.distances.lg.df
    }
    
  }
  mean.distances.lg.all <- do.call("rbind", mean.distances.lg.all)
  mean.distances.lg.all$study_com = x$study_com[1]
  rownames(mean.distances.lg.all) <- NULL
  return(mean.distances.lg.all)
  
}

split.df <- split(snp.file, snp.file$study_com)  
null.distances <- mclapply(split.df, permute_distances_snp_list, num.samples = 10000, mc.cores = 8, mc.silent = FALSE)

null.distances <- do.call("rbind", null.distances)
rownames(null.distances) <- NULL
write.table(null.distances, file="null_snp_distances.txt", row.names = FALSE)

#### START

null.distances <- fread("null_snp_distances.txt")
#null.distances <- null.distances[,-1,with=FALSE]
#col.names <- c("mean.distance", "lg", "study_com")
#setnames(null.distances, col.names)

null.distances %>%
  filter(study_com == "benlim_lq") %>%
  ggplot(aes(x = mean.distance))+
  geom_histogram()+
  facet_wrap(~lg)

################# calc distances between outliers

snp.outliers <- snp.file %>%
  filter(fst.outlier == TRUE) 

dist.outliers <- data.frame()

for (i in unique(snp.outliers$study_com)){
  print(paste0("processing ", i,"..."))
  snp.outlier.sub <- subset(snp.outliers, snp.outliers$study_com == i)
  diff.sub <- data.frame()
  
  for (j in unique(snp.outlier.sub$lg)){
    snp.outlier.sub.lg <- snp.outlier.sub %>%
      filter(lg == j) %>%
      arrange(gen.pos)
    
    diff.lg <- data.frame(outlier.dist = diff(snp.outlier.sub.lg$gen.pos))
    if (length(diff.lg$outlier.dist)>0){
      diff.lg$lg <- j
      diff.sub <- rbind(diff.sub, diff.lg)
    }
    
  }
  
  diff.sub$study <- snp.outlier.sub$study[1]
  diff.sub$comparison <- snp.outlier.sub$comparison[1]
  diff.sub$study_com <- snp.outlier.sub$study_com[1]
  
  dist.outliers <- rbind(dist.outliers, diff.sub)
  
}

# inject geography

add.geo <- function (x) {
  
  geography<-"parapatric.d"
  
  if(grepl("allopatricD",x)==TRUE){
    geography<-"allopatric.d"
  }
  
  if(grepl("allopatricS",x)==TRUE){
    geography<-"allopatric.s"
  }
  
  if(grepl("parapatricS",x)==TRUE){
    geography<-"parapatric.s"
  }
  
  if(grepl("japan",x)==TRUE){
    geography<-"allopatric.d"
  }
  
  return(geography)
  
}

names(dist.outliers)[1]<- "mean.distance"

outlier.mean.distances <- dist.outliers %>%
  select(study_com, lg, mean.distance) %>%
  group_by(study_com, lg) %>%
  summarise(outlier.mean.distance = mean(mean.distance, na.rm = TRUE)) %>%
  ungroup

null.distances %>%
  filter(study_com == "whtcmn_group1")%>%
   ggplot(aes(x=mean.distance))+
   geom_histogram()+
    facet_wrap(~lg)

snp.file %>%
  filter(study_com == "japan_group1")%>%
  ggplot(aes(x=pos, y=fst, color=fst.outlier))+
  geom_point()+
  facet_wrap(~lg)

# calculate z scores and p values for each outlier mean

dist.df <- null.distances %>%
  group_by(study_com, lg) %>%
  summarise(null.05 = quantile(mean.distance,na.rm=TRUE,probs=0.05)[1], null.mean = mean(mean.distance)) %>%
  left_join(outlier.mean.distances, copy=TRUE) %>%
  mutate(dist.diff = null.mean - outlier.mean.distance, clustered = outlier.mean.distance <= null.05)
  

dist.df$geography <- sapply(dist.df$study_com, add.geo)
  
dist.df %>%
  ggplot(aes(x = study_com, y = log(dist.diff+1), fill = geography))+
  geom_boxplot()+
  facet_grid(~geography, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dist.df %>%
  ggplot(aes(x = study_com, y = log(dist.diff+1), fill = geography))+
  geom_boxplot()+
  facet_grid(~geography, scales = "free_x")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylim(0,0.1)

dist.df %>%
  ggplot(aes(x = geography, y = dist.diff, fill = geography))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dist.summary <- dist.df %>%
  group_by(study_com) %>%
  summarise(avg.clustering = mean(dist.diff))
  
dist.summary$geography <- sapply(dist.summary$study_com, add.geo)

dist.summary %>%
  ggplot(aes(x = geography, y = avg.clustering, fill = geography))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cluster.tally <- dist.df %>%
  group_by(geography, study_com)%>%
  tally(as.numeric(clustered))

cluster.lgs <- dist.df %>%
  group_by(geography, study_com) %>%
  mutate(num.lg = length(unique(lg))) %>%
  ungroup %>%
  select(geography, study_com, num.lg)
  
cluster.tally <- left_join(cluster.tally, cluster.lgs)

cluster.tally <- cluster.tally %>%
  mutate(p.clustered = n/num.lg)

cluster.tally %>%
  ggplot(aes(x = geography, y = p.clustered, color = geography))+
  geom_point()

cluster.tally %>%
  with(lm(p.clustered ~ geography)) %>%
  anova


#############

dist.df  <- read.table(file = "clustering_distances.txt")
cluster.tally<- read.table(file = "clustering_tally.txt")


dist.summary %>%
  ggplot(aes(x = geography, y = avg.clustering, fill = geography))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

  