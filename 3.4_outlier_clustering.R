# calculate nearest neighbour distances

library(data.table)
library(dplyr)
library(ggplot2)


# read in snp file
snp.file <- list.files(file.path("stats/snp_filtered"),pattern=".txt$", full.names = TRUE)
snp.file <- fread(snp.file)

# remove weird outliers
snp.file[!is.na(snp.file$fst) & snp.file$fst<0, ]$fst <- NA
hist(snp.file$fst)

# find outliers

is.outlier<-function(x){
  return(x >= quantile(x,na.rm=TRUE,probs=0.95)[1])
}

snp.file<-snp.file%>%
  group_by(study_com)%>%
  mutate(fst.outlier = is.outlier(fst)) %>% 
  ungroup

snp.file%>%
  filter(study_com == "allopatricS_boos.cons") %>%
  ggplot(aes(x = fst)) +
  geom_histogram()

snp.outliers <- snp.file %>%
  filter(fst.outlier == TRUE) 

snp.not.outliers <- snp.file %>%
  filter(fst.outlier == FALSE) 

# calc distances between outliers

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
  
  return(geography)
  
}

dist.outliers$geography <- sapply(dist.outliers$study_com,add.geo)

dist.outliers$outlier <- TRUE

################### NOT OUTLIERS

dist.not.outliers <- data.frame()

for (i in unique(snp.outliers$study_com)){
  print(paste0("processing ", i,"..."))
  snp.outlier.sub <- subset(snp.not.outliers , snp.not.outliers$study_com == i)
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
  
  dist.not.outliers <- rbind(dist.not.outliers, diff.sub)
  
}

# inject geography

dist.not.outliers$geography <- sapply(dist.not.outliers$study_com,add.geo)

dist.not.outliers$outlier <- FALSE

dist.snps <- rbind(dist.outliers, dist.not.outliers)

dist.snps <- dist.snps %>%
  arrange(study_com, lg)

# bind snp distance types


##################### PLOTS

  
dist.snps  %>%
  ggplot(aes(x = study_com, y = outlier.dist, fill = geography, color = outlier)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_grid(~geography, scale="free_x")

snp.dist.summary %>%
  ggplot(aes(x = study_com, y = dispersion, color = geography)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

dist.snps  %>%
  ggplot(aes(x = geography, y = log(outlier.dist), fill = outlier)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



############### NEW APPROACH

dist.distributions <- data.frame()

for (i in unique(snp.outliers$study_com)){
  print(paste0("processing ", i,"..."))
  snp.outlier.sub <- subset(snp.not.outliers , snp.not.outliers$study_com == i)
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
  
  dist.not.outliers <- rbind(dist.not.outliers, diff.sub)
  
}


permute_distances <- function(x, num.samples){
  
  mean.distances.df <- data.frame()
  
  for (i in unique(x$study_com)){
    
    print(paste0("processing ", i,"..."))
    dist.sub <- subset(x, x$study_com == i)
    mean.distances.lg.all <- data.frame()
    for (j in unique(snp.outlier.sub$lg)){
      
      dist.sub.lg <- dist.sub %>%
        filter(lg == j) 

      num.outliers <- sum(as.numeric(dist.sub.lg$outlier))
      mean.distances.lg <- replicate(num.samples, mean(sample(dist.sub.lg$outlier.dist[!is.na(dist.sub.lg$outlier.dist)], num.outliers)))
      
      mean.distances.lg.df <- data.frame(mean.distance = mean.distances.lg)
      mean.distances.lg.df$lg <- j
      
      mean.distances.lg.all <- rbind(mean.distances.lg.all, mean.distances.lg.df)
      
    }
    mean.distances.lg.all$study_com = i
    mean.distances.df <- rbind(mean.distances.df, mean.distances.lg.all)
    
  }
  
  return(mean.distances.df)
  
}


sample_sites <- . %>% sample_n(num.outliers) %>% arrange(gen.pos) %>% select(gen.pos) %>% with(as.numeric(gen.pos)) %>% mean(na.rm = TRUE)
  
permute_distances_snp <- function(x, num.samples){
  
  mean.distances.df <- data.frame()
  
  for (i in unique(x$study_com)){
    
    print(paste0("processing ", i,"..."))
    dist.sub <- subset(x, x$study_com == i)
    mean.distances.lg.all <- data.frame()
    for (j in unique(snp.outlier.sub$lg)){
      
      dist.sub.lg <- dist.sub %>%
        filter(lg == j) %>%
        filter(!is.na(gen.pos)) %>%
        arrange(gen.pos)
      
      num.outliers <- sum(as.numeric(dist.sub.lg$fst.outlier), na.rm = TRUE)
      
      if (num.outliers > 0){
        
        mean.distances.lg <- replicate(num.samples, sample_sites(dist.sub.lg))
        
        mean.distances.lg.df <- data.frame(mean.distance = mean.distances.lg)
        mean.distances.lg.df$lg <- j
        
        mean.distances.lg.all <- rbind(mean.distances.lg.all, mean.distances.lg.df)
      }
      
    }
    mean.distances.lg.all$study_com = i
    mean.distances.df <- rbind(mean.distances.df, mean.distances.lg.all)
    
  }
  
  return(mean.distances.df)
  
}

distance.null <- permute_distances_snp(snp.file, 10)

distance.null %>%
  filter(study_com == "benlim_lq") %>%
  ggplot(aes(x = mean.distance)) +
  geom_histogram()+
  facet_wrap(~lg)
