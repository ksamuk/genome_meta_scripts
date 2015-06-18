# calculate nearest neighbour distances

library(data.table)
library(dplyr)
library(ggplot2)
library(microbenchmark)
library(parallel)


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

# Old distance method -----------------------------------------------------

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


############### NEW APPROACH: PERMUTATION FTW


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


sample_sites <- function (gen.pos, num.outliers) {
  
  site.sample <- gen.pos %>% sample_n(num.outliers) %>% arrange(gen.pos) %>% select(gen.pos) %>% with(as.numeric(gen.pos)) %>% mean(na.rm = TRUE)
  return(site.sample)
  
}
  
permute_distances_snp <- function(x, num.samples){
  
  mean.distances.df <- list()
  
  for (i in unique(x$study_com)){
    
    print(paste0("processing ", i,"..."))
    dist.sub <- subset(x, x$study_com == i)
    mean.distances.lg.all <- list()
    
    for (j in unique(snp.outlier.sub$lg)){
      
      dist.sub.lg <- dist.sub %>%
        filter(lg == j) %>%
        filter(!is.na(gen.pos)) %>%
        arrange(gen.pos)
      
      num.outliers <- sum(as.numeric(dist.sub.lg$fst.outlier), na.rm = TRUE)
      
      if (num.outliers > 0){
        
        mean.distances.lg <- replicate(num.samples, sample_sites(dist.sub.lg, num.outliers))
        mean.distances.lg.df <- data.frame(mean.distance = mean.distances.lg)
        mean.distances.lg.df$lg <- j
        mean.distances.lg.all[[j]] <- mean.distances.lg.df
      }
      
    }
    mean.distances.lg.all <- do.call("rbind", mean.distances.lg.all)
    mean.distances.lg.all$study_com = i
    mean.distances.df[[i]] <- mean.distances.lg.all
    
  }
  mean.distances.df <- do.call("rbind", mean.distances.df)
  rownames(mean.distances.df) <- NULL
  return(mean.distances.df)
  
}

distances.null <- permute_distances_snp(snp.file, 1000)

permute_distances_snp_list <- function(x, num.samples){
  
  mean.distances.lg.all <- list()  
  
    for (j in unique(x$lg)){
      
      dist.sub.lg <- x %>%
        filter(lg == j) %>%
        filter(!is.na(gen.pos)) %>%
        arrange(gen.pos)
      
      num.outliers <- sum(as.numeric(dist.sub.lg$fst.outlier), na.rm = TRUE)
      
      if (num.outliers > 0){
        
        mean.distances.lg <- replicate(num.samples, sample_sites(dist.sub.lg, num.outliers))
        mean.distances.lg.df <- data.frame(mean.distance = mean.distances.lg)
        mean.distances.lg.df$lg <- j
        mean.distances.lg.all[[j]] <- mean.distances.lg.df
      }
      
    }
    mean.distances.lg.all <- do.call("rbind", mean.distances.lg.all)
    mean.distances.lg.all$study_com = i
    return(mean.distances.lg.all)
  
}

split.df <- split(snp.file, snp.file$study_com)  
mclapply(split.df, permute_distances_snp_list, num.samples = 10)


