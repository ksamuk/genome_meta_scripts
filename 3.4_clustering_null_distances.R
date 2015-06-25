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
  
  site.sample <- gen.pos %>% 
    sample_n(num.outliers)%>% 
    arrange(gen.pos) %>% 
    select(gen.pos) %>%
    mutate(dist.1 = c(NA,diff(gen.pos))) %>%
    mutate(dist.2 = c(diff(sort(gen.pos)),NA))
  
  nn.dist <- rep(NA, length(site.sample$genpos))
  for (i in 1:length(site.sample$gen.pos)){
    
    if(!is.na(site.sample$dist.1[i]) & !is.na(site.sample$dist.2[i])){
      nn.dist[i] <- min(c(site.sample$dist.1[i],site.sample$dist.2[i]))
    }else if(is.na(site.sample$dist.1[i])){
      nn.dist[i] <- site.sample$dist.2[i]
    } else if(is.na(site.sample$dist.2[i])){
      nn.dist[i] <- site.sample$dist.1[i]
    }
  }
  
  return(mean(nn.dist))
  
}

permute_distances_snp_list <- function(x, num.samples){
  
  cat(paste0("permuting null distances for ",x$study_com[1]))
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
null.distances <- mclapply(split.df, permute_distances_snp_list, num.samples = 10000, mc.cores = 6, mc.silent = FALSE)

null.distances <- do.call("rbind", null.distances)
rownames(null.distances) <- NULL
write.table(null.distances, file="null_snp_distances.txt", row.names = FALSE)

  