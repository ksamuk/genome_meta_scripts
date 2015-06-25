# calculate nearest neighbour distances

library(data.table)
library(dplyr)
library(parallel)
library(ggplot2)

# read in snp file
snp.file <- "stats/snp_filtered/stats_master_variant_2015-06-11.filt.txt"
snp.file <- fread(snp.file)

# remove weird outliers
snp.file[!is.na(snp.file$fst) & snp.file$fst<0, ]$fst <- NA

# find outliers

is.outlier<-function(x){
  return(x >= quantile(x,na.rm=TRUE,probs=0.95)[1])
}

snp.file$fst <- as.numeric(snp.file$fst)

snp.file<-snp.file%>%
  group_by(study_com)%>%
  mutate(fst.outlier = is.outlier(fst)) %>% 
  ungroup

str(snp.file)

############### NEW APPROACH: PERMUTATION FTW

nndist.lg <- function (lg.df) {
  
  site.sample <- lg.df %>%
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


for (i in unqiue(snp.file$study_com)){
  snp.sub <- subset(snp.file)
}

split.df <- split(snp.file, snp.file$study_com)  
null.distances <- mclapply(split.df, permute_distances_snp_list, num.samples = 10000, mc.cores = 6, mc.silent = FALSE)

null.distances <- do.call("rbind", null.distances)
rownames(null.distances) <- NULL
write.table(null.distances, file="null_snp_distances.txt", row.names = FALSE)

# GRAB EMMPIRICAL NNDs

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
    geography<-"parapatric.d"
  }
  
  return(geography)
  
}

nndist.df<- read.table("test.distances.txt",stringsAsFactors = FALSE,header=TRUE)
nndist.df$geography <- sapply(nndist.df$study_com, add.geo)

melt(nndist, id.vars = c("fst.outlier"))

nndist.df <- nndist.df %>%
  filter(!is.na(fst.outlier))

nn.true <- subset(nndist.df, nndist.df$fst.outlier ==TRUE) %>% select(study_com,lg,nndist)
names(nn.true)[3]<-"nn.true"
nn.false <- subset(nndist.df, nndist.df$fst.outlier ==FALSE) %>% select(study_com,lg,nndist)
names(nn.false)[3]<-"nn.false"

nn.matched <- left_join(nn.true,nn.false)
nn.matched$diff <- nn.matched$nn.false/nn.matched$nn.true

nn.matched$geography <- sapply(nn.matched$study_com, add.geo)

nn.matched %>%
  filter(diff < 100) %>%
  ggplot(aes(x = study_com, y = log(diff+1), fill = geography))+
  geom_boxplot()+
  facet_grid(~geography, scales = "free_x")
  
  


  