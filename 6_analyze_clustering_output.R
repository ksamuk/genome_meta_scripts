## permutation test to assess significance of clustering data

library("ggplot2")
library("dplyr")
library("devtools")
#install_github("karthik/wesanderson")
#evtools::install_github("sjmgarnier/viridis")
library("wesanderson")
library("grid")
library("gridExtra")
library("viridis")
library("ggthemes")

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
select <- dplyr::select

## libraries
library("dplyr")

## the clustering file
cluster.df <- read.table(file = "analysis_ready/snp_clustering_metrics.txt", header = TRUE, stringsAsFactors = FALSE)
cluster.df <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp)

#### Do the groups differ in the number of significantly NND clustered chromosomes?
#### STRICT

#permute means
coeff.dat.small <- cluster.df %>%
  select(group, nnd.diff)

names(coeff.dat.small)[1] <- "group"

# function that shuffles groups, calculates their mean recom coeff, and returns the latter
permute_means <- function(data) {
  data$group <- sample(data$group, length(data$group))
  mean <- data %>% group_by(group) %>% summarise(mean.nnd.diff = mean(nnd.diff, na.rm = TRUE)) %>% ungroup
  return(mean)
}

# run the funciton above for 10,000 interations and bind into df
permuted.means.list <- replicate(50000, permute_means(coeff.dat.small), simplify = FALSE)
permuted.means.df <- bind_rows(permuted.means.list)
observed.means <- coeff.dat.small  %>% group_by(group) %>% summarise(mean.nnd.diff = mean(nnd.diff, na.rm = TRUE)) %>% ungroup

ecdf.df <- permuted.means.df %>% 
  group_by(group) %>%
  do(ecdf = ecdf(.$mean.nnd.diff)) 

ecdf.df$obs <- observed.means$mean.nnd.diff
ecdf.df$p[1] <- ecdf.df$ecdf[1][[1]](ecdf.df$obs[1])
ecdf.df$p[2] <- ecdf.df$ecdf[2][[1]](ecdf.df$obs[2])
ecdf.df$p[3] <- ecdf.df$ecdf[3][[1]](ecdf.df$obs[3])
ecdf.df$p[4] <- ecdf.df$ecdf[4][[1]](ecdf.df$obs[4])

# plots

# where do the empirical means fall in the permuted distributions
permuted.means.df %>%
  ggplot(aes(x = mean.nnd.diff, fill = group)) +
  geom_histogram() +
  geom_segment(data=ecdf.df,aes(x = ecdf.df$obs, xend = ecdf.df$obs, y = 0,yend = 7500,show_guide = F), size = 1, color = "black")+
  scale_fill_manual(values = pal)+
  theme_classic(base_size = 16)+
  theme(strip.text.x = element_blank(), 
        strip.background = element_blank(), 
        axis.title.x = element_text(vjust=-1.5), 
        axis.title.y = element_text(vjust=1.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  facet_wrap(~group, scales = "free_x")+
  xlab("Mean Expected NND - Outlier NND (cM)")+
  ylab("Frequency")


#### Do the groups differ in the number of significantly NND clustered chromosomes?
#### RELAXED

#permute means
coeff.dat.small <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.emp - nnd.mean.null) %>%
  select(group2, nnd.diff)

names(coeff.dat.small)[1] <- "group"

# function that shuffles groups, calculates their mean recom coeff, and returns the latter
permute_means <- function(data) {
  data$group <- sample(data$group, length(data$group))
  mean <- data %>% group_by(group) %>% summarise(mean_zscore = mean(nnd.diff, na.rm = TRUE)) %>% ungroup
  return(mean)
}

# run the funciton above for 10,000 interations and bind into df
permuted.means.list <- replicate(10000, permute_means(coeff.dat.small), simplify = FALSE)
permuted.means.df <- bind_rows(permuted.means.list)
observed.means <- coeff.dat.small  %>% group_by(group) %>% summarise(mean_zscore = mean(nnd.diff, na.rm = TRUE)) %>% ungroup

ecdf.df <- permuted.means.df %>% 
  group_by(group) %>%
  do(ecdf = ecdf(.$mean_zscore)) 

ecdf.df$obs <- observed.means$mean_zscore
ecdf.df$p[1] <- ecdf.df$ecdf[1][[1]](ecdf.df$obs[1])
ecdf.df$p[2] <- ecdf.df$ecdf[2][[1]](ecdf.df$obs[2])
ecdf.df$p[3] <- ecdf.df$ecdf[3][[1]](ecdf.df$obs[3])
ecdf.df$p[4] <- ecdf.df$ecdf[4][[1]](ecdf.df$obs[4])

# plots

# where do the empirical means fall in the permuted distributions
permuted.means.df %>%
  ggplot(aes(x = mean_zscore)) +
  geom_histogram() +
  geom_segment(data=ecdf.df,aes(x = ecdf.df$obs, xend = ecdf.df$obs, y = 0,yend = 2000,show_guide = F), size = 1, color = "red")+
  facet_wrap(~group, scales = "free_x")


#grouped 
grouped.df <- cluster.df %>%
  mutate(group = paste0(geography, "_", ecology)) %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01) %>%
  group_by(group2) %>%
  summarise(num.cluster = mean(nnd.sig, na.rm = TRUE)) 

# dem plots?

#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
size <- 16
theme.all <- theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_text(vjust=1.5))

# counts of clustered lgs
prop.clustered <- cluster.df %>%
  mutate(more.clustered = nnd.emp.percentile*2 < 0.05) %>%
  mutate(less.clustered = nnd.emp.percentile*2 > 0.95) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  select(comparison, lg, group, group2, more.clustered, less.clustered) %>%
  group_by(group2, comparison) %>%
  summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
  mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
  filter(num.lg > 1) %>%
    ggplot(aes(x = group2, y = more.clustered, fill = group2, color = group2))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
      geom_jitter(size = 1)+
      geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, weight=1, color = 1, linetype = 1)+
      #theme(panel.grid = element_blank())+
      theme_classic(base_size = size)+
      theme.all+
      ylab("Prop. chromosomes clustered")

# coeff of disp
coeff.dispersion <- cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2, comparison) %>%
  summarise(disp.out = mean(disp.out))%>%
    ggplot(aes(x = group2, y = disp.out, fill = group2, color = group2))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    geom_jitter(size = 1)+
    geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, weight=1, color = 1, linetype = 1)+
    theme_classic(base_size = size)+
    coord_cartesian(ylim=c(0,25))+
    theme.all+
    ylab("Outlier dispersion coefficient")

# nnd.diff

nnd.diff <- cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2, comparison) %>%
  summarise(nnd.diff = mean(nnd.diff))%>%
    ggplot(aes(x = group2, y = nnd.diff, fill = group2, color = group2))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    geom_jitter(size = 1)+
    geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, weight=1, color = 1, linetype = 1)+
    theme_classic(base_size = size) +
    coord_cartesian(ylim=c(-5,5))+
    theme.all+
    ylab("Expected NND - Outlier NND (cM)")

grid.arrange(prop.clustered, nnd.diff, coeff.dispersion, ncol = 2)


cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  #group_by(group2, comparison) %>%
  #summarise(nnd.diff = mean(nnd.diff), num.outliers = num.outliers)%>%
  filter(group2 == "para_S", lg == 4) %>% data.frame

### representative chromosome plot

# read in rep snp files

is.outlier <- function(x){
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  return(x >= x95)
}

rep.files <- list.files("stats/snp_representative", full.names = TRUE)

source("3_process_snp_clustering.R")

extract_lg <- function (file){
  name.split <- strsplit(file,split = "/") 
  name.split <- name.split[[1]][3]
  name.split <- strsplit(name.split,split = "[.]") %>% unlist
  comparison <- paste0(name.split[1],".",name.split[2])
  group <- paste0(name.split[3],"_",name.split[4])
  dat.tmp <- read.table(file = file, header = TRUE, stringsAsFactors = FALSE)
  
  dat.tmp <- dat.tmp %>%
    filter(CHROM == "groupIV") %>%
    filter(!is.infinite(Fst), !is.na(Fst), Fst > 0) %>%
    select(CHROM, POS, Fst)
  
  names(dat.tmp)[1:3] <- c("lg","pos","fst")
  dat.tmp$comparison <- comparison
  dat.tmp$group <- group
  dat.tmp$lg <- 4
  dat.tmp$fst.outlier <- is.outlier(dat.tmp$fst)
  dat.tmp <- dat.tmp %>%
    filter(fst.outlier == TRUE)
  
  #add map distances (generates some warnings, but works as intended)
  dat.tmp <- add_map_distance(dat.tmp)
  return(dat.tmp)
}

rep.df <- lapply(rep.files, extract_lg)
rep.df <- bind_rows(rep.df)

# calculate empirical nnd 

calc_emp_nnd_dist <- function(stats.file.lg){
  
  stats.file.lg <- stats.file.lg %>%
    filter(!is.na(gen.pos))
  
  site.sample <- stats.file.lg %>%
    filter(!is.na(gen.pos)) %>%
    select(gen.pos) %>%
    arrange(gen.pos) %>%
    mutate(dist.1 = c(NA,diff(gen.pos))) %>%
    mutate(dist.2 = c(diff(sort(gen.pos)),NA))
  
  nn.dist <- rep(NA, length(site.sample$genpos))
  for (k in 1:length(site.sample$gen.pos)){
    
    if(!is.na(site.sample$dist.1[k]) & !is.na(site.sample$dist.2[k])){
      nn.dist[k] <- min(c(site.sample$dist.1[k],site.sample$dist.2[k]))
    }else if(is.na(site.sample$dist.1[k])){
      nn.dist[k] <- site.sample$dist.2[k]
    } else if(is.na(site.sample$dist.2[k])){
      nn.dist[k] <- site.sample$dist.1[k]
    }
  }
  
  stats.file.lg$nnd <- nn.dist
  
  return(stats.file.lg)
}

rep.nnd.df <- lapply(rep.df, calc_emp_nnd_dist)
rep.nnd.df <- bind_rows(rep.nnd.df)

rep.nnd.df$group2 <- ifelse(rep.nnd.df$comparison == "cr_stream.wc_stream", "para_S", rep.nnd.df$group)

rep.nnd.df %>%
  #filter(nnd > 0.01) %>%
  ggplot(aes(x = gen.pos, y = 1, color = nnd)) +
  scale_color_viridis()+
  #geom_histogram(binwidth = 0.01)  +
  geom_jitter(position = position_jitter(height = .1), size = 3)+
    facet_grid(group2~.) +
    coord_cartesian(ylim = c(0.8,1.2))+
    theme_classic(base_size = size) +
    theme(strip.text.x = element_blank(), 
        strip.background = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x = element_text(vjust=1.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"))



