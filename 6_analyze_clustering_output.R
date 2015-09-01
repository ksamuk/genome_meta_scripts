## permutation test to assess significance of clustering data

library("ggplot2")
library("dplyr")
#library("devtools")
#install_github("karthik/wesanderson")
library("wesanderson")

select <- dplyr::select

## libraries
library("dplyr")

## the clustering file
cluster.df <- read.table(file = "analysis_ready/snp_clustering_metrics.txt", header = TRUE, stringsAsFactors = FALSE)

#### Do the groups differ in the number of significantly NND clustered chromosomes?
#### STRICT

#permute means
coeff.dat.small <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.emp - nnd.mean.null) %>%
  select(group, nnd.diff)

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
  geom_segment(data=ecdf.df,aes(x = ecdf.df$obs, xend = ecdf.df$obs, y = 0,yend = 5000,show_guide = F), size = 1, color = "red")+
  facet_wrap(~group, scales = "free_x")


# elaborating groupings
cluster.df <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp)

#grouped 
grouped.df <- cluster.df %>%
  mutate(group = paste0(geography, "_", ecology)) %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01) %>%
  group_by(group2) %>%
  summarise(num.cluster = mean(nnd.sig, na.rm = TRUE)) 

# dem plots?
cluster.df %>%
  ggplot(aes(x = group2, y = log(nnd.diff+1))) +
  geom_boxplot()

cluster.df %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01)%>%
  filter(nnd.sig == TRUE) %>%
  ggplot(aes(x = group, y = log(nnd.diff+1))) +
  geom_jitter()

cluster.df %>%
  ggplot(aes(x = group2, y = disp.out+1)) +
  geom_boxplot()

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
  geom_segment(data=ecdf.df,aes(x = ecdf.df$obs, xend = ecdf.df$obs, y = 0,yend = 5000,show_guide = F), size = 1, color = "red")+
  facet_wrap(~group, scales = "free_x")


# elaborating groupings
cluster.df <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp)

#grouped 
grouped.df <- cluster.df %>%
  mutate(group = paste0(geography, "_", ecology)) %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01) %>%
  group_by(group2) %>%
  summarise(num.cluster = mean(nnd.sig, na.rm = TRUE)) 

# dem plots?

# counts of clustered lgs
cluster.df %>%
  mutate(more.clustered = nnd.emp.percentile*2 < 0.05) %>%
  mutate(less.clustered = nnd.emp.percentile*2 > 0.95) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  select(comparison, lg, group, group2, more.clustered, less.clustered) %>%
  group_by(group2, comparison) %>%
  summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
  mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
  filter(num.lg > 1) %>%
    ggplot(aes(x = group2, y = more.clustered))+
      geom_boxplot()

cluster.df %>% 
  mutate(disp.out = nnd.emp.percentile*2 < 0.05) %>%
  mutate(less.clustered = nnd.emp.percentile*2 > 0.95) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  select(comparison, lg, group, group2, more.clustered, less.clustered) %>%
  group_by(group2, comparison) %>%
  summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
  mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
  filter(num.lg > 1) %>%
  ggplot(aes(x = group2, y = more.clustered))+
  geom_boxplot()


cluster.df %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2, comparison) %>%
  summarise(disp.out = mean(disp.out))%>%
  ggplot(aes(x = group2, y = disp.out)) +
  geom_boxplot()


cluster.df %>%
  filter(num.outliers >= 2) %>%
  #filter(group2 == "para_D")
  ggplot(aes(x = log(num.outliers), y = nnd.emp.zscore, color = group2)) +
  geom_jitter(position = position_jitter(width = .1, height = 0)) +
  geom_smooth(se = FALSE, size = 3)+
  theme_classic()+
  coord_cartesian(xlim=c(0,4))

cluster.df %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01)%>%
  filter(nnd.sig == TRUE) %>%
  ggplot(aes(x = group, y = nnd.diff+)) +
  geom_jitter()

cluster.df %>%
  ggplot(aes(x = group2, y = disp.out+1)) +
  geom_boxplot()