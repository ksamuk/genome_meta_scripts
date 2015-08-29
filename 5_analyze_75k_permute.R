#washy washy
rm(list=ls())

library("ggplot2")
library("dplyr")
#library("devtools")
#install_github("karthik/wesanderson")
library("wesanderson")

## can start fresh here
rm(list=ls())
coeff.dat <- read.table(file = "analysis_ready/75k_stats_model_dxy_fits.txt", header = TRUE, stringsAsFactors = FALSE)
coeff.dat <- read.table(file = "analysis_ready/75k_stats_model_fits.txt", header = TRUE, stringsAsFactors = FALSE)

## add in region data (for looser geography)
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short pop codes
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

# associate pops in coeff dat with regions in region dat
region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop1","reg1")
coeff.dat$reg1 <- region.sub$reg1[match(coeff.dat$pop1, region.sub$pop1)]

region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop2","reg2")
coeff.dat$reg2 <- region.sub$reg2[match(coeff.dat$pop2, region.sub$pop2)]

# make new geographic categories
coeff.dat$geography2 <- ifelse(coeff.dat$reg1==coeff.dat$reg2, "para", "allo")

#make new groups :o
coeff.dat$group2 <- paste0(coeff.dat$geography2,"_",coeff.dat$ecology)


############################## DO ANY GROUPS DIFFER FROM THE "STRONG NULL"

#permute means
coeff.dat.small <- coeff.dat %>%
  select(group2, recomb_rate)

names(coeff.dat.small)[1] <- "group"

# function that shuffles groups, calculates their mean recom coeff, and returns the latter
permute_means <- function(data) {
  data$group <- sample(data$group, length(data$group))
  mean <- data %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup
  return(mean)
}

# run the funciton above for 100,000 interations and bind into df
permuted.means.list <- replicate(10000, permute_means(coeff.dat.small), simplify = FALSE)
permuted.means.df <- bind_rows(permuted.means.list)
observed.means <- coeff.dat.small  %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup

ecdf.df <- permuted.means.df %>% 
  group_by(group) %>%
  do(ecdf = ecdf(.$mean_recomb)) 

ecdf.df$obs <- observed.means$mean_recomb
ecdf.df$p[1] <- ecdf.df$ecdf[1][[1]](ecdf.df$obs[1])
ecdf.df$p[2] <- ecdf.df$ecdf[2][[1]](ecdf.df$obs[2])
ecdf.df$p[3] <- ecdf.df$ecdf[3][[1]](ecdf.df$obs[3])
ecdf.df$p[4] <- ecdf.df$ecdf[4][[1]](ecdf.df$obs[4])

# plots

# where do the empirical means fall in the permuted distributions
permuted.means.df %>%
  ggplot(aes(x = mean_recomb)) +
  geom_histogram() +
  geom_segment(data=ecdf.df,aes(x = ecdf.df$obs, xend = ecdf.df$obs, y = 0,yend = 5000,show_guide = F), size = 1, color = "red")+
  facet_wrap(~group, scales = "free_x")


############################## END DO ANY GROUPS DIFFER FROM THE "STRONG NULL"

############################## DO GROUPS DIFFER FROM ONE ANOTHER?

#permute means
coeff.dat.small <- coeff.dat %>%
  select(group, recomb_rate)

permute_mean_differences <- function(data) {
  data$group <- sample(data$group, length(data$group))
  mean <- data %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup
  mean.diffs <- data.frame(group1 = c(rep("allo_D", 3), rep("allo_S", 2), rep("para_D" ,1)), 
                           group2 = c("allo_S", "para_D", "para_S", "para_D", "para_S", "para_S"))
  mean.diffs$mean1 <- mean$mean_recomb[match(mean.diffs$group1, mean$group)]
  mean.diffs$mean2 <- mean$mean_recomb[match(mean.diffs$group2, mean$group)]
  mean.diffs$diff <- mean.diffs$mean1 - mean.diffs$mean2 
  return(mean.diffs[,c(1,2,5)])
}

empirical_mean_differences <- function(data) {
  mean <- data %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup
  mean.diffs <- data.frame(group1 = c(rep("allo_D", 3), rep("allo_S", 2), rep("para_D" ,1)), 
                           group2 = c("allo_S", "para_D", "para_S", "para_D", "para_S", "para_S"))
  mean.diffs$mean1 <- mean$mean_recomb[match(mean.diffs$group1, mean$group)]
  mean.diffs$mean2 <- mean$mean_recomb[match(mean.diffs$group2, mean$group)]
  mean.diffs$diff <- mean.diffs$mean1 - mean.diffs$mean2
  return(mean.diffs[,c(1,2,5)])
}

# run the funciton above for 100,000 interations and bind into df
permuted.means.list <- replicate(100000, permute_mean_differences(coeff.dat.small), simplify = FALSE)
permuted.means.df <- bind_rows(permuted.means.list)

permuted.means.df <- permuted.means.df %>%
  mutate(comparison = paste0(group1,"_",group2))

empirical.means.df <- empirical_mean_differences(coeff.dat.small)

empirical.means.df  <- empirical.means.df  %>%
  mutate(comparison = paste0(group1,"_",group2))

permuted.means.df %>%
  ggplot(aes(x = diff))+
  geom_histogram(binwidth = 0.01)+
  geom_segment(data=empirical.means.df,aes(x = empirical.means.df$diff, xend = empirical.means.df$diff, y = 0,yend = 500,show_guide = F), size = 1, color = "red")+
  facet_wrap(~comparison, scales = "free_y")

observed.means <- coeff.dat.small  %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup
permute_means <- function(data) {
  data$group <- sample(data$group, length(data$group))
  mean <- data %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup
  return(mean)
}

# run the funciton above for 100,000 interations and bind into df
permuted.means.list <- replicate(10000, permute_means(coeff.dat.small), simplify = FALSE)
permuted.means.df <- do.call("rbind",permuted.means.list)
observed.means <- coeff.dat.small  %>% group_by(group) %>% summarise(mean_recomb = mean(recomb_rate)) %>% ungroup

############################## DO GROUPS DIFFER FROM ONE ANOTHER?