## permutation test to assess significance of clustering data

########################################
# Libraries and initalizing variables
########################################

library("ggplot2")
library("dplyr")
library("devtools")
library("wesanderson")
library("grid")
library("gridExtra")
library("viridis")
library("ggthemes")
library("lazyeval")
library("Hmisc")

list.files("shared_functions", full.names = TRUE) %>% sapply(source)

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
select <- dplyr::select

########################################
# Read in and format input files
########################################

## the clustering file
cluster.df <- read.table(file = "analysis_ready/snp_clustering_metrics.txt", header = TRUE, stringsAsFactors = FALSE)
cluster.df <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp)

## add in region data (for looser geography)
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short pop codes
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

# associate pops in coeff dat with regions in region dat
region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop1","reg1")
cluster.df$reg1 <- region.sub$reg1[match(cluster.df$pop1, region.sub$pop1)]

region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop2","reg2")
cluster.df$reg2 <- region.sub$reg2[match(cluster.df$pop2, region.sub$pop2)]

# make new geographic categories
cluster.df$geography2 <- ifelse(cluster.df$reg1==cluster.df$reg2, "para", "allo")

#make new groups 
cluster.df$group2 <- paste0(cluster.df$geography2,"_",cluster.df$ecology)

group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
group.rename <- c("Allopatry\nDivergent", "Allopatry\nParallel", "Gene Flow\nDivergent", "Gene Flow\nParallel")
cluster.df$group2.new <-group.rename[match(cluster.df$group2, group.old.names)]
cluster.df$group.new <- group.rename[match(cluster.df$group, group.old.names)]

########################################
# Permutation tests
########################################

##### Nearest neighbour distance: Relaxed

stat <- "nnd.diff"
group_type <- "group2.new" # relaxed
n_permutations <- 10000

# run the permutation function above for 10,000 interations and bind into df
permuted.means.list <- replicate(n_permutations, 
																 permute_means(cluster.df, stat, group_type), simplify = FALSE)
permuted.means.df <- bind_rows(permuted.means.list)

observed.means <- cluster.df %>% 
	group_by(group) %>% 
	summarise_(mean_stat = interp(~ mean(var, na.rm = TRUE), var = as.name(stat))) %>% 
	ungroup

# ecdfs for plotting
ecdf.df <- permuted.means.df %>% 
	group_by_(interp(group_type)) %>% 
	do_(ecdf = interp(~ ecdf(.$var), var = as.name(stat)))

ecdf.df$obs <- observed.means[,"mean_stat"] %>% unlist
ecdf.df$p[1] <- ecdf.df$ecdf[1][[1]](ecdf.df$obs[1])
ecdf.df$p[2] <- ecdf.df$ecdf[2][[1]](ecdf.df$obs[2])
ecdf.df$p[3] <- ecdf.df$ecdf[3][[1]](ecdf.df$obs[3])
ecdf.df$p[4] <- ecdf.df$ecdf[4][[1]](ecdf.df$obs[4])

# function for two-sided, monte carlo style pvalue


# calculate pvalues (used in paper)
pvals <- list()
group_names <- unlist(unique(permuted.means.df[,1]))
for (i in 1:length(group_names)){
	df <- subset(permuted.means.df, permuted.means.df[,1] == group_names[i])
	pvals[[i]] <- data.frame(pvalue = two_side_p(unlist(df[,2]), observed.means$mean_stat[i]), group = group_names[i])
} 

# plots

# where do the empirical means fall in the permuted distributions
names(permuted.means.df)[1] <- "group"
names(ecdf.df)[1] <- "group"

dev.off()
dev.new()
cairo_pdf(file = "figures/nnd_diff_relaxed_raw.pdf", width =8.5, height = 6)
#quartz(width = 6, height = 6)

permuted.means.df %>%
  ggplot(aes(x = nnd.diff, fill = group)) +
  geom_histogram(binwidth = 0.02) +
  geom_segment(data = ecdf.df,aes(x = obs, xend = obs, y = 400,yend = 0, show_guide = F), 
  						 size = 1, color = "black", arrow = arrow(length = unit(0.3, "cm"), type = "closed"))+
  scale_fill_manual(values = pal)+
  theme_classic(base_size = 16)+
  theme(strip.text.x = element_blank(), 
        strip.background = element_blank(), 
        axis.title.x = element_text(vjust=-1.5), 
        axis.title.y = element_text(vjust=1.5),
        plot.margin=unit(c(1,1,1.5,1.2),"cm"))+
  facet_grid(~group)+
  xlab("Mean Expected NND - Outlier NND (cM)")+
  ylab("Frequency")

#quartz.save("figures/nnd_diff_relaxed_raw.pdf", type = "pdf")
dev.off()




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

cluster.df$group2.old <- cluster.df$group2
group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
group.rename <- c("Allopatric D", "Allopatric S", "Parapatric D", "Parapatric S")
cluster.df$group2 <- group.rename[match(cluster.df$group2, group.old.names)]


#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
size <- 16
theme.all <- theme(legend.position="none", 
                   axis.title.x = element_blank(), 
                   axis.title.y = element_text(vjust=1.5),
                   axis.text.x = element_blank(),
                   axis.line.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text = element_blank(), 
                   strip.background = element_blank(),
                   legend.title=element_blank())

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
      coord_cartesian(ylim=c(0,0.8))+
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
    geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, width=1, color = 1, linetype = 1)+
    theme_classic(base_size = size) +
    coord_cartesian(ylim=c(-5,5))+
    theme.all+
    ylab("Expected NND - Outlier NND (cM)")

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
#rep.df <- bind_rows(rep.df)

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
rep.nnd.df$group2 <- group.rename[match(rep.nnd.df$group2, group.old.names)]

nnd.rep <- rep.nnd.df %>%
  ggplot(aes(x = group2, y = gen.pos, color = group2, fill = group2)) +
  geom_jitter(position = position_jitter(height = 2), size = 1.5)+
    scale_color_manual(values = pal) +   
    scale_fill_manual(values = pal) +
    theme_classic(base_size = size) +
    theme.all+
    coord_cartesian(ylim = c(0,80))+
    ylab("Map position on LG 4 (cM)")

grid.arrange(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep, ncol = 2)
grid_arrange_shared_legend(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep)
