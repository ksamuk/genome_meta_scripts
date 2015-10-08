#washy washy
rm(list=ls())

library("ggplot2")
library("dplyr")
#library("devtools")
#install_github("karthik/wesanderson")
library("wesanderson")
library("ggthemes")
library("gridExtra")
library("Hmisc")

# choose either fst or fst/dxy 
coeff.dat <- read.table(file = "analysis_ready/75k_stats_model_dxy_fits.txt", header = TRUE, stringsAsFactors = FALSE)
#coeff.dat <- read.table(file = "analysis_ready/75k_stats_model_fst_fits.txt", header = TRUE, stringsAsFactors = FALSE)

#initialize graphics device
dev.off()
dev.new()
quartz(width = 6, height = 6)

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

group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
group.rename <- c("Allopatry\nDivergent", "Allopatry\nParallel", "Gene Flow\nDivergent", "Gene Flow\nParallel")
coeff.dat$group2.new <-group.rename[match(coeff.dat$group2, group.old.names)]
coeff.dat$group.new <- group.rename[match(coeff.dat$group, group.old.names)]


############################## DO ANY GROUPS DIFFER FROM THE "STRONG NULL"

#permute means

#choose group2.new (relaxed), or group.new (strict)
coeff.dat.small <- coeff.dat %>%
  select(group.new, recomb_rate)

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

# create ecdf functions (for plott)
ecdf.df <- permuted.means.df %>% 
  group_by(group) %>%
  do(ecdf = ecdf(.$mean_recomb)) 

ecdf.df$obs <- observed.means$mean_recomb
ecdf.df$p[1] <- ecdf.df$ecdf[1][[1]](ecdf.df$obs[1])
ecdf.df$p[2] <- ecdf.df$ecdf[2][[1]](ecdf.df$obs[2])
ecdf.df$p[3] <- ecdf.df$ecdf[3][[1]](ecdf.df$obs[3])
ecdf.df$p[4] <- ecdf.df$ecdf[4][[1]](ecdf.df$obs[4])

# function for two-sided, monte carlo style pvalue
two_side_p <- function(dist, mean){
	p1 <- (sum(dist > mean)+1) / (length(dist)+1)
	p2 <- (sum(dist < mean)+1) / (length(dist)+1)
	p <- min(p1, p2)*2
	return(p)
}

# calculate pvalues (used in paper)
pvals <- list()
for (i in 1:length(unique(permuted.means.df$group))){
	df <- permuted.means.df %>% filter(group == unique(permuted.means.df$group)[i])
	pvals[[i]] <- data.frame(pvalue = two_side_p(df$mean_recomb, observed.means$mean_recomb[i]), group = unique(permuted.means.df$group)[i])
} 

pvals 

# figure s3

permuted.means.df %>%
  ggplot(aes(x = mean_recomb)) +
  geom_histogram(binwidth = 0.001) +
  geom_segment(data = ecdf.df,aes(x = ecdf.df$obs, xend = ecdf.df$obs, 
  																y = 0, yend = 200,
  																show_guide = F), 
  						 										size = 1, color = "red")+
  facet_wrap(~group, scales = "free_y") +
	theme_hc(base_size = 14) +
	theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"),
				axis.title.y = element_text(vjust=1),
				axis.title.x = element_text(vjust=-1),
				axis.text.x = element_text(hjust = 0.8, vjust = 0.8, angle = 45))+
	xlab("Mean permuted recombination rate coefficient") +
	ylab("Frequency")

#fst relaxed
#quartz.save("figures/figureS5.pdf", type = "pdf")
#fst/dxy relaxed
#quartz.save("figures/figureS6.pdf", type = "pdf")
#fst strict
#quartz.save("figures/figureS7.pdf", type = "pdf")
#fst/dxy strict
quartz.save("figures/figureS8.pdf", type = "pdf")

dev.off()
