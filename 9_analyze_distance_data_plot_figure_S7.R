############################################################
# inspect fst vs. geographic distance relationship
# KS Aug 2015
############################################################

rm(list =ls())

#########################
# LIBRARIES & FUNCTIONS
#########################

library("dplyr")
library("tidyr")
library("broom")
library("ggplot2")
library("ggthemes")
library("Hmisc")
library("cowplot")
library("lazyeval")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

#########################
# INPUT FILES
#########################

# df containing mean fst values for all populations
fst.df <- read.table("meta_data/fst_df.txt", header = TRUE, stringsAsFactors = FALSE)

# dfs containing linear model coefficients for outliers vs. recomb rate
fst.model.fits <- initialize_coeff_dat_files()

# df containing distances (least cost, great circle) between populations
# also measures of isolation from ocean (dist to coast, i.e. how inland, for each population)
dist.df <- read.table("meta_data/pop_geo_distances.txt", header = TRUE, stringsAsFactors = FALSE)

#########################
# PRE PROCESSING
#########################

# bind the three data files together
all.df <- cbind(fst.model.fits, dist.df[,7:10], fst.df$fst)
names(all.df)[length(all.df)] <- "fst"

# add in combined measures of allopatry (simple transformations of distance data)
all.df <- all.df %>% 
	mutate(isolation.sum = dist.to.coast1 + dist.to.coast2)%>%
	mutate(comparison = paste0(pop1, pop2))%>%
	rowwise() %>%
	mutate(isolation.max = max(dist.to.coast1, dist.to.coast2))

all.df$ecology <- gsub("D", "Divergent", all.df$ecology)
all.df$ecology <- gsub("S", "Parallel", all.df$ecology)

##################################
# Size / theme for figs s7 - s9
##################################

size_all <- 2
alpha_all <- 0.2

theme_all <- theme_hc(base_size = 16) + 
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.text.x = element_text(angle = 45, hjust = 1)) 

##################################
# Figures S7 - S9
##################################

# S7
pdf(file = "figures/figureS7.pdf", height = 8.5, width = 8.5, onefile = FALSE)

coeff_type <- "recomb_rate_fst"
coeff_name <- "Recombination rate vs. \nFST outlier coeffficient"
plot_coef_vs_distance(all.df, coeff_type, coeff_name,size_all, alpha_all, theme_all)

dev.off()

# S8
pdf(file = "figures/figureS8.pdf", height = 8.5, width = 8.5, onefile = FALSE)

coeff_type <- "recomb_rate_dxy"
coeff_name <- "Recombination rate vs. \nDXY outlier coeffficient"
plot_coef_vs_distance(all.df, coeff_type, coeff_name,size_all, alpha_all, theme_all)

dev.off()

# S9
pdf(file = "figures/figureS9.pdf", height = 8.5, width = 8.5, onefile = FALSE)

coeff_type <- "recomb_rate_fst_dxy"
coeff_name <- "Recombination rate vs. \nFST/DXY outlier coeffficient"
plot_coef_vs_distance(all.df, coeff_type, coeff_name,size_all, alpha_all, theme_all)

dev.off()

#########################
# Permutations re: above
# NOT RUN
#########################

# options(warn = -1)
# 
# <- permute_means(all.df, "recomb_rate", "ecology")
# 
# options(warn = 0)
# 
# per1 <- replicate(1000, perm(all.df))
# 
# 	
# 
# all.df %>%
# 	lm(data = ., recomb_rate ~ euc.distance *ecology) %>%
# 	anova %>%
# 	tidy
# 
# all.df %>%
# 	lm(data = ., recomb_rate ~ least.cost.distance *ecology) %>%
# 	tidy
# 
# all.df %>%
# 	lm(data = ., recomb_rate ~ fst *ecology) %>%
# 	tidy
