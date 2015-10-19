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

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

#########################
# INPUT FILES
#########################

# df containing mean fst values for all populations
fst.df <- read.table("meta_data/fst_df.txt", header = TRUE, stringsAsFactors = FALSE)

# df containing linear model coefficients for outliers vs. recomb rate
fst.model.fits <- read.table("analysis_ready/75k_stats_model_fits_fst.txt", header = TRUE, stringsAsFactors = FALSE)

# df containing distances (least cost, great circle) between populations
# also measures of isolation from ocean (dist to coast, i.e. how inland, for each population)
dist.df <- read.table("meta_data/pop_geo_distances.txt", header = TRUE, stringsAsFactors = FALSE)

#########################
# PRE PROCESSING
#########################

# bind the three data files together
all.df <- cbind(fst.model.fits, dist.df[,7:10], fst.df$fst)
names(all.df)[length(all.df)] <- "fst"

# create a 
fst.iso.resid <- all.df %>%
	mutate(isolation.sum = dist.to.coast1 + dist.to.coast2) %>%
	lm(fst ~ isolation.sum, data = .) %>%
	residuals

all.df$fst.resid <- fst.iso.resid 

all.df <- all.df %>% 
	mutate(isolation.sum = dist.to.coast1 + dist.to.coast2)%>%
	mutate(comparison = paste0(pop1, pop2))%>%
	rowwise() %>%
	mutate(isolation.max = max(dist.to.coast1, dist.to.coast2)) 

#########################
# Figure S7
#########################

pdf(file = "figures/figureS7.pdf", height = 8.5, width = 8.5)
	
size_all <- 2
alpha_all <- 0.2

theme_all <- theme_hc(base_size = 16) + 
	theme(legend.position = "none",
				strip.background = element_blank(),
				strip.text.x = element_blank(),
				axis.text.x = element_text(angle = 45, hjust = 1)) 

isolation_sum_plot <- all.df %>%
	ggplot(aes(x = isolation.sum, y = recomb_rate, color = ecology, label = comparison)) +
	geom_point(size = size_all, alpha = alpha_all)+
	scale_color_manual(values = pal[3:4])+
	geom_smooth(method = "lm", size = size_all)+
	facet_wrap(~ ecology)+
	theme_all+
	xlab("Non-oceanic distance (km)")+
	ylab("Recombination coeffficient")
	

gc_distance_plot <- all.df %>%
	ggplot(aes(x = euc.distance, y = recomb_rate, color = ecology, label = comparison))+
	geom_point(size = size_all, alpha = alpha_all)+
	scale_color_manual(values = pal[3:4])+
	geom_smooth(method = "lm", size = size_all)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab("Great circle distance (km)")+
	ylab("Recombination coeffficient")

lc_distance_plot <- all.df %>%
	ggplot(aes(x = least.cost.distance, y = recomb_rate, color = ecology, label = comparison)) +
	geom_point(size = size_all, alpha = alpha_all)+
	scale_color_manual(values = pal[3:4])+
	geom_smooth(method = "lm", size = size_all)+
	theme_all+
	facet_wrap(~ ecology)+
	xlab("Least cost distance (km)")+
	ylab("Recombination coeffficient")

fst_plot <- all.df %>%
	ggplot(aes(x = fst, y = recomb_rate, color = ecology, label = comparison)) +
	geom_point(size = size_all, alpha = alpha_all)+
	#geom_text(size = 3)+
	scale_color_manual(values = pal[3:4])+
	geom_smooth(method = "lm", size = size_all)+
	theme_all +
	facet_wrap(~ ecology)+
	xlab("Average FST")+
	ylab("Recombination coeffficient")
	
plot_grid(isolation_sum_plot, gc_distance_plot, lc_distance_plot, fst_plot)

dev.off()

#########################
# Permutations re: above
#########################

shuffle_all_df <- function(all.df){
	all.df$ecology <- sample(all.df$ecology, size = length(all.df$ecology))
	all.df
}

perm <- function(all.df) {
	shuffle_all_df(all.df) %>% 
	lm(data = ., recomb_rate ~ isolation.sum *ecology) %>% 
		tidy %>% 
		.[4,2] %>%
	  return
}

per1 <- replicate(1000, perm(all.df))

	

all.df %>%
	lm(data = ., recomb_rate ~ euc.distance *ecology) %>%
	anova %>%
	tidy

all.df %>%
	lm(data = ., recomb_rate ~ least.cost.distance *ecology) %>%
	tidy

all.df %>%
	lm(data = ., recomb_rate ~ fst *ecology) %>%
	tidy
