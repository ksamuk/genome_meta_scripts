############################################################
# inspect fst vs. geographic distance relationship
# KS Aug 2015
############################################################

rm(list =ls())

#########################
# LIBRARIES & FUNCTIONS
#########################

library("dplyr")
library("ggplot2")
library("visreg")

#########################
# INPUT FILES
#########################

fst.df <- read.table("meta_data/fst_df.txt", header = TRUE, stringsAsFactors = FALSE)
fst.model.fits <- read.table("analysis_ready/75k_stats_model_fst_fits.txt", header = TRUE, stringsAsFactors = FALSE)
dist.df <- read.table("meta_data/pop_geo_distances.txt", header = TRUE, stringsAsFactors = FALSE)


#########################
# BODY
#########################

all.df <- cbind(fst.model.fits, dist.df[,7:10], fst.df$fst)
names(all.df)[length(all.df)] <- "fst"

all.df %>%
	mutate(isolation = dist.to.coast1) %>%
	ggplot(aes(x = log(isolation), y = recomb_rate, color = group2)) +
	geom_point(size = 3)

recomb.lm <- all.df %>%
	mutate(isolation = dist.to.coast1 + dist.to.coast2) %>%
	lm(recomb_rate ~ fst + euc.distance + isolation + ecology, data = .)

