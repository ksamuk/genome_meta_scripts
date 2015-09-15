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

#########################
# INPUT FILES
#########################

fst.df <- read.table("meta_data/fst_df.txt", header = TRUE, stringsAsFactors = FALSE)
dist.df <- read.table("meta_data/pop_geo_distances.txt", header = TRUE, stringsAsFactors = FALSE)

#########################
# BODY
#########################

fst.df <- cbind(fst.df, dist.df[,7:10])

fst.df %>%
	ggplot(aes(x = log(euc.distance), y = fst, label = paste0(pop1,pop2))) +
	geom_point() + 
	#geom_smooth(method = "lm")+
	#facet_wrap(~group2)
