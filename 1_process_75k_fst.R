####Match EVs to 75k stats files, fit models, write model results to file

rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("IRanges")
library("ggplot2")
library("dplyr")
library("robustbase")
library("devtools")
library("wesanderson")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Input file locations
################################################################################

#ev dir location and file list
ev.dir <- file.path(getwd(),"evs")
ev.files <- list.files(ev.dir, pattern = "txt",full.names = TRUE)

#stats files location and list
stats.dir <- file.path(getwd(),"stats/75k_all")
stats.files <- list.files(stats.dir,"*sliding*", full.names = TRUE)

################################################################################
# Add evs, call outliers, and  fit linear models
################################################################################

# how to filter the data & call outliers
filt_outliers_function <- function (matched.all){
	
	matched.all <- matched.all %>%
		filter(!is.na(fst)) %>%
		filter(!is.infinite(fst))
	
	matched.all$fst[matched.all$fst < 0] <- 0 
	matched.all$fst[matched.all$dxy < 0] <- 0 
	
	matched.all <- matched.all %>%
		mutate(fst.outlier = is.outlier(fst))%>%
		mutate(fst.outlier = is.outlier(dxy))%>%
		mutate(hs = (hexp1+hexp2)/2)
	
	return(matched.all)
}


# the type of linear model to fit
linear_model_function <- function (matched.all){
	
	model <- glm(as.numeric(fst.outlier) ~ recomb_rate + ds + gene_count,
		na.action = "na.omit",
		family = binomial,
		data = matched.all)
		
	return(model)
}


# apply the above function to the list of stats.files
coeff.df <- lapply(stats.files[1], match_evs, linear_model_function = linear_model_function, filt_outliers_function = filt_outliers_function)

#bind into a data frame
coeff.dat <- do.call("rbind", coeff.df)

#remove nas
coeff.dat <- coeff.dat[!is.na(coeff.dat$recomb_rate),]
coeff.dat <- coeff.dat %>%
  mutate(group = paste0(geography,"_",ecology))

################################################################################
# Add in region data (for relaxed geograpic classification)
################################################################################

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

#make new (relaxed) groups 
coeff.dat$group2 <- paste0(coeff.dat$geography2,"_",coeff.dat$ecology)

################################################################################
# Write to file
################################################################################

# write to file
write.table(coeff.dat, file = "analysis_ready/75k_stats_model_fst_fits.txt", row.names = FALSE, quote = FALSE)


