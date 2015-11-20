##### FIGURES 1 & S1
##### This script generates the two above figures. Note that the Quartz graphic device requires a Mac OS X.

# figure 1
rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("ggplot2")
library("dplyr")
#library("devtools")
#install_github("karthik/wesanderson")
library("lazyeval")
library("IRanges")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Raw data
################################################################################

#ev dir location and file list
ev.dir <- file.path(getwd(),"evs")
ev.files <- list.files(ev.dir, pattern = "txt",full.names = TRUE)

#stats files location and list
stats.dir <- file.path(getwd(),"stats/75k_all")
stats.files <- list.files(stats.dir,"*sliding*", full.names = TRUE)

coeff.df <- mclapply(stats.files, bind_stats_files, 
										 mc.cores = 6, mc.silent = FALSE, mc.preschedule = FALSE)

out.df <- mclapply(coeff.df, filter_data_call_outliers_function, 
										 mc.cores = 6, mc.silent = FALSE, mc.preschedule = FALSE)

coeff.dat <- do.call("rbind", out.df)

coeff.dat <- coeff.dat[!is.na(coeff.dat$recomb_rate),]
coeff.dat <- coeff.dat %>%
	mutate(group = paste0(geography,"_",ecology))

################################################################################
# Add in region data (for relaxed geograpic classification)
################################################################################

coeff.dat <- add_region_data(coeff.dat)

################################################################################
# Write to file
################################################################################

# write to file
dir.create("analysis_ready") %>% invisible
file.name <- "analysis_ready/75k_stats_combined.txt"
write.table(coeff.dat, file = file.path(file.name), row.names = FALSE, quote = FALSE)
