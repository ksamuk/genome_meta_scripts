# figure s3
rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("ggplot2")
library("dplyr")
library("wesanderson")
library("Hmisc")
library("ggthemes")
library("gridExtra")
library("cowplot")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible
#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

################################################################################
# load raw data
################################################################################

# load raw data
coeff.dat <- initialize_coeff_dat_files()

################################################################################
# create plots of recombination coefficients
################################################################################

# recombination rate vs fst -- relaxed groupings
group <- "group2.new"
stat <- "recomb_rate_fst"

fst_relaxed <- coeff.dat %>%
	filter(n_windows_fst > 100) %>%
	filter(!is.na(recomb_rate_fst)) %>%
	plot_dot_line_plot(., group, stat, label = "", 
								pal = pal, y_lab = "", theme_all = NULL, 
								point_size = 1, line_size = 2)

# recombination rate vs dxy -- relaxed groupings
group <- "group2.new"
stat <- "recomb_rate_dxy"

dxy_relaxed <- coeff.dat %>%
	filter(!is.na(intercept_dxy)) %>%
	plot_dot_line_plot(., group, stat, label = "", 
																	pal = pal, y_lab = "", theme_all = NULL, 
																	point_size = 1, line_size = 2)

# recombination rate vs fst -- strict groupings
group <- "group.new"
stat <- "recomb_rate_fst"
fst_strict <- coeff.dat %>%
	filter(!is.na(recomb_rate_fst)) %>%
	plot_dot_line_plot(., group, stat, label = "", 
																 pal = pal, y_lab = "", theme_all = NULL, 
																 point_size = 1, line_size = 2)

# recombination rate vs fst -- strict groupings
group <- "group.new"
stat <- "recomb_rate_dxy"

dxy_strict <- coeff.dat %>%
	filter(!is.na(recomb_rate_dxy)) %>%
	plot_dot_line_plot(., group, stat, label = "", 
																 pal = pal, y_lab = "", theme_all = NULL, 
																 point_size = 1, line_size = 2)

################################################################################
# create unified plot
################################################################################

# plot both graphs to pdf

pdf(file = "figures/figureS3.pdf", height = 8.5, width = 8.5)

labels <- c("fst_relaxed", "fst_strict",
						"dxy_strict","dxy_relaxed")

plot_grid(fst_relaxed, fst_strict, 
					dxy_relaxed, dxy_strict,
					ncol = 2, labels = labels, align = "hv")
