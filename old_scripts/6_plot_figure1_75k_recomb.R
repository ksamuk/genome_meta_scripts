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
library("wesanderson")
library("lazyeval")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Global figure options
################################################################################

# define color palatte for plotting
pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]

# define group order (for matching to palatte)
groups <- c("para_S","allo_S", "allo_D","para_D")

# read in linear model output
coeff.dat <- initialize_coeff_dat_files()


################################################################################
# Figure 1: FST & FST/DXY, Relaxed Groupings
################################################################################

# relaxed groups
group_variable <- c("group","group2")

# initialize graphics device
dev.off()
pdf(file = "figures/Figure1.pdf", height = 8.5, width = 12, onefile = FALSE)

# plot : fst outliers, relaxed

label.scaling <-  1.5
line.weight <-  8
layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(6, 6))
par(mgp=c(2.5,1,0), mar=c(5,4,2,1), mex = 1, bty = "l")

plot_averaged_regression_lines(coeff.dat, groups, pal, group_variable[2],
															 "fst", adjust = 0.0005, label = "A",
															 ylim = c(0,0.16), xlim = c(0,50), 
															 ylab = c(expression('F'["ST"]*" Outlier Probability")))
plot_averaged_regression_lines(coeff.dat, groups, pal, group_variable[2], 
															 "dxy", adjust = 0.0001, label ="B",
															 ylim = c(0,0.1), xlim = c(0,50), 
															 ylab = c(expression('Joint F'["ST"]*'/D'["XY"]*" Outlier Probability")))	

legend(20,0.023, # places a legend at the appropriate place 
			 c("Gene Flow Parallel","Allopatry Parallel","Allopatry Divergent","Gene Flow Divergent"), 
			 lty = c(1,1), # gives the legend appropriate symbols (lines)
			 lwd = line.weight, col = pal, box.lty = 0, cex = label.scaling, y.intersp = 1,xjust = 0, seg.len = 0.75)

dev.off()

################################################################################
# Figure S4: FST & FST/DXY, Strict Groupings
################################################################################

# relaxed groups
group_variable <- c("group","group2")

# initialize graphics device
dev.off()
pdf(file = "figures/FigureS4.pdf", height = 8.5, width = 12, onefile = FALSE)

# plot : fst outliers, relaxed

label.scaling = 1.5
line.weight = 8
layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(6, 6))
par(mgp=c(2.5,1,0), mar=c(5,4,2,1), mex = 1, bty = "l")

plot_averaged_regression_lines(coeff.dat, groups, pal, group_variable[1],
															 "fst", adjust = 0.0005, label = "A",
															 ylim = c(0,0.16), xlim = c(0,50), 
															 ylab = c(expression('F'["ST"]*" Outlier Probability")))
plot_averaged_regression_lines(coeff.dat, groups, pal, group_variable[1], 
															 "fst_dxy", adjust = 0.0001, label ="B",
															 ylim = c(0,0.025), xlim = c(0,50), 
															 ylab = c(expression('Joint F'["ST"]*'/D'["XY"]*" Outlier Probability")))	

legend(20,0.025, # places a legend at the appropriate place 
			 c("Gene Flow Parallel","Allopatry Parallel","Allopatry Divergent","Gene Flow Divergent"), 
			 lty = c(1,1), # gives the legend appropriate symbols (lines)
			 lwd = line.weight, col = pal, box.lty = 0, cex = label.scaling, y.intersp = 1,xjust = 0, seg.len = 0.75)

dev.off()

################################################################################
# NOT RUN: DXY only
################################################################################

plot_averaged_regression_lines(coeff.dat, groups, pal, group_variable[1], 
															 "dxy",
															 ylim = c(0,0.1), xlim = c(0,50), 
															 ylab = c(expression('D'["XY"]*" Outlier Probability")))
coeff.dat %>%
ggplot(aes(x = group2, y = recomb_rate_dxy)) +
	geom_point()

coeff.dat %>%
	filter(group2 == "para_S") %>%
	ggplot(aes(x = group2, y = recomb_rate_dxy)) +
	geom_point()

	summarise(intercept_dxy_mean = mean(intercept_dxy, na.rm = TRUE), slope_dxy_mean = mean(recomb_rate_dxy, na.rm = TRUE))


coeff.dat %>% 
	summarise(intercept_fst_mean = mean(intercept_fst, na.rm = TRUE), slope_fst_mean = mean(recomb_rate_fst, na.rm = TRUE))


