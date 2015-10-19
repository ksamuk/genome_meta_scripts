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
pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
#pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

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

fst_relaxed <- plot_dot_line_plot(coeff.dat, group, stat, label = "A")

# recombination rate vs fst -- strict groupings
group <- "group.new"
stat <- "recomb_rate_fst"
fst_strict <- plot_dot_line_plot(coeff.dat, group, stat, label = "B")

# recombination rate vs fst -- relaxed groupings
group <- "group2.new"
stat <- "recomb_rate_dxy"

fst_dxy_relaxed <- plot_dot_line_plot(coeff.dat, group, stat, label = "C")

# recombination rate vs fst -- strict groupings
group <- "group.new"
stat <- "recomb_rate_dxy"

fst_dxy_strict <- plot_dot_line_plot(coeff.dat, group, stat, label = "D")

################################################################################
# create unified plot
################################################################################

# plot both graphs to pdf

pdf(file = "figures/figureS3.pdf", height = 8.5, width = 8.5)

plot_grid(fst_relaxed, fst_strict, fst_dxy_relaxed, fst_dxy_strict, 
					ncol = 1)

dev.off()


