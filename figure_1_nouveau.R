# figure 1 
rm(list=ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("ggplot2")
library("readr")
library("dplyr")
library("wesanderson")
library("Hmisc")
library("ggthemes")
library("gridExtra")
library("cowplot")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

################################################################################
# plot color scheme/theme
################################################################################

#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

theme_all <- theme_hc(base_size = 16) +
	theme(axis.title.y = element_blank(),
				axis.text.x = element_blank(),
				axis.ticks.x = element_blank(),
				plot.margin = unit(c(0,0,0,0),"cm"))

################################################################################
# load raw data
################################################################################

coeff.dat <- initialize_coeff_dat_files()

stats_df <- read_delim(file = "analysis_ready/75k_stats_combined.txt", delim = " ")
stats_df <- add_region_data(stats_df)

stats_df <- stats_df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = "."))

rep.comparisons <- c("cp.marine.sk.marine", 
										 "constance.stream.joes.stream", 
										 "mariager.marine.joes.lake", 
										 "pri.limnetic.pri.benthic")

################################################################################
# Figure 1 A: Averaged Fitted model coefficients 
################################################################################

# create functions based on the averaged the regression coefficients for each selection/gene flow group
avg_regression_functions <- coeff.dat %>% filter(n_windows_fst > 200) %>% create_average_regression_functions_fst
#avg_regression_functions <- coeff.dat %>% filter(n_windows_dxy > 200) %>% create_average_regression_functions_dxy

# plot the functions 

# an empty plot
avg_reg_plot <- ggplot()+
	scale_x_continuous(limits=c(0, 25), expand=c(0,0))+
	scale_y_continuous(limits=c(0, 0.2), expand = c(0,0))

# add in the avg reg function plots 
for(i in 1:length(avg_regression_functions)){
	avg_reg_plot <- avg_reg_plot + 
		stat_function(aes(y = 0), fun = avg_regression_functions[[i]], colour = pal[i], size = 3)
}
	
avg_reg_plot 

################################################################################
# Figure 1 B: Jitter plot of model coefficients 
################################################################################

fst_relaxed <- coeff.dat %>% 
	filter(!is.na(recomb_rate_fst)) %>%
	filter(n_windows_fst > 100) %>%
	plot_dot_line_plot(., group = "group2", stat = "recomb_rate_fst", label = "", 
										 pal = pal, y_lab = "", theme_all = NULL, 
										 point_size = 1, line_size = 2)

dxy_relaxed <- coeff.dat %>% 
	filter(n_windows_dxy > 100) %>% 
	filter(!is.na(recomb_rate_dxy)) %>%
	plot_dot_line_plot(., group = "group2.new", stat = "recomb_rate_dxy", label = "", 
										 pal = pal, y_lab = "", theme_all = NULL, 
										 point_size = 1, line_size = 2)

hs_relaxed <- coeff.dat %>% 
	filter(!is.na(recomb_rate_hs)) %>%
	filter(n_windows_fst > 100) %>%
	plot_dot_line_plot(., group = "group2", stat = "recomb_rate_hs", label = "", 
										 pal = pal, y_lab = "", theme_all = NULL, 
										 point_size = 1, line_size = 2)

################################################################################
# Figure 1 C: Representitive FST distributions
################################################################################
