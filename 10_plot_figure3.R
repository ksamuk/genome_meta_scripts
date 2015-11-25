##### FIGURE 3

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
library("readr")
library("hexbin")
library("viridis")
library("ggthemes")
library("scales")

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Raw data
################################################################################

# combined stats file

stats_df <- read_delim(file = "analysis_ready/75k_stats_combined.txt", delim = " ")
stats_df <- add_region_data(stats_df)

stats_df <- stats_df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = "."))

rep.comparisons <- c("cp.marine.sk.marine", 
										 "constance.stream.joes.stream", 
										 "mariager.marine.joes.lake", 
										 "pri.limnetic.pri.benthic")

################################################################################
# Main plot: FST
################################################################################

#lazer beam plot

stats_df %>%
	filter(comparison %in% rep.comparisons) %>%
	filter(var.sites >= 3) %>%
	group_by(comparison) %>%
	filter(lg %in% c(4)) %>%
	mutate(fst_outlier_value = ifelse(fst.outlier == TRUE, fst, NA)) %>%
	mutate(fst_non_outlier_value = ifelse(fst.outlier == FALSE, fst, NA)) %>%
	ungroup %>%
	filter(recomb_rate <= 30 ) %>%
	ggplot(aes(x = midpos/1000000, y = fst))+
	geom_point(aes(x = midpos/1000000, y = fst_outlier_value, color = "zzzzred"))+
	geom_point(aes(x = midpos/1000000, y = fst_non_outlier_value, color = "zzzzgrey"))+
	#geom_line(color = "grey")+
	#geom_point(aes(x = midpos/1000000, y = fst_outlier), shape = "|", size = 4)+
	stat_smooth(method = "loess", aes(color = group2), size = 1.5, se = FALSE)+
	stat_smooth(method = "loess", aes(color = "zzz", x = midpos/1000000, y = recomb_rate/10), size = 1.5, se = FALSE, linetype ="11111111")+
	theme_hc_border()+
	ylab(expression('F'["ST"]))+
	xlab("Chromosomal position (MB)")+
	ylim(0, 1)+
	theme(panel.grid = element_blank(),
				strip.background = element_blank(),
				legend.position = "none",
				panel.border = element_rect(color="grey", fill=NA))+
	facet_grid(group2~lg)+
	scale_color_manual(values = c(pal, 1, "grey" ,1))+
	scale_x_continuous(labels = comma)
ggsave(filename = "figures/Figure3.pdf", height = 8.5, width = 4.25)

################################################################################
# Get recomb_rate axes: FST
################################################################################

#lazer beam plot
stats_df %>%
	group_by(comparison) %>%
	filter(lg  == 4) %>%
	mutate(fst_outlier_value = ifelse(fst.outlier == TRUE, fst, NA)) %>%
	mutate(fst_non_outlier_value = ifelse(fst.outlier == FALSE, fst, NA)) %>%
	ungroup %>%
	ggplot(aes(x = midpos/1000000, y = recomb_rate))+
	stat_smooth(method = "loess", aes(color = "zzz", x = midpos/1000000, y = recomb_rate), size = 1.5, se = FALSE, linetype = 3)+
	theme_hc()+
	ylim(0, 10)+
	ylab("Recombination rate cM/MB")+
	xlab("Chromosomal position (MB)")+
	theme(panel.grid = element_blank(),
				strip.background = element_blank(),
				legend.position = "none")+
	facet_grid(group2~lg)+
	scale_color_manual(values = c(pal,1))+
	scale_x_continuous()
ggsave(filename = "figures/Figure3_rep_axis_hack.pdf", height = 8.5, width = 4.5)

################################################################################
# DXY
################################################################################

stats_df %>%
	filter(comparison %in% rep.comparisons) %>%
	filter(sites > 4) %>%
	filter(var.sites > 4) %>%
	group_by(comparison) %>%
	filter(lg %in% c(4,7)) %>%
	#mutate(dxy_outlier = ifelse(dxy.outlier == TRUE, 0.02, NA)) %>%
	ungroup %>%
	filter(recomb_rate <= 30) %>%
	ggplot(aes(x = midpos/1000000, y = dxy/sites))+
	geom_line(color = "grey")+
	#geom_point(aes(x = midpos/1000000, y = dxy_outlier), shape = "|", size = 4)+
	stat_smooth(method = "loess", aes(color = group2), size = 1.5, se = FALSE)+
	#stat_smooth(method = "loess", aes(color = "zzgreen", x = midpos/1000000, y = abs(hexp1-hexp2)), size = 1.5, se = FALSE, linetype = 1)+
	stat_smooth(method = "loess", aes(color = "zzz", x = midpos/1000000, y = recomb_rate/100000), size = 1.5, se = FALSE, linetype ="11111111")+
	theme_hc_border()+
	#ylim(0, 0.02)+
	ylab(expression('D'["XY"]))+
	xlab("Chromosomal position (MB)")+
	theme(panel.grid = element_blank(),
				strip.background = element_blank(),
				legend.position = "none",
				panel.border = element_rect(color="grey", fill=NA))+
	facet_grid(group2~lg)+
	scale_color_manual(values = c(pal,1))+
	scale_x_continuous(labels = comma)



stats_df %>%
	group_by(comparison) %>%
	#filter(sites > 20) %>%
	filter(recomb_rate <= 25 )%>%
	filter(var.sites >= 3) %>%
	mutate(low_site_window = sites < 500) %>%
	#filter(dxy < 0.05)%>%
	#sample_frac(0.25) %>%
	#filter(lg %in% c(1,4,7)) %>%
	ggplot(aes(x = sites))+
	#ggplot(aes(x = recomb_rate , y = as.numeric(dxy_avg_outlier)))+
  #ggplot(aes(y = as.numeric(low_site_window), x = recomb_rate))+
	#ggplot(aes(x = sites, y = var.sites))+
	#ggplot(aes(x = recomb_rate, y = dxy/(var.sites/sites)))+
	#geom_point(alpha = 0.01)+
	geom_histogram(binwidth = 100)+
	xlim(0,6000)+
  #geom_smooth()+
	facet_wrap(~group2, scales = "free_y")
	
	stats_df %>%
		filter(var.sites > 10) %>%
		filter(hs == 0) %>%
		View()
		ggplot(aes(x = dxy))