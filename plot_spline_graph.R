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
library("readr")
library("hexbin")
library("viridis")
library("ggthemes")

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

################################################################################
# Raw data
################################################################################

# combined stats file

stats_df <- read_delim(file = "analysis_ready/75k_stats_combined.txt", delim = " ")

stats_df <- stats_df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = "."))

stats_df_filt <- stats_df %>%
	filter(!is.na(fst)) %>%
	filter(!is.infinite(fst)) %>%
	filter(recomb_rate <= 25) %>% 
	filter(var.sites >= 2) %>%
	filter(ds <= 3) %>%
	filter(gene_count <= 15) %>%
	filter(lg!=19) %>%
	mutate(dxy_adj = ifelse(sites >= 500, dxy, NA)) %>% 
	group_by(comparison) %>%
	mutate(fst.outlier = is.outlier(fst))%>%
	mutate(dxy.outlier = is.outlier(dxy)) %>%
	mutate(both.outlier = fst.outlier & dxy.outlier)%>%
	mutate(hs = (hexp1+hexp2)/2) %>%
	ungroup

stats_df_filt %>%
	group_by(comparison) %>%
	sample_n(100, replace = TRUE) %>%
	#mutate(low_recomb = ifelse(recomb_rate <= 5, TRUE, FALSE)) %>%
	ggplot(aes(x = group2, y = hs, color = fst.outlier))+
	#ggplot(aes(x = recomb_rate, y = log(dxy_adj+1)))+
	geom_boxplot()+
	#geom_point()+
	#geom_smooth()+
	facet_wrap(~group2)

# fst "scaled" plot
stats_df %>%
	group_by(comparison) %>%
	mutate(fst_scale = scale(fst, center = FALSE)) %>%
	ungroup %>%
	mutate(comparison = reorder(comparison, factor(group2) %>% as.numeric)) %>%
	filter(recomb_rate <= 20 ) %>%
	sample_frac(0.2) %>%
	ggplot(aes(x = recomb_rate, y = fst))+
	#geom_hex(size = 20)+
	geom_density2d()+
	#geom_point(aes(color = group2), size = 2, alpha = 0.1)+
	#stat_smooth(method= "lm", color = "black", size = 2)+
	#theme_hc()+
	facet_wrap(~group2)+
	#scale_color_manual(values = pal)
ggsave(filename = "scaled_fst_plot.png", height = 8.5, width = 11)

stats_df_filt %>%
	mutate(comparison = reorder(comparison, factor(group2) %>% as.numeric)) %>%
	filter(recomb_rate <= 25 ) %>%
	group_by(comparison) %>%
	sample_n(200, replace = TRUE) %>%
	ggplot(aes(x = recomb_rate, y = fst.outlier %>% as.numeric %>% jitter(0.5), color = hs)) +
	geom_point(alpha = 0.5, size = 2)+
	geom_smooth(method = "loess")+
	facet_wrap(~group2)


stat_function(aes(y = 0), fun = avg_regression_functions[[i]], colour = pal[i], size = 3)

stat_density2d(aes(alpha=..level.., fill=..level..), size=2, 
							 bins = 30, geom="polygon") + 
	scale_fill_gradient() +
	scale_alpha(range = c(0.00, 1), guide = FALSE) +
	#geom_density2d(colour="black", bins = 30)+




#lazer beam plot
rep.comparisons <- list.files("stats/snp_representative") %>% 
	gsub(".stats.txt.gz", "", .) %>% 
	gsub("\\.allo|\\.para|\\.S|\\.D", "", .) %>%
	gsub("_", "\\.", .)
rep.comparisons
stats_df %>%
	filter(comparison %in% rep.comparisons) %>%
	group_by(comparison) %>%
	filter(lg %in% c(1,4,7)) %>%
	mutate(fst_scale = scale(fst, center = FALSE)) %>%
	ungroup %>%
	filter(recomb_rate <= 30 ) %>%
	ggplot(aes(x = midpos/1000000, y = recomb_rate))+
	#ggplot(aes(x = midpos/1000000, y = fst))+
	stat_smooth(method = "loess", aes(color = "zzz", x = midpos/1000000, y = recomb_rate), size = 1.5, se = FALSE, linetype = 3)+
	
	#ggplot(aes(x = midpos, y = log(fst_scale+1)))+
	#geom_hex()+
	#geom_point(size = 1, alpha = 1)+
	#geom_line(color = "grey")+
	#stat_smooth(method = "loess", aes(color = group2), size = 1.5, se = FALSE)+
	stat_smooth(method = "loess", aes(color = "zzz", x = midpos/1000000, y = recomb_rate), size = 1.5, se = FALSE, linetype = 3)+
	theme_hc()+
	ylab(expression('F'["ST"]))+
	xlab("Chromosomal position (MB)")+
	theme(panel.grid = element_blank(),
				strip.background = element_blank(),
				legend.position = "none")+
	facet_grid(group2~lg)+
	scale_color_manual(values = c(pal,1))+
	scale_x_continuous(labels = comma)
ggsave(filename = "scaled_fst_plot.pdf", height = 4, width = 4)


get_smooth_values <- function(comp_df, stat = "fst", required_recomb_obs = 20){ 
	
	if (is.numeric(required_recomb_obs)) {
		comp_df <- comp_df %>%
		filter(recomb_rate < 50) %>%
			mutate(recomb_bin = cut(recomb_rate, breaks = 0:50)) 
		
		recomb_bin_counts <- comp_df %>%
			group_by(recomb_bin) %>%
			tally
		
		comp_df <- left_join(comp_df, recomb_bin_counts, by = "recomb_bin") %>%
			filter(n > required_recomb_obs)
	}
	
	
	plot <- comp_df %>%
		filter(recomb_rate < 50) %>%
		mutate(fst.outlier = as.numeric(fst.outlier)) %>%
		mutate(dxy.outlier = as.numeric(dxy.outlier)) %>%
	#ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier))) +
	ggplot(aes_string(x = "recomb_rate", y = stat)) +
		geom_point(alpha = 0.5)+
		stat_smooth(method = "loess", se = FALSE, n = 4)
		#stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, n = 5)
	
	smooth_dat <- ggplot_build(plot)$data[[2]]
	
	return(data.frame(comparison = unique(comp_df$comparison), 
										group = unique(comp_df$group),
										group2 = unique(comp_df$group2),
										smooth_x = smooth_dat$x, smooth_y = smooth_dat$y))
}
# fst
smooth_fst_df <- stats_df %>%
	group_by(comparison) %>%
	do(smooth_dat = get_smooth_values(.))

smooth_fst_df <- bind_rows(smooth_fst_df$smooth_dat)

# fst outliers
smooth_fst_out_df <- stats_df %>%
	group_by(comparison) %>%
	do(smooth_dat = get_smooth_values(., stat = "fst.outlier"))
	
smooth_fst_out_df <- bind_rows(smooth_fst_out_df$smooth_dat)

# dxy
smooth_dxy_df <- stats_df %>%
	group_by(comparison) %>%
	do(smooth_dat = get_smooth_values(., stat = "dxy"))

smooth_dxy_df <- bind_rows(smooth_dxy_df$smooth_dat)

# dxy outliers
smooth_dxy_out_df <- stats_df %>%
	group_by(comparison) %>%
	do(smooth_dat = get_smooth_values(., stat = "dxy.outlier"))

smooth_dxy_out_df <- bind_rows(smooth_dxy_out_df$smooth_dat)

# hs
smooth_hs_df <- stats_df %>%
	group_by(comparison) %>%
	do(smooth_dat = get_smooth_values(., stat = "hs"))

smooth_hs_df <- bind_rows(smooth_hs_df$smooth_dat)


# yo dawg i heard you like smooths 
# so i smoothed your smooth so you can smooth while you smooth
smooth_dxy_df %>%
	ggplot(aes(x = smooth_x, y = smooth_y, color = group2)) +
	geom_point()+
	#geom_smooth(size = 3)+
	theme_hc()+
	facet_wrap(~group2)+
	scale_color_manual(values = pal)

smooth_fst_out_df %>%
	#filter(!grepl("camp", comparison)) %>%
	ggplot(aes(x = smooth_x, y = smooth_y, color = group2)) +
	geom_smooth(size = 3)+
	theme_hc()+
	scale_color_manual(values = pal)

# 
smooth_fst_df %>%
	filter(!grepl("camp", comparison)) %>%
	ggplot(aes(x = smooth_x, y = smooth_y, color = group2)) +
	geom_smooth(method = "loess", size = 3)+
	theme_hc()+
	scale_color_manual(values = pal)


stats_df %>%
	filter(recomb_rate < 50) %>%
	group_by(fst.outlier) %>%
	sample_frac(0.1) %>%
	#ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier))) +
	ggplot(aes(x = recomb_rate, y = as.numeric(fst.outlier))) +
	geom_point(alpha = 0.01, size = 0.5)+
	stat_smooth(se = FALSE, size = 2, aes(color = group2)) +
	scale_color_manual(values = pal)

ggplot_build(p)$data[[3]]

stat_smooth(data = pax_dat, aes(x = recomb_rate, y = fst, outfit = fit <<-..y..))

outfit=fit<<-..y..



