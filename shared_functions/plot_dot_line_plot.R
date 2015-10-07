plot_dot_line_plot <- function(data, group, stat, label = ""){
	
	plot <- data %>%
		ggplot(aes_string(x = group, y = stat))+
		geom_point(position = position_jitter(width = 0.3), size = 1)+
		geom_errorbar(stat = "hline", 
									yintercept = "mean",
									width = 0.8, size = 2, 
									aes(ymax = ..y.. , ymin = ..y..),
									color = c(pal[3], pal[2], pal[4], pal[1]))+
		#geom_text(mapping = NULL, aes_string(label = label, hjust = 0, vjust = 0), color = "black" ,alpha = 0.25, size = 6)+
		theme_hc(base_size = 16) +
		theme(axis.title.y = element_blank(),
					axis.text.x = element_blank(),
					axis.ticks.x = element_blank(),
					plot.margin = unit(c(0,0,0,0),"cm"))+
		ylab(NULL)+
		xlab(NULL)
	
	return(plot)
	
}