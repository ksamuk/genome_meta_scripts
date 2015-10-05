# figure s3
rm(list=ls())

library("ggplot2")
library("dplyr")
library("wesanderson")
library("Hmisc")
library("ggthemes")
library("gridExtra")

# initialize graphic device
dev.off()
dev.new()
quartz(width = 6, height = 6)

# load raw data

coeff.dat.fst <- read.table(file = "analysis_ready/75k_stats_model_fst_fits.txt", header = TRUE, stringsAsFactors = FALSE)

# rename grouping variables for plotting
group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
group.rename <- c("Allopatry\n Divergent", "Allopatry\n Parallel", "Gene Flow\n Divergent", "Gene Flow\n Parallel")
coeff.dat.fst$group2.new <-group.rename[match(coeff.dat.fst$group2, group.old.names)]
coeff.dat.fst$group.new <- group.rename[match(coeff.dat.fst$group, group.old.names)]

# palatte town
pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]

# recombination rate coefficients -- relaxed groupings
relaxed <- coeff.dat.fst %>%
		ggplot(aes(x = group2.new, y = recomb_rate))+
			geom_point(position = position_jitter(width = 0.3), size = 1)+
			geom_errorbar(stat = "hline", 
										yintercept = "mean",
										width = 0.8, size = 2, 
										aes(ymax = ..y.. , ymin = ..y..),
										color = c(pal[3], pal[2], pal[4], pal[1]))+
				geom_text(mapping = NULL, label = "A", x = 4.5, y= 0.15, color = "black" ,alpha = 0.25, size = 4)+
				theme_hc(base_size = 16) +
				theme(axis.title.y = element_text(vjust=4, hjust = -3),
							axis.text.x = element_blank(),
							axis.ticks.x = element_blank(),
							plot.margin = unit(c(1,1,1,1),"cm"))+
				ylab("Recombination rate coefficient")+
				xlab(NULL)

# recombination rate coefficients -- strict groupings

strict <- coeff.dat.fst %>%
	ggplot(aes(x = group.new, y = recomb_rate))+
	geom_point(position = position_jitter(width = 0.3), size = 1)+
	geom_errorbar(stat = "hline", 
								yintercept = "mean",
								width = 0.8, size = 2, 
								aes(ymax = ..y.. , ymin = ..y..),
								color = c(pal[3], pal[2], pal[4], pal[1]))+
	geom_text(mapping = NULL, label = "B", x = 4.5, y= 0.15, color = "black" ,alpha = 0.25, size = 4)+
	theme_hc(base_size = 16) +
	theme(axis.title.y = element_text(vjust=4, hjust = -1.25),
				axis.title.x = element_text(vjust=-3),
				axis.ticks.x = element_blank(),
				plot.margin = unit(c(0,1,1,1),"cm"))+
	ylab("Recombination Rate Coefficient")+
	xlab(NULL)

# plot both graphs to pdf
grid.arrange(relaxed, strict)
quartz.save("figures/figureS3.pdf", type = "pdf")



