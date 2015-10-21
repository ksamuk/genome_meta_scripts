# plot figure 2

rm(list=ls())

################################################################################
# Libraries and shared functions
################################################################################

library("ggplot2")
library("dplyr")
library("wesanderson")
library("grid")
library("gridExtra")

list.files("shared_functions", full.names = TRUE) %>% sapply(source)

select <- dplyr::select
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

################################################################################
# Read in and format input files
################################################################################

# the clustering file
cluster.df <- initialize_clustering_output()

# filter out linakge groups with super low data
cluster.df <- cluster.df %>%
	filter(n.sites > 30) %>%
	filter(num.outliers > 5)

################################################################################
# Figure 2 (clustering data, relaxed)
################################################################################

# theme settings
point_size <- 1
line_size <- 3
size <- 16
theme.all <- theme(legend.position="none", 
                   axis.title.x = element_blank(), 
                   axis.title.y = element_text(vjust=1.5),
                   axis.text.x = element_blank(),
                   axis.line.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text = element_blank(), 
                   strip.background = element_blank(),
                   legend.title=element_blank())

# counts of clustered lgs
prop.clustered <- cluster.df %>%
  mutate(more.clustered = nnd.emp.percentile < 0.05) %>%
  mutate(less.clustered = nnd.emp.percentile > 0.95) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  select(comparison, lg, group, group2, more.clustered, less.clustered, group2.new) %>%
  group_by(group2.new, comparison) %>%
  summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
  mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
  filter(num.lg > 1) %>%
		plot_dot_line_plot(data = ., group = "group2.new", 
											 stat = "more.clustered", label = "", 
											 pal = pal, y_lab = "Prop. chromosomes clustered",
											 theme_all = theme.all, point_size = point_size, line_size = line_size)


# coeff of disp
coeff.dispersion <- cluster.df %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2.new, comparison) %>%
  summarise(disp.out = mean(disp.out)) %>%
	plot_dot_line_plot(data = ., group = "group2.new", 
										 stat = "disp.out", label = "", 
										 pal = pal, y_lab = "Outlier dispersion coefficient",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)

# nnd.diff

nnd.diff <- cluster.df %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2.new, comparison) %>%
  summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	filter(!is.na(nnd.diff.sd)) %>%
	filter(nnd.diff.sd > -4)%>% # filters out one particular impossible value
	plot_dot_line_plot(data = ., group = "group2.new", 
										 stat = "nnd.diff.sd", label = "", 
										 pal = pal, y_lab = "Expected NND - Outlier NND (cM)",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)

### representative chromosome plot

# read in rep snp files
if (!file.exists("meta_data/representative_nnd_df.txt")){
  
  rep.files <- list.files("stats/snp_representative", full.names = TRUE)
  rep.df <- lapply(rep.files, extract_lg)
  # calculate empirical nnd 
  rep.nnd.df <- lapply(rep.df, calc_emp_nnd_dist)
  rep.nnd.df <- bind_rows(rep.nnd.df)
  rep.nnd.df$group2 <- ifelse(rep.nnd.df$comparison == "cr_stream.wc_stream", "para_S", rep.nnd.df$group)
  write.table(rep.nnd.df, "meta_data/representative_nnd_df.txt", quote = FALSE, row.names = FALSE)
  
}else{
  rep.nnd.df <- read.table("meta_data/representative_nnd_df.txt", header = TRUE, stringsAsFactors = FALSE)
}

nnd.rep <- rep.nnd.df %>%
  ggplot(aes(x = group2, y = gen.pos, color = group2, fill = group2)) +
  geom_jitter(position = position_jitter(height = 2), size = point_size)+
  #geom_text(mapping = NULL, label = "D", x = 4.25, y= 78, color = "black" ,alpha = 0.25)+
  scale_color_manual(values = pal) +   
  scale_fill_manual(values = pal) +
  theme_classic(base_size = size) +
  theme.all+
  coord_cartesian(ylim = c(0,80))+
  ylab("Map position on LG 4 (cM)")


# plot all the things
dev.off()
pdf(file = "figures/Figure2.pdf", height = 8.5, width = 8.5, onefile = FALSE)
plot_grid(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep, align = "hv")
dev.off()

################################################################################
# Figure S2 (clustering data, strict)
################################################################################

# counts of clustered lgs
prop.clustered <- cluster.df %>%
	mutate(more.clustered = nnd.emp.percentile < 0.05) %>%
	mutate(less.clustered = nnd.emp.percentile > 0.95) %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	select(comparison, lg, group, group2, more.clustered, less.clustered, group.new) %>%
	group_by(group.new, comparison) %>%
	summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
	mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
	filter(num.lg > 1) %>%
	plot_dot_line_plot(data = ., group = "group.new", 
										 stat = "more.clustered", label = "", 
										 pal = pal, y_lab = "Prop. chromosomes clustered",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)


# coeff of disp
coeff.dispersion <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group.new, comparison) %>%
	summarise(disp.out = mean(disp.out)) %>%
	plot_dot_line_plot(data = ., group = "group.new", 
										 stat = "disp.out", label = "", 
										 pal = pal, y_lab = "Outlier dispersion coefficient",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)

# nnd.diff

nnd.diff <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group.new, comparison) %>%
	summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	filter(!is.na(nnd.diff.sd)) %>%
	filter(nnd.diff.sd > -4)%>% # filters out one particular impossible value
	plot_dot_line_plot(data = ., group = "group.new", 
										 stat = "nnd.diff.sd", label = "", 
										 pal = pal, y_lab = "Expected NND - Outlier NND (cM)",
										 theme_all = theme.all, point_size = point_size, line_size = line_size)


# plot all the things
dev.off()
pdf(file = "figures/FigureS2.pdf", height = 8.5, width = 8.5, onefile = FALSE)
plot_grid(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep, align = "hv")
dev.off()
