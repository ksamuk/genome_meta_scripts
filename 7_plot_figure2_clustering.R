# plot figure 2

rm(list=ls())

# required libraries

library("ggplot2")
library("dplyr")
library("wesanderson")
library("grid")
library("gridExtra")

# some custom tweaks
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
select <- dplyr::select

# source in custom functions
list.files("shared_functions", full.names = TRUE) %>% sapply(source)

## the clustering file
cluster.df <- read.table(file = "analysis_ready/snp_clustering_metrics.txt", header = TRUE, stringsAsFactors = FALSE)
cluster.df <- cluster.df %>%
  mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp)

# dem plots?
group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
group.rename <- c("Allopatry Divergent", "Allopatry Parallel", "Gene Flow Divergent", "Gene Flow Parallel")
cluster.df$group2.new <-group.rename[match(cluster.df$group2, group.old.names)]
cluster.df$group.new <- group.rename[match(cluster.df$group, group.old.names)]

################# relaxed groupings

dev.off()
pdf(file = "figures/figure2.pdf", width = 8.5, height = 8.5, onefile=FALSE)
#quartz(width = 8.5, height = 8.5)

#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
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
  mutate(more.clustered = nnd.emp.percentile*2 < 0.05) %>%
  mutate(less.clustered = nnd.emp.percentile*2 > 0.95) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  select(comparison, lg, group, group2, more.clustered, less.clustered, group2.new) %>%
  group_by(group2.new, comparison) %>%
  summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
  mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
  filter(num.lg > 1) %>%
    ggplot(aes(x = group2.new, y = more.clustered, fill = group2.new, color = group2.new))+
    geom_jitter(size = 1)+
    geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, weight=1, color = 1, linetype = 1)+
    geom_text(mapping = NULL, label = "A", x = 4.25, y= 0.75, color = "black" ,alpha = 0.25)+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    #theme(panel.grid = element_blank())+
    theme_classic(base_size = size)+
    coord_cartesian(ylim=c(0,0.8))+
    theme.all+
    ylab("Prop. chromosomes clustered")

# coeff of disp
coeff.dispersion <- cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2.new, comparison) %>%
  summarise(disp.out = mean(disp.out))%>%
    ggplot(aes(x = group2.new, y = disp.out, fill = group2.new, color = group2.new))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    geom_jitter(size = 1)+
    geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, weight=1, color = 1, linetype = 1)+
    geom_text(mapping = NULL, label = "C", x = 4.25, y= 24, color = "black" ,alpha = 0.25)+
    theme_classic(base_size = size)+
    coord_cartesian(ylim=c(0,25))+
    theme.all+
    ylab("Outlier dispersion coefficient")

# nnd.diff

nnd.diff <- cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group2.new, comparison) %>%
  summarise(nnd.diff = mean(nnd.diff))%>%
  ggplot(aes(x = group2.new, y = nnd.diff, fill = group2.new, color = group2.new))+
  scale_color_manual(values = pal)+
  scale_fill_manual(values = pal) + 
  geom_jitter(size = 1)+
  geom_boxplot(notch = TRUE, outlier.size = 0, notchwidth = 0.75, width=1, color = 1, linetype = 1)+
  geom_text(mapping = NULL, label = "B", x = 4.25, y= 4.75, color = "black" ,alpha = 0.25)+
  theme_classic(base_size = size) +
  coord_cartesian(ylim=c(-5,5))+
  theme.all+
  ylab("Expected NND - Outlier NND (cM)")

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

rep.nnd.df$group2.new <-group.rename[match(rep.nnd.df$group2, group.old.names)]
rep.nnd.df$group.new <- group.rename[match(rep.nnd.df$group, group.old.names)]

nnd.rep <- rep.nnd.df %>%
  ggplot(aes(x = group2.new, y = gen.pos, color = group2.new, fill = group2.new)) +
  geom_jitter(position = position_jitter(height = 2), size = 1.5)+
  geom_text(mapping = NULL, label = "D", x = 4.25, y= 78, color = "black" ,alpha = 0.25)+
  scale_color_manual(values = pal) +   
  scale_fill_manual(values = pal) +
  theme_classic(base_size = size) +
  theme.all+
  coord_cartesian(ylim = c(0,80))+
  ylab("Map position on LG 4 (cM)")


# plot all the things
#grid.arrange(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep, ncol = 2)
grid_arrange_shared_legend(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep)
dev.off()
#quartz.save("figures/figure2.pdf", type = "pdf")


################# strict groupings

dev.off()
dev.new()
quartz(width = 8.5, height = 8.5)

#pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]
pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")
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
  mutate(more.clustered = nnd.emp.percentile*2 < 0.05) %>%
  mutate(less.clustered = nnd.emp.percentile*2 > 0.95) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  select(comparison, lg, group, group2, more.clustered, less.clustered, group.new) %>%
  group_by(group.new, comparison) %>%
  summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
  mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
  filter(num.lg > 1) %>%
    ggplot(aes(x = group.new, y = more.clustered, fill = group.new, color = group.new))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    geom_jitter(size = 1)+
    geom_boxplot(outlier.size = 0, weight=1, color = 1, linetype = 1)+
    #theme(panel.grid = element_blank())+
    theme_classic(base_size = size)+
    coord_cartesian(ylim=c(0,0.8))+
    theme.all+
    ylab("Prop. chromosomes clustered")

# coeff of disp
coeff.dispersion <- cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group.new, comparison) %>%
  summarise(disp.out = mean(disp.out))%>%
    ggplot(aes(x = group.new, y = disp.out, fill = group.new, color = group.new))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    geom_jitter(size = 1)+
    geom_boxplot(outlier.size = 0, weight=1, color = 1, linetype = 1)+
    theme_classic(base_size = size)+
    coord_cartesian(ylim=c(0,25))+
    theme.all+
    ylab("Outlier dispersion coefficient")

# nnd.diff

nnd.diff <- cluster.df %>%
  filter(num.outliers >= 3) %>%
  mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
  group_by(group.new, comparison) %>%
  summarise(nnd.diff = mean(nnd.diff))%>%
    ggplot(aes(x = group.new, y = nnd.diff, fill = group.new, color = group.new))+
    scale_color_manual(values = pal)+
    scale_fill_manual(values = pal) + 
    geom_jitter(size = 1)+
    geom_boxplot(outlier.size = 0, width=1, color = 1, linetype = 1)+
    theme_classic(base_size = size) +
    coord_cartesian(ylim=c(-5,5))+
    theme.all+
    ylab("Expected NND - Outlier NND (cM)")

# plot all the things
#grid.arrange(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep, ncol = 2)
grid_arrange_shared_legend(prop.clustered, nnd.diff, coeff.dispersion, nnd.rep)

quartz.save("figures/figureS2.pdf", type = "pdf")
