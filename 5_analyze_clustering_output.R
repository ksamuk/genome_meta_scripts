## permutation test to assess significance of clustering data

rm(list = ls())

################################################################################
# Libraries and shared functions
################################################################################

library("ggplot2")
library("dplyr")
library("wesanderson")
library("grid")
library("gridExtra")
library("ggthemes")
library("lazyeval")
library("Hmisc")

list.files("shared_functions", full.names = TRUE) %>% sapply(source)

select <- dplyr::select

################################################################################
# Read in and format input files
################################################################################

cluster.df <- initialize_clustering_output()

# filter out linakge groups with super low data
# clustering metrics behave strangely when number of sites / number of outliers is very low
cluster.df <- cluster.df %>%
	filter(n.sites > 30) %>%
	filter(num.outliers > 5)

pvals <- list()
plots <- list()

theme_all <- theme_classic(base_size = 12)+
		theme(strip.text.x = element_blank(), 
				strip.background = element_blank(), 
				axis.title.x = element_blank(), 
				axis.title.y = element_blank(),
				#plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
				plot.margin = unit(c(0,0,0,0),"cm"),
				legend.position = "none"
				#axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
				)

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

n_permutations <- 10000

################################################################################
# Permutation test : nnd.diff, relaxed
################################################################################

stat <- "nnd.diff.sd"
group_type <- "group2.new" # relaxed

permutation_output <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group2.new, comparison) %>%
	summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	filter(!is.na(nnd.diff.sd)) %>%
	filter(nnd.diff.sd > -4)%>% # remove single crazy low value (some kind of error)
	run_cluster_permutations(., stat, group_type, n_permutations, collapsed = TRUE)

pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# Permutation test : nnd.diff, strict
################################################################################

stat <- "nnd.diff.sd"
group_type <- "group.new" # relaxed

permutation_output <- cluster.df %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group.new, comparison) %>%
	summarise(nnd.diff.sd = mean(nnd.diff.sd, na.rm = TRUE))%>%
	filter(!is.na(nnd.diff.sd)) %>%
	filter(nnd.diff.sd > -4)%>% # remove single crazy low value (some kind of error)
		run_cluster_permutations(., stat, group_type, n_permutations, collapsed = TRUE)

pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# Permutation test : dispersion, relaxed
################################################################################

stat <- "disp.out"
group_type <- "group2.new" # relaxed

permutation_output <- run_cluster_permutations(cluster.df, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# Permutation test : dispersion, strict
################################################################################

stat <- "disp.out"
group_type <- "group.new" # strict

permutation_output <- run_cluster_permutations(cluster.df, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# Permutation test : number of cluster chromosomes, relaxed
################################################################################

stat <- "more.clustered"
group_type <- "group2.new" # strict

# calculate the number of over clustered and under clustered lgs 
# using null distribution of site differences (all snps), finds the percentile where the 
# empirical nnd falls ("nnd.emp.percentile"), then sums up the number of those in the upper and lower %5
# lower %5 = outliers have shorter nnd than expected, upper %5 = the opposite

cluster.df.collapsed <- cluster.df %>%
	mutate(more.clustered = nnd.emp.percentile < 0.05) %>%
	mutate(less.clustered = nnd.emp.percentile > 0.95) %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group2.new, comparison) %>%
	summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
	mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
	filter(num.lg > 1)

permutation_output <- run_cluster_permutations(cluster.df.collapsed, stat, group_type, n_permutations, collapsed = TRUE)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)


################################################################################
# Permutation test : number of cluster chromosomes, relaxed
################################################################################

stat <- "more.clustered"
group_type <- "group.new" # strict

# calculate the number of over clustered and under clustered lgs 
# using null distribution of site differences (all snps), finds the percentile where the 
# empirical nnd falls ("nnd.emp.percentile"), then sums up the number of those in the upper and lower %5
# lower %5 = outliers have shorter nnd than expected, upper %5 = the opposite

cluster.df.collapsed <- cluster.df %>%
	mutate(more.clustered = nnd.emp.percentile < 0.05) %>%
	mutate(less.clustered = nnd.emp.percentile > 0.95) %>%
	mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
	group_by(group.new, comparison) %>%
	summarise(more.clustered = sum(more.clustered), less.clustered = sum(less.clustered), num.lg = sum(lg <= 21)) %>%
	mutate(more.clustered = more.clustered / num.lg, less.clustered = less.clustered / num.lg) %>%
	filter(num.lg > 1)

permutation_output <- run_cluster_permutations(cluster.df.collapsed, stat, group_type, n_permutations, collapsed = TRUE)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output, stat, group_type)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)


################################################################################
# make unified plot
################################################################################

pdf(file = "figures/FigureS6.pdf", height = 8.5, width = 12, onefile = FALSE)
labels <- c("num.clustered_relaxed", "num.clustered_strict", 
						"nnd.diff_relaxed", "nnd.diff_strict",
						"disp_relaxed", "disp_strict"
						)
plot_grid(plots[[5]], plots[[6]], plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
					ncol = 1, labels = labels, align = "hv")
dev.off()

################################################################################
# format and write pvalue table to file
################################################################################

pval.out <- pvals %>% bind_rows
pval.out$group <- gsub("\\n"," ",pval.out$group)
pval.out[pval.out$group_type=="group2.new",]$group_type <- "relaxed"
pval.out[pval.out$group_type=="group.new",]$group_type <- "strict"
pval.out <- pval.out %>%
	mutate(effect.size = observed.means - permuted.means)

write.table(pval.out, file = "meta_data/clustering_permutation_pvalues_fst.txt", row.names = FALSE, quote = FALSE)


