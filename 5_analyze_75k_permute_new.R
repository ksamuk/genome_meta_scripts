
rm(list=ls())

################################################################################
# Libraries and initalizing variables
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
# plotting parameters + global variables
################################################################################

theme_all <- theme_classic(base_size = 16)+
	theme(strip.text.x = element_blank(), 
				strip.background = element_blank(), 
				axis.title.x = element_blank(), 
				axis.title.y = element_blank(),
				#plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
				plot.margin=unit(c(0,0,0,0),"cm"),
				legend.position = "none")

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

n_permutations <- 10000
pvals <- list()
plots <- list()

################################################################################
# initialize coefficient data files
################################################################################

coeff.dat <- initialize_coeff_dat_files()

################################################################################
# perform permutations + plot results
################################################################################

stat <- "recomb_rate_fst"
group_type <- "group2.new" # relaxed

permutation_output <- run_cluster_permutations(coeff.dat, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

stat <- "recomb_rate_fst"
group_type <- "group.new" # strict

permutation_output <- run_cluster_permutations(coeff.dat, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

stat <- "recomb_rate_dxy"
group_type <- "group2.new" # relaxed

permutation_output <- run_cluster_permutations(coeff.dat, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

stat <- "recomb_rate_dxy"
group_type <- "group.new" # strict

permutation_output <- run_cluster_permutations(coeff.dat, stat, group_type, n_permutations)
pvals[[paste(stat,group_type)]] <- save_pvals(permutation_output)
plots[[length(plots)+1]] <- plot_permutation_output(permutation_output, stat, pal = pal, theme_all = theme_all)

################################################################################
# plot figure S5
################################################################################

pdf(file = "figures/Figure S5.pdf", height = 8.5, width = 12)
plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]] , 
					ncol = 1)
dev.off()

################################################################################
# write pvalue table to file
################################################################################

pval.out <- pvals %>% bind_rows
pval.out$group <- gsub("\\n"," ",pval.out$group)
pval.out[pval.out$group_type=="group2.new",]$group_type <- "relaxed"
pval.out[pval.out$group_type=="group.new",]$group_type <- "strict"
pval.out <- pval.out %>%
	mutate(effect.size = observed.means - permuted.means)

write.table(pval.out, file = "meta_data/recomb_permutation_pvalues.txt", row.names = FALSE, quote = FALSE)

