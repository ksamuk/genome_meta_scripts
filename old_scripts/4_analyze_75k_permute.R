
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
library("cowplot")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

select <- dplyr::select

################################################################################
# plotting parameters + global variables
################################################################################

theme_all <- theme_classic(base_size = 12)+
	theme(strip.text.x = element_blank(), 
				strip.background = element_blank(), 
				axis.title.x = element_blank(), 
				axis.title.y = element_blank(),
				#plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
				plot.margin = unit(c(0,0,0,0),"cm"),
				legend.position = "none",
				axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

pal <- c("#E7C11A", "#9BBD95", "#F21A00", "#3B9AB2")

#number of permutations for  tessts
n_permutations <- 10000

################################################################################
# initialize coefficient data files
################################################################################

coeff.dat <- initialize_coeff_dat_files()

################################################################################
# perform permutations + plot results
################################################################################

coeff.dat <- coeff.dat %>%
	filter(n_windows_fst > 200)
	
# build a list of combinations of groups (strict, relaxed) 
# and fit types (fst, fst_dxy, dxy) to permute over (and plot / save pvalues)

stat <- c("recomb_rate_fst", "recomb_rate_dxy", "recomb_rate_hs")
stat <- rep(stat, each = 2)
group_type <- c("group2.new", "group.new")
group_type <- rep(group_type, 3)
combo.df <- data.frame(stat, group_type, stringsAsFactors = FALSE)

# interate over combo.df (permute, extract p-values, plot)

# initialize results lists
pvals <- list()
plots <- list()
permutation_output <- list()

for (i in 1:nrow(combo.df)){
	
	permutation_output[[i]] <- run_cluster_permutations(coeff.dat, combo.df$stat[i], combo.df$group_type[i], n_permutations)
	pvals[[i]] <- save_pvals(permutation_output[[i]], combo.df$stat[i], combo.df$group_type[i])

}

for (i in 1:nrow(combo.df)){
	
	plots[[i]] <- plot_permutation_output(permutation_output[[i]], combo.df$stat[i], pal = pal, theme_all = theme_all)
	
}
################################################################################
# plot figure S5
################################################################################

# note this plots the "raw" figure, annotations seen in paper were made manually

labels <- c("fst_relaxed", "fst_strict",
						"dxy_relaxed", "dxy_strict",
						"hs_relaxed", "hs_strict")
plot_grid(plots[[1]], plots[[2]], 
					plots[[3]], plots[[4]], 
					plots[[5]], plots[[6]],
					ncol = 1, labels = labels, align = "hv")



################################################################################
# write pvalue table to file
################################################################################

pval.out <- pvals %>% bind_rows 
pval.out$group <- gsub("\\n"," ",pval.out$group)
pval.out$group_type <- gsub("group2.new", "relaxed", pval.out$group_type)
pval.out$group_type <- gsub("group.new", "strict", pval.out$group_type)

pval.out <- pval.out %>%
	mutate(effect.size = observed.means - permuted.means)

write.table(pval.out, file = "meta_data/recomb_permutation_pvalues.txt", row.names = FALSE, quote = FALSE)

