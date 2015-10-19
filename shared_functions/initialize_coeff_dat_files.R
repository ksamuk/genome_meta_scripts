initialize_coeff_dat_files<- function(){

	# choose either fst or fst/dxy 
	coeff.dat.fst <- read.table(file = "analysis_ready/75k_stats_model_fits_fst.txt", header = TRUE, stringsAsFactors = FALSE)
	coeff.dat.fst.dxy <- read.table(file = "analysis_ready/75k_stats_model_fits_fst_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
	coeff.dat.dxy <- read.table(file = "analysis_ready/75k_stats_model_fits_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
	
	coeff.dat.fst <- coeff.dat.fst %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		rename(recomb_rate_fst = recomb_rate) %>%
		select(comparison, pop1, ecotype1, pop2, ecotype2, 
					 reg1, reg2, geography, geography2, ecology, 
					 group, group2, recomb_rate_fst)
	
	coeff.dat.fst.dxy <- coeff.dat.dxy %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		rename(recomb_rate_fst_dxy = recomb_rate) %>%
		select(comparison, recomb_rate_fst_dxy)
	
	coeff.dat.dxy <- coeff.dat.dxy %>%
		mutate(comparison = paste(pop1, ecotype1, pop2, ecotype2, sep = ".")) %>%
		rename(recomb_rate_dxy = recomb_rate) %>%
		select(comparison, recomb_rate_dxy)
	
	coeff.dat <- left_join(coeff.dat.fst, coeff.dat.dxy)
	coeff.dat <- left_join(coeff.dat, coeff.dat.fst.dxy)
	
	## add in region data (for looser geography)
	region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)
	
	# make short pop codes
	region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist
	
	# associate pops in coeff dat with regions in region dat
	region.sub <- region.dat[,c(1,4)]
	names(region.sub) <- c("pop1","reg1")
	coeff.dat$reg1 <- region.sub$reg1[match(coeff.dat$pop1, region.sub$pop1)]
	
	region.sub <- region.dat[,c(1,4)]
	names(region.sub) <- c("pop2","reg2")
	coeff.dat$reg2 <- region.sub$reg2[match(coeff.dat$pop2, region.sub$pop2)]
	
	# make new geographic categories
	coeff.dat$geography2 <- ifelse(coeff.dat$reg1==coeff.dat$reg2, "para", "allo")
	
	#make new groups :o
	coeff.dat$group2 <- paste0(coeff.dat$geography2,"_",coeff.dat$ecology)
	
	group.old.names <- c("allo_D","allo_S", "para_D", "para_S")
	group.rename <- c("Allopatry\nDivergent", "Allopatry\nParallel", "Gene Flow\nDivergent", "Gene Flow\nParallel")
	coeff.dat$group2.new <- group.rename[match(coeff.dat$group2, group.old.names)]
	coeff.dat$group.new <- group.rename[match(coeff.dat$group, group.old.names)]
	return(coeff.dat)
}