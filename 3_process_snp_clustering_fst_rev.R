# process a folder of stats files, compute clustering, perform permutation test
# ks august 2015

#
rm(list = ls())

################################################################################
# Libraries and initalizing variables
################################################################################

library("dplyr")
library("IRanges")
library("data.table")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

################################################################################
# input files and output locations
################################################################################

# find raw stats files
stats.folder <- "stats/snp_all"
stats.files <- list.files(stats.folder, full.names = TRUE)

################################################################################
# Main clustering analysis function
################################################################################

calculate_coeff_dispersion_stats_file <- function(stats.filename, num_permutations = 10000){
	
	#### read in file and filter 
	
	print(paste0("Processing ",stats.filename,"..."))
	
	stats.file <- data.table(read.table(stats.filename, stringsAsFactors = FALSE, header=TRUE))
	
	setnames(stats.file, tolower(names(stats.file)))
	stats.file <- stats.file[ , chrom.pos:=NULL]
	setnames(stats.file, 1, "lg")
	stats.file$lg <- chrom.to.num(stats.file$lg)
	
	# filter FST values and call
	
		stats.file <- stats.file %>%
			filter(!is.na(fst)) %>%
			filter(!is.infinite(fst))
			
		stats.file$fst[stats.file$fst < 0] <- 0
	
		stats.file <- stats.file %>%
			mutate(fst.outlier = is.outlier(fst))

	# downsample wgs datasets to 100000 markers (otherwise they massively overload all clustering permutations)
		
	if(nrow(stats.file) > 100000){
		stats.file <- stats.file %>% sample_n(100000)
	}
	
	### add map distances
	stats.file <- add_map_distance(stats.file)
	
	### calculate coefficients of dispersion for each lg
	dispersion.stats <- list()
	for (j in unique(stats.file$lg)){
		
		stats.file.lg <- stats.file %>%
			filter(lg == j)
		
		disp.lg.all <- calculate_dispersion(stats.file.lg)
		disp.lg.outl <- calculate_dispersion_outliers(stats.file.lg)
		
		dispersion.stats[[j]] <- data.frame(lg = j, disp.all = disp.lg.all, disp.out = disp.lg.outl)
	}
	
	#### bind dispersion estimates into a df
	disp.df <- do.call("rbind", dispersion.stats)
	nnd.df <- calculate_nndist_all_lg(stats.file, num_permutations)
	
	#### format the cluster df for output
	cluster.df <- left_join(nnd.df, disp.df, by = "lg")
	
	file.name.stripped <- sapply(strsplit(stats.filename, split = "/"), function(x)gsub(".txt","",x[length(x)]))
	file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
	cluster.df$pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
	cluster.df$ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
	cluster.df$pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
	cluster.df$ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
	cluster.df$geography <- file.name.split[3]
	cluster.df$ecology <- file.name.split[4]
	
	cluster.df <- cluster.df %>% select(pop1, ecotype1, pop2, ecotype2, geography, ecology, everything())
	date.stamp <- paste("_", format(Sys.time(), "%Y-%m-%d"), sep="")
	out.file.name <- file.path(out.folder, paste(file.name.stripped, date.stamp, ".clustered.txt", sep=""))
	
	#### write to file
	write.table(cluster.df, file = out.file.name, row.names = FALSE, quote = FALSE)
}

################################################################################
# run clustering stats function
################################################################################

out.folder <- "analysis_ready/clustering_fst_new"
out.files <- list.files(out.folder)
dir.create(out.folder)

out.files.exist <- out.files %>% gsub("\\d","",.) %>% gsub("-","",.) %>% gsub(".gz_.clustered.txt",".txt.gz",.)
stats.reformat <- list.files(stats.folder)

stats.files <- rev(stats.files[!stats.reformat %in% out.files.exist])

################################################################################
# run clustering stats function
################################################################################

#cluster.master <- lapply(stats.files, calculate_coeff_dispersion_stats_file) 

mclapply(stats.files, calculate_coeff_dispersion_stats_file, mc.cores = 2, mc.silent = FALSE, mc.preschedule = FALSE) 
#cluster.master <- do.call("rbind", cluster.master)

#$date.stamp <- paste("_", format(Sys.time(), "%Y-%m-%d"), sep="")
#out.file.name <- file.path("analysis_ready", paste("clustering_master", date.stamp, ".txt", sep=""))
#write.table(cluster.master,file = out.file.name, row.names = FALSE, quote = FALSE)


