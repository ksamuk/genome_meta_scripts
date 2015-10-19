############################################################
# calculate euclidian & 'least cost' distance between populations
# KS Aug 2015
############################################################

rm(list =ls())

#########################
# LIBRARIES & FUNCTIONS
#########################

library(parallel)
library(marmap)
library(dplyr)
library(fossil)

source("shared_functions/lc_path_to_dist.R")
source("shared_functions/getNOAA_bathy_prefix.R")
source("shared_functions/lc_dist_no_bar.R")
source("shared_functions/lon_to_region.R")
source("shared_functions/point_to_nearest_coastline.R")

select <- dplyr::select

#########################
# INPUT FILES
#########################

pop.dat <- read.csv("meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
stats.file.names <- list.files("stats/75k_all")

#########################
# BODY
#########################

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

# preload and transform bathy data
# this is a big bottleneck if placed in the function
bat1 <- getNOAA_bathy_prefix(0, -45, 30, 80, res = 8, keep = TRUE, antimeridian = TRUE, prefix = "meta_data/")
bat2 <- getNOAA_bathy_prefix(180, -180, 30, 80, res = 8, keep = TRUE, prefix = "meta_data/")

trans.1.file <- "meta_data/bathy.1.trans"
trans.2.file <- "meta_data/bathy.2.trans"

if(!file.exists(trans.1.file)){
	tr1 <- trans.mat(bat1, min.depth = -10, max.depth = -800)
	save(tr1, file = trans.1.file)
} else{
	load(trans.1.file )
}

if(!file.exists(trans.2.file)){
	tr2 <- trans.mat(bat2, min.depth = -10, max.depth = -800)
	save(tr2, file = trans.2.file)
} else{
	load(trans.2.file )
}

compute_lc_distance_stats_file <- function (stats.file.name, bathy1, bathy2, trans1, trans2, plot.map = TRUE){
	start.time <- Sys.time()
	print(paste0("processing ", stats.file.name, "..."))
	
  # parse pop names
  file.name.stripped <- sapply(strsplit(stats.file.name, split = "/"), function(x)gsub(".txt","",x[length(x)]))
  file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
  pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
  pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
  geography <- file.name.split[3]
  ecology <- file.name.split[4]
  
  # find coordinates
  
  pop1.coord <- pop.dat %>%
    filter(pop_alias == pop1, ecotype == ecotype1) %>%
    select(latitude, longitude) %>%
    unique %>%
    unlist %>% 
    as.numeric
  pop1.lat <- pop1.coord[1]
  pop1.lon <- pop1.coord[2]
  
  pop2.coord <- pop.dat %>%
    filter(pop_alias == pop2, ecotype == ecotype2) %>%
    select(latitude, longitude) %>%
    unique %>%
    unlist %>% 
    as.numeric
  pop2.lat <- pop2.coord[1]
  pop2.lon <- pop2.coord[2]
  
  # set region from lon dat
  pop1.region <- lon_to_region(pop1.lon) 
  pop2.region <- lon_to_region(pop2.lon) 
  
  
  # import NOAA coast data
	# choose a map based on the locations, to prevent Inf distances due to wrapping
  if (sum(c(pop1.region, pop2.region) %in% c("jp","na")) >1) {
  	
      # na vs. na, na vs. jp, jp vs. jp
	  	bat <- bathy1
	  	tr <- trans1
	  	
	  	if (pop1.lon < 0){
	  		pop1.lon.adj <- pop1.lon + 360
	  	} else{
	  		pop1.lon.adj <- pop1.lon
	  	}
	  	if (pop2.lon < 0){
	  		pop2.lon.adj <- pop2.lon + 360
	  	} else{
	  		pop2.lon.adj <- pop2.lon
	  	}

	  	loc <- data.frame( x = c(pop1.lon.adj, pop2.lon.adj), y = c(pop1.lat, pop2.lat))
	  	
	  	# find point to nearest coastline for both pops
	  	nearest.coastline <- point_to_nearest_coastline(bat, loc, mode = 1)
	  	
	  	loc$x[1] <- nearest.coastline$loc.x.1.[1] 
	  	loc$x[2] <- nearest.coastline$loc.x.2.[1]
	  	loc$y[1] <- nearest.coastline$loc.y.1.[1]
	  	loc$y[2] <- nearest.coastline$loc.y.2.[1] 
	  	dist.to.coast1 <- nearest.coastline$dist.to.coast1
	  	dist.to.coast2 <- nearest.coastline$dist.to.coast2
  	
  	} else if (sum(c(pop1.region, pop2.region) %in% c("eu", "na")) >1){
	   
  		 #eu vs. eu, eu vs. na	
	  	bat <- bathy2
	  	tr <- trans2
	  	
	  	loc <- data.frame(x = c(pop1.lon, pop2.lon), y = c(pop1.lat, pop2.lat))
	  	
	  	# find point to nearest coastline for both pops
	  	nearest.coastline <- point_to_nearest_coastline(bat, loc, mode = 2)
	  	
	  	loc$x[1] <- nearest.coastline$loc.x.1.[1] 
	  	loc$x[2] <- nearest.coastline$loc.x.2.[1]
	  	loc$y[1] <- nearest.coastline$loc.y.1.[1]
	  	loc$y[2] <- nearest.coastline$loc.y.2.[1] 
	  	dist.to.coast1 <- nearest.coastline$dist.to.coast1
	  	dist.to.coast2 <- nearest.coastline$dist.to.coast2
	  	
  	} else if (sum(c(pop1.region, pop2.region) %in% c("jp","eu")) > 1){
  		
  		# jp vs. eu
  		bat <- bathy1
  		tr <- trans1
  		
  		if (pop1.lon < 0){
  			pop1.lon.adj <- pop1.lon + 360
  		} else{
  			pop1.lon.adj <- pop1.lon
  		}
  		if (pop2.lon < 0){
  			pop2.lon.adj <- pop2.lon + 360
  		} else{
  			pop2.lon.adj <- pop2.lon
  		}
  		
  		loc <- data.frame( x = c(pop1.lon.adj, pop2.lon.adj), y = c(pop1.lat, pop2.lat)) 
  		
  		# find point to nearest coastline for both pops
  		nearest.coastline <- point_to_nearest_coastline(bat, loc, mode =1)
  		
  		loc$x[1] <- nearest.coastline$loc.x.1.[1] 
  		loc$x[2] <- nearest.coastline$loc.x.2.[1]
  		loc$y[1] <- nearest.coastline$loc.y.1.[1]
  		loc$y[2] <- nearest.coastline$loc.y.2.[1] 
  		dist.to.coast1 <- nearest.coastline$dist.to.coast1
  		dist.to.coast2 <- nearest.coastline$dist.to.coast2
  		
  	}

  # find the least cost path between coastline points
  least.cost.path <- NA
  least.cost.distance <- NA
  
  try(least.cost.path <- lc_dist_no_bar(tr, loc), silent = TRUE)
  try(least.cost.distance <- lc_path_to_dist(least.cost.path), silent = TRUE)

  # find the great circle distances between the *original* points
  loc.orig <- data.frame( x = c(pop1.lon, pop2.lon), y = c(pop1.lat, pop2.lat)) 
  euc.distance <- earth.dist(loc.orig)[1] %>% round

  if (plot.map){

	  # plot bathy map with lc path + points
  	png(paste0("ev_prep_scripts/map_out/",stats.file.name,".map.png"))
  	plot(bat, image = TRUE, asp = 2, 
	  		 land = TRUE, deep= -100000, 
	  		 shallow=-100, step=100000, 
	  		 drawlabels = FALSE, 
	  		 bpal = list(c(min(bat,na.rm = TRUE), 0, blues), c(0, max(bat, na.rm = TRUE), greys)), 
	  		 lwd = 0.0)
  	
	  if (is.na(least.cost.path)){
	  	dummy <- lapply(loc, lines, col = col2alpha("orange", 0.5), lwd =5, lty = 1)
	  }else{
	  	dummy <- lapply(least.cost.path, lines, col = col2alpha("orange", 0.5), lwd =5, lty = 1)
	  }
	  points(loc, bg = "orange", cex = 1, pch = 16)
	  dev.off()
  }
  
  
  row.out <- data.frame(pop1, ecotype1, pop2, ecotype2, geography, ecology, 
  											least.cost.distance, euc.distance, dist.to.coast1, dist.to.coast2)
  
  time.elapsed <- (Sys.time() - start.time) %>% as.numeric
  print(paste0("Processing took ", time.elapsed, " seconds"))
  return(row.out)
  
  
}

# run the above function for all files, on three cores 
# might not work in *NIX?

cl <- makeCluster(getOption("cl.cores", 3))

clusterEvalQ(cl, {
	library(parallel)
	library(marmap)
	library(dplyr)
	library(fossil)
	source("shared_functions/lc_path_to_dist.R")
	source("shared_functions/getNOAA_bathy_prefix.R")
	source("shared_functions/lc_dist_no_bar.R")
	source("shared_functions/lon_to_region.R")
	source("shared_functions/point_to_nearest_coastline.R")
	pop.dat <- read.csv("meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
	stats.file.names <- list.files("stats/75k_all")
	blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
	greys <- c(grey(0.6), grey(0.93), grey(0.99))
})

distances.df <- parLapply(cl = cl, stats.file.names, compute_lc_distance_stats_file, 
													bathy1 = bat1, bathy2 = bat2, trans1 = tr1, trans2 = tr2,
													plot.map = TRUE)

stopCluster(cl)

# single core version
# works for all OSs
#distances.df <- lapply(stats.file.names, compute_lc_distance_stats_file, 
#											 bathy1 = bat1, bathy2 = bat2, trans1 = tr1, trans2 = tr2,
#											 plot.map = TRUE)

# clean up results and write to file
distances.df <- bind_rows(distances.df)

distances.df <- distances.df %>%
	unique

write.table(distances.df, file = "meta_data/pop_geo_distances.txt", 
						quote = FALSE, row.names = FALSE)
