############################################################
# calculate euclidian distance between populations
# KS Aug 2015
############################################################

###############
# LIBRARIES
###############

library(marmap)
library(dplyr)

###############
# INPUT FILES
###############

pop.dat <- read.csv("meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
stats.file.names <- list.files("stats/75k_all")[1:10]

###############
# BODY
###############

# Create nice looking color palettes
blues <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2", "lightsteelblue1")
greys <- c(grey(0.6), grey(0.93), grey(0.99))

stats.file.name <- stats.file.names[1]

compute_lc_distance_stats_file <- function (stats.file.name){
  
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
    filter(pop_alias == pop1) %>%
    select(latitude, longitude) %>%
    unique %>%
    unlist %>% 
    as.numeric
  pop1.lat <- pop1.coord[1]
  pop1.lon <- pop1.coord[2]
  
  pop2.coord <- pop.dat %>%
    filter(pop_alias == pop2) %>%
    select(latitude, longitude) %>%
    unique %>%
    unlist %>% 
    as.numeric
  pop2.lat <- pop2.coord[1]
  pop2.lon <- pop2.coord[2]
  
  # lit cam
  pop1.lat <- 49.012722
  pop1.lon <- -122.778247
  
  # japan
  pop2.lat <- 43.05375
  pop2.lon <- 144.778247
  
  43.05375	144.89443
  
  
  # import NOAA coast data :o

  # north america
  bat.1 <- getNOAA.bathy(0, -45, 30, 80, res = 10, keep = TRUE, antimeridian = TRUE)
  bat.2 <- getNOAA.bathy(-180, 100, 30, 80, res = 10, keep = TRUE, antimeridian = FALSE)
  
  loc.1 <- data.frame( x = c(pop1.lon+360, pop2.lon), y = c(pop1.lat, pop2.lat) )
  loc.2 <- data.frame( x = c(pop1.lon, pop2.lon), y = c(pop1.lat, pop2.lat) )
  tr.1 <- trans.mat(bat.1, min.depth = -10, max.depth = -300)
  tr.2 <- trans.mat(bat.2, min.depth = -10, max.depth = -300)
  cost.1 <- lc.dist(tr.1, loc.1, res="dist")
  cost.2 <- lc.dist(tr.2, loc.2, res="path")
  plot(bat.1, image = TRUE, asp = 2, land = TRUE, deep=-3000, shallow=-100, step=10000, drawlabels = FALSE, bpal = list(c(min(bat.1,na.rm=TRUE), 0, blues), c(0, max(bat.1, na.rm=TRUE), greys)), lwd = 0.0)
  plot(bat.2, image = TRUE, asp=2, land = TRUE, deep=-4000, shallow=-1000, step=1000, drawlabels = FALSE, bpal = list(c(min(bat.2,na.rm=TRUE), 0, blues), c(0, max(bat.2, na.rm=TRUE), greys)), lwd = 0.1)
  
  dummy <- lapply(cost.1, lines, col = col2alpha("orange", 0.5), lwd =5, lty = 1) 
  points(loc, bg = "orange", cex = 2, pch = 2)
  
  dummy <- lapply(cost.2, lines, col = col2alpha("orange", 0.5), lwd =5, lty = 1) 
  points(loc, bg = "orange", cex = 2, pch = 2)
  
  
  row.out <- data.frame(pop1, ecotype1, pop2, ecotype2, geography, ecology, lc.distance, euc.distance)
  return(row.out)
  
}



# Import bathymetry
bat <- getNOAA.bathy(-100, -80, 22, 31, res = 1, keep = TRUE)

# Load location of individuals (these are NOT from Viricel 2012)
loc <- data.frame( x = c(-96.92707, -96.60861, -96.86875, -96.14351, -92.82518, -90.86053, -90.14208, -84.64081, -83.81274, -81.13277, -80.33498, -88.52732, -94.46049), y = c(25.38657, 25.90644, 26.57339, 27.63348, 29.03572, 28.16380, 28.21235, 26.71302, 25.12554, 24.50031, 24.89052, 30.16034, 29.34550) )

# Compute least cost paths between -5m and -300m. 
# Beware! Computation takes time with high resolution bathymetries!
tr <- trans.mat(bat, min.depth = -1, max.depth = -300)
cost <- lc.dist(tr, loc, res="path")

# Add least cost paths and the position of individuals
