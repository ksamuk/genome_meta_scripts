library(mapdata)
library(maps)
library(maptools)
library(rworldmap)
library(dplyr)


map("worldHires",  xlim = c(-130, -100), ylim = c(50,60), resolution = 0)

# input files
pop.dat <- read.csv(file = "meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
pop.dat <- pop.dat %>%
	select(pop, ecotype, latitude, longitude, Region, N) %>%
	filter(N > 1) %>%
	filter(Region != "Russia")

region.dat <-  pop.dat %>%
	group_by(Region) %>%
	summarise(lat.min = min(latitude), lat.max = max(latitude), 
						lon.min = min(longitude), lon.max = max(longitude))

region.dat$adjust.lat <- c(2.05, 1.25, 2.55, 17, 1.9, 2, 5)
region.dat$adjust.lon <- c(6, 1, 4.5, 25, 2.5, 2, 8)



plot_pop_map("Alaska", pop.dat)
plot_pop_map("Denmark", pop.dat)
plot_pop_map("BC", pop.dat)
plot_pop_map("Japan", pop.dat)
plot_pop_map("Nova Scotia", pop.dat)

par(mfrow = c(3,3), mar = c(0,0,0,0))
lapply(region.dat$Region, plot_pop_map, pop.dat = pop.dat)

plot_pop_map <- function (region, pop.dat, manual = FALSE, adjust.lat = 0, adjust.lon = 0){

	region.dat.sub <- region.dat %>%
		filter(Region == region)
	
	if (manual == FALSE){
		adjust.lat <- region.dat.sub$adjust.lat
		adjust.lon <- region.dat.sub$adjust.lon
	}

	lat <- c(region.dat.sub$lat.min-adjust.lat, region.dat.sub$lat.max+adjust.lat) 
	lon <- c(region.dat.sub$lon.min-adjust.lon, region.dat.sub$lon.max+adjust.lon)
	# the map itself
	map("worldHires", ylim = lat, 
										xlim = lon,
										fill = TRUE,
										col = "grey",
										resolution = 0, mar = c(1,1,1,1))
	box()
	
	# population points
	
	region.dat <- pop.dat %>% 
		filter(Region == region) %>%
		filter(N > 1)
	
	region.coords <- list(x = region.dat$longitude, y = region.dat$latitude)
	points(region.coords, pch = 20, cex = 3, col = as.factor(region.dat$ecotype))
	map.centre <- c(mean(region.dat.sub$lat.min, region.dat.sub$lat.max), mean(region.dat.sub$lon.min, region.dat.sub$lon.max))
	map.scale(ratio = FALSE)
	#map.scale(map.centre[2], map.centre[1], units = "Km")
	
}






# a nothern hemisphere map (not run)

# map parameters
ocean.col <- rgb(213,229,255, max = 255)
grid.col <-  rgb(147,157,172, max = 255)
margins <- c(0, 0, 0, 0)
lat.extent <- c(0,70)
long.extent <- c(0, 360)
resolution <- 1

par(mar = margins)
m <- map("worldHires",plot=FALSE)
map("worldHires", 
    projection="azequalarea", 
    ylim = lat.extent, 
    fill = FALSE,
    col = NA,
    bg = ocean.col,
    resolution = 1)

map.grid(lim = c(0,360,20,90),
         pretty = TRUE,
         lty = 1, 
         nx = 18,
         ny = 4,
         col = grid.col,
         labels = FALSE)


map.dat <- map("worldHires", 
    projection="azequalarea", 
    ylim = lat.extent + c(-2, 2), 
    add = TRUE,
    fill = TRUE,
    col = "white",
    bg = NULL,
    resolution = 1)



# adding points (not run)

pop.dat <- read.csv(file = "meta_data/populations_summary_table.csv", header = TRUE, stringsAsFactors = FALSE)
pop.dat <- pop.dat %>%
  select(pop, ecotype, latitude, longitude, Region, N) %>%


coords <- list(x = pop.dat$longitude, y = pop.dat$latitude)
points(mapproject(coords), pch=20, cex=1.2, col = as.factor(pop.dat$ecotype))



map.grid(m,
         pretty = TRUE,
         ylim = lat.extent + c(-2, 2),
         lty = 1, 
         nx = 1,
         ny = 20,
         col = "red",
         labels = FALSE)
