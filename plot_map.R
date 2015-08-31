library(mapproj)
library(rgdal)
library(mapdata)
library(dplyr)

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

lat <- c(49) 
lon <- c(-123)

coords <- list(x = -123, y = 49)
points(mapproject(coords), pch=20, cex=1.2, col = "red")



map.grid(m,
         pretty = TRUE,
         ylim = lat.extent + c(-2, 2),
         lty = 1, 
         nx = 1,
         ny = 20,
         col = "red",
         labels = FALSE)
