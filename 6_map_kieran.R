library(ggplot2)
library(sp)
library(rgeos)
library(maptools)

data(wrld_simpl) # use this since it is already a SpatialPolygons object

# make a bounding box for the northern hemisphere
norths <- 35 # the minimum latitude?
ndiscr <- 50 # does nothing?
easts <- 0 # also does nothing
tweak <- 90

g <- gridlines(wrld_simpl, easts=easts, norths=norths , ndiscr = ndiscr)
x <- c(coordinates(g)[[1]][[1]][,1], rev(coordinates(g)[[1]][[1]][,1]), coordinates(g)[[1]][[1]][1])
y <- c(rep(tweak, ndiscr), rep(norths, ndiscr), tweak)
bnds <- cbind(x=x, y=y)
SP <- SpatialPolygons(list(Polygons(list(Polygon(bnds)), "1")), proj4string=CRS(proj4string(wrld_simpl)))

# the code snippet below is courtesy of Roger Bivand https://stat.ethz.ch/pipermail/r-sig-geo/2012-June/015340.html

gI <- gIntersects(wrld_simpl, SP, byid=TRUE)
out <- vector(mode="list", length=length(which(gI)))
ii <- 1
for (i in seq(along=gI)) if (gI[i]) {
  out[[ii]] <- gIntersection(wrld_simpl[i,], SP)
  row.names(out[[ii]]) <- row.names(wrld_simpl)[i]; ii <- ii+1
}
out1 <- do.call("rbind", out)

# fortify SP object to plot wiht ggplot
mymap <- fortify(out1, region='id') 
mymap <- mymap[order(mymap$order), ]

ggplot(data = mymap, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill='black') +
  #geom_path(color='gray60',aes(x=long,y=lat),size=0.1) +
  coord_map("orthographic")