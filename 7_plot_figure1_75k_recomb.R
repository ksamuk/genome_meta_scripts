##### FIGURES 1 & S1
##### This script generates the two above figures. Note that the Quartz graphic device requires a Mac OS X.

# figure 1
rm(list=ls())

library("ggplot2")
library("dplyr")
#library("devtools")
#install_github("karthik/wesanderson")
library("wesanderson")

############## FIG 1: STRICT GROUPINGS 

dev.off()
dev.new()
quartz(width = 8.5, height = 6)

##### what do the empirical regression lines look like? - FST ONLY

coeff.dat.fst <- read.table(file = "analysis_ready/75k_stats_model_fst_fits.txt", header = TRUE, stringsAsFactors = FALSE)

groups <- c("para_S","allo_S", "allo_D","para_D")
pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]

label.scaling = 0.8
line.weight = 6
layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(6, 6))
par(mgp=c(2.5,1,0), mar=c(5,4,2,1), mex = 1, bty = "l")
plot(1, ylim = c(0,0.13), xlim=c(0,50), 
     xlab = "Recombination Rate (cM/Mb)", ylab = c(expression('F'["ST"]*" Outlier Probability")), 
     cex.lab = label.scaling, cex.axis = label.scaling, yaxs= "i", xaxs = "i")
text(49,0.127, "A")

for (i in 1:length(groups)){
  
  coeff.dat.group <- coeff.dat.fst %>%
    filter(group == groups[i])
  
  intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
  slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
  
  #plot the first function
  
  eq <- function(x){
    ylog <- intercept + slope*x 
    return(1 - 1/(1 + exp(ylog))+0.0005)
  } 
  
  
  col.line <- pal[i]
  curve(eq, from = 0, to = 50, ylim = c(0,0.13), add=TRUE, lwd = line.weight, col = col.line)
  
}

##### what do the empirical regression lines look like? - FST/DXY

coeff.dat.dxy <- read.table(file = "analysis_ready/75k_stats_model_dxy_fits.txt", header = TRUE, stringsAsFactors = FALSE)

plot(1, ylim = c(0,0.02), xlim=c(0,50), 
     xlab = "Recombination Rate (cM/Mb)", ylab = c(expression('Joint F'["ST"]*'/D'["XY"]*" Outlier Probability")), 
     cex.lab = label.scaling, cex.axis = label.scaling, yaxs= "i", xaxs = "i")
text(49,0.0195, "B")

for (i in 1:length(groups)){
  
  coeff.dat.group <- coeff.dat.dxy %>%
    filter(group == groups[i])
  
  intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
  slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
  
  #plot the first function
  
  eq <- function(x){
    ylog <- intercept + slope*x 
    return(1 - 1/(1 + exp(ylog))+0.00005)
  } 
  
  col.line <- pal[i] 
  
  curve(eq, from = 0, to = 50, ylim = c(0,0.02), add=TRUE, lwd = line.weight, col = col.line)
  
}


legend(23,0.02, # places a legend at the appropriate place 
       c("Gene Flow Parallel","Allopatry Parallel","Allopatry Divergent","Gene Flow Divergent"), 
       lty = c(1,1), # gives the legend appropriate symbols (lines)
       lwd = line.weight, col = pal, box.lty = 0, cex = label.scaling, y.intersp = 1.5, xjust = 0, seg.len = 0.75)

quartz.save("figures/figure1.pdf", type = "pdf")

############## FIG S1: RELAXED GROUPINGS 

dev.off()
dev.new()
quartz(width = 8.5, height = 6)

##### what do the empirical regression lines look like? - FST ONLY

coeff.dat.fst <- read.table(file = "analysis_ready/75k_stats_model_fst_fits.txt", header = TRUE, stringsAsFactors = FALSE)

groups <- c("para_S","allo_S", "allo_D","para_D")
pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]

label.scaling = 0.8
line.weight = 6
layout(matrix(c(1,2), ncol=2, byrow=TRUE), heights=c(6, 6))
par(mgp=c(2.5,1,0), mar=c(5,4,2,1), mex = 1, bty = "l")
plot(1, ylim = c(0,0.13), xlim=c(0,50), 
     xlab = "Recombination Rate (cM/Mb)", ylab = c(expression('F'["ST"]*" Outlier Probability")), 
     cex.lab = label.scaling, cex.axis = label.scaling, yaxs= "i", xaxs = "i")
text(49,0.127, "A")
for (i in 1:length(groups)){
  
  coeff.dat.group <- coeff.dat.fst %>%
    filter(group2 == groups[i])
  
  intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
  slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
  
  #plot the first function
  
  eq <- function(x){
    ylog <- intercept + slope*x 
    return(1 - 1/(1 + exp(ylog))+0.0005)
  } 
  
  
  col.line <- pal[i]
  curve(eq, from = 0, to = 50, ylim = c(0,0.13), add=TRUE, lwd = line.weight, col = col.line)
  
}

##### what do the empirical regression lines look like? - FST/DXY

coeff.dat.dxy <- read.table(file = "analysis_ready/75k_stats_model_dxy_fits.txt", header = TRUE, stringsAsFactors = FALSE)

plot(1, ylim = c(0,0.02), xlim=c(0,50), 
     xlab = "Recombination Rate (cM/Mb)", ylab = c(expression('Joint F'["ST"]*'/D'["XY"]*" Outlier Probability")), 
     cex.lab = label.scaling, cex.axis = label.scaling, yaxs= "i", xaxs = "i")
text(49,0.0195, "B")
for (i in 1:length(groups)){
  
  coeff.dat.group <- coeff.dat.dxy %>%
    filter(group2 == groups[i])
  
  intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
  slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
  
  #plot the first function
  
  eq <- function(x){
    ylog <- intercept + slope*x 
    return(1 - 1/(1 + exp(ylog))+0.00005)
  } 
  
  col.line <- pal[i] 
  
  curve(eq, from = 0, to = 50, ylim = c(0,0.02), add=TRUE, lwd = line.weight, col = col.line)
  
}


legend(23,0.02, # places a legend at the appropriate place 
       c("Gene Flow Parallel","Allopatry Parallel","Allopatry Divergent","Gene Flow Divergent"), 
       lty = c(1,1), # gives the legend appropriate symbols (lines)
       lwd = line.weight, col = pal, box.lty = 0, cex = label.scaling, y.intersp = 1.5,xjust = 0, seg.len = 0.75)

quartz.save("figures/figureS1.pdf", type = "pdf")

