# figure 1

# what do the empirical regression lines look like?

groups <- c("para_S","allo_S", "allo_D","para_D")
pal <- wes_palette("Zissou", 50, type = "continuous")[c(1,17,30,50)]

dev.off()
dev.new()
plot.new()
par(mgp=c(3.5,1,0), mar=c(6,6,6,6), bty = "l")
plot(1, ylim = c(0,0.12), xlim=c(0,50), 
     xlab = "Recombination Rate (cM/Mb)", ylab = c(expression('F'["ST"]*" Outlier Probability")), 
     cex.lab = 1.5, cex.axis = 1.25, yaxs= "i", xaxs = "i")

for (i in 1:length(groups)){
  
  coeff.dat.group <- coeff.dat %>%
    filter(group == groups[i])
  
  intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
  slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
  
  #plot the first function
  
  eq <- function(x){
    ylog <- intercept + slope*x 
    return(1 - 1/(1 + exp(ylog))+0.0005)
  } 
  
  col.line <- pal[i]
  
  curve(eq, from = 0, to = 50, ylim = c(0,0.12), add=TRUE, lwd = 10, col = col.line)
  
}

legend(30,0.12, # places a legend at the appropriate place 
       c("Parapatric Same","Allopatric Same","Allopatric Different","Parapatric Different"), # puts text in the legend
       lty = c(1,1), # gives the legend appropriate symbols (lines)
       lwd = 10, col = pal, box.lty = 0, cex = 1.25)