plot_averaged_regression_lines <- function(coeff.dat, group_order, pal, group_variable, ylim = c(0,1), xlim = c(0,1), ylab = ylab){
	
	# set up the plot
	label.scaling <-  0.8
	line.weight <- 6
	par(mgp = c(2.5, 1, 0), mar = c(5, 4, 2, 1), mex = 1, bty = "l")
	plot(1, ylim = ylim, xlim = xlim, 
			 xlab = "Recombination Rate (cM/Mb)", ylab = ylab, 
			 cex.lab = label.scaling, cex.axis = label.scaling, yaxs= "i", xaxs = "i")
	text(49, 0.127, "A")
	
	# plot each line individually
	for (i in 1:length(group_order)){
		
		# subset data by group
		filter_criteria <- interp(~ which_column == group_order[i], which_column = as.name(group_variable))
		coeff.dat.group <- coeff.dat %>%
			filter_(filter_criteria)
		
		# average intcepts and slopes
		intercept <- mean(coeff.dat.group$intercept, na.rm = TRUE)
		slope <- mean(coeff.dat.group$recomb_rate, na.rm = TRUE)
		
		#define an equation function
		#0.0005 is for style only (keeps line from overplotting axis)
		eq <- function(x){
			ylog <- intercept + slope*x 
			return(1 - 1/(1 + exp(ylog)) + 0.0005)
		} 
		
		col.line <- pal[i]
		curve(eq, from = 0, to = 50, ylim = ylim, add=TRUE, lwd = line.weight, col = col.line)
		
	}
	
}