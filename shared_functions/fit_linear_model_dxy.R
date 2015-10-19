fit_linear_model_dxy <- function (matched.all){
	
	model <- glm(as.numeric(dxy.outlier) ~ recomb_rate + ds + gene_count,
							 na.action = "na.omit",
							 family = binomial,
							 data = matched.all)
	
	return(model)
}