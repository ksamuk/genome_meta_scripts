filter_data_call_outliers_function <- function (matched.all){
	
	matched.all <- matched.all %>%
		filter(!is.na(fst)) %>%
		filter(!is.infinite(fst))
	
	matched.all$fst[matched.all$fst < 0] <- 0 
	matched.all$fst[matched.all$dxy < 0] <- 0 
	
	matched.all <- matched.all %>%
		mutate(fst.outlier = is.outlier(fst))%>%
		mutate(dxy.outlier = is.outlier(dxy))%>%
		mutate(both.outlier = fst.outlier & dxy.outlier)%>%
		mutate(hs = (hexp1+hexp2)/2)
	
	return(matched.all)
}