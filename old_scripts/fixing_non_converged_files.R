library("IRanges")
library("dplyr")

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible


coeff.dat <- initialize_coeff_dat_files()

files_sub <- coeff.dat %>% filter(recomb_rate_fst < -5) %>% select(comparison) %>% unlist %>% as.character

files_slug <- list.files("stats/75k_all") %>% 
	strsplit(split = ".(allo|para)") %>%
	lapply(function(x)x[1]) %>%
	unlist %>%
	gsub("_","\\.",.)

ev.dir <- file.path(getwd(),"evs")
ev.files <- list.files(ev.dir, pattern = "txt",full.names = TRUE)

stats_files_all <- list.files("stats/75k_all", full.names = TRUE)

stats_files_good <- list.files("stats/75k_all", full.names = TRUE)[!(files_slug %in% files_sub)] 
stats_matched <- match_evs(stats_files_good[1:20], linear_model_function = fit_linear_model_dxy)

tmp <- lapply(stats_files_all[8:12], match_evs, linear_model_function = fit_linear_model_dxy)
#linear_model_warnings <- lapply(tmp, function(x)x[12]) %>% as.character %>% unlist
coeff.dat <- bind_rows(tmp)


matched.all %>%
	ggplot(aes (x = recomb_rate, y = as.numeric(fst.outlier))) +
	geom_point()+
	geom_smooth()