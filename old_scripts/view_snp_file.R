library("readr")
library("dplyr")
library("ggplot2")

list.files("shared_functions", full.names = TRUE) %>% sapply(source) %>% invisible

stat_df <- read_delim(file = "stats/snp_all/pri_limnetic.pri_benthic.para.D.stats.txt.gz", delim = "\t", 
											col_types = cols("CHROM-POS" = "c", CHROM = "c", POS = "d", N1 ="d", N2 ="d", NTotal ="d", Dxy ="d", FstNum = "d", FstDenom = "d", Fst = "d", Hexp1 = "d", Hexp2 ="d"))


stat_df <- stat_df %>%
	filter(!is.infinite(Fst)) %>%
	filter(!is.na(Fst)) %>%
	mutate(hs = (Hexp1 + Hexp2) / 2)

stat_df$Fst[stat_df$Fst < 0] <- 0

stat_df %>%
	mutate(fst.outlier = is.outlier(Fst)) %>%
	sample_frac(0.1) %>%
	ggplot(aes(x = Dxy, y = FstDenom, color = fst.outlier))+
	geom_point()+
	geom_smooth()