# verify if there is a bias in MAF representation in the data set

################################################################################
# Libraries and initalizing variables
################################################################################

#library("dplyr")
library("ggplot2")
library('readr')
library('chunked')
library('hexbin')

list.files("shared_functions", full.names = TRUE) %>% sapply(.,source, verbose = FALSE, echo = FALSE) %>% invisible

file_name <- "meta_data/maf.csv"

for (i in 1:21){
	file_name_lg <- paste0("meta_data/MAF_lg",i,".csv")
	
	if (!file.exists(file_name_lg)){
		read_chunkwise(file_name, chunk_size = 10000) %>% 
			#mutate(lg = chrom.to.num(chrom)) %>% 
			mutate(lg = chrom.to.num(chrom))%>%
			select(-chrom) %>%
			select(lg, everything()) %>%
			filter(lg == i) %>%
			write_chunkwise(file_name_lg)
	}
	
	maf_dat <- read_csv(file_name_lg)
	
	names(maf_dat)[5] <- "ecology"
	maf_dat$ecotype1 <- lapply(strsplit(maf_dat$pop1, split = "_"), function(x) x[2])
	maf_dat$ecotype2 <- lapply(strsplit(maf_dat$pop2, split = "_"), function(x) x[2])
	maf_dat$pop1 <- lapply(strsplit(maf_dat$pop1, split = "_"), function(x) x[1])
	maf_dat$pop2 <- lapply(strsplit(maf_dat$pop2, split = "_"), function(x) x[1])
	
	maf_dat <- add_region_data(maf_dat)
	maf_dat <- data.frame(maf_dat)
	
	maf_dat %>%
		select(group2, MAF) %>%
		filter(MAF >= 0.05) %>%
			ggplot(aes(x = MAF)) +
			geom_histogram(binwidth = 0.01) +
			facet_wrap(~group2, scales = "free_y")

		ggsave(filename = paste0("meta_data/MAF_lg_hist",i,".pdf"))
	
}




