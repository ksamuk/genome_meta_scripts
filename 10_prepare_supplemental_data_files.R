# make supplemental datafiles 1 & 2

dat.fst <- read.table(file = "analysis_ready/75k_stats_model_fits_fst.txt", header = TRUE, stringsAsFactors = FALSE)
dat.dxy <- read.table(file = "analysis_ready/75k_stats_model_fits_dxy.txt", header = TRUE, stringsAsFactors = FALSE)
dat.fst.dxy <- read.table(file = "analysis_ready/75k_stats_model_fits_fst_dxy.txt", header = TRUE, stringsAsFactors = FALSE)

dat.fst <- with(dat.fst, data.frame(dat.fst[,1:4],
																		region1 = reg1,
																		region2 = reg2,
												            geography_strict = geography, 
																		geography_relaxed = geography2,
																		group_strict = group,
																		group_relaxed = group2,
																		dat.fst[,7:10]))

dat.dxy <- with(dat.fst, data.frame(dat.dxy[,1:4],
																		region1 = reg1,
																		region2 = reg2,
																		geography_strict = geography, 
																		geography_relaxed = geography2,
																		group_strict = group,
																		group_relaxed = group2,
																		dat.fst[,7:10]))

dat.fst.dxy <- with(dat.fst.dxy, data.frame(dat.fst.dxy[,1:4],
																		region1 = reg1,
																		region2 = reg2,
																		geography_strict = geography, 
																		geography_relaxed = geography2,
																		group_strict = group,
																		group_relaxed = group2,
																		dat.fst.dxy[,7:10]))

write.table(dat.fst, file = "figures/supplemental_data_1_model_fits_fst.txt", quote = FALSE, row.names = FALSE)
write.table(dat.dxy, file = "figures/supplemental_data_2_model_fits_dxy.txt", quote = FALSE, row.names = FALSE)
write.table(dat.fst.dxy, file = "figures/supplemental_data_3_model_fits_fst_dxy.txt", quote = FALSE, row.names = FALSE)

