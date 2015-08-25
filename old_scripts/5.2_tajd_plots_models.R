# analyze a tajima's D file

# the tajd file

tajd.dat <- list.files(file.path(getwd(), "analysis_ready"), pattern = "taj", full.names = TRUE)[1]
tajd.dat <- read.table(tajd.dat, header = TRUE, stringsAsFactors = FALSE)

tajd.dat$pop <- paste0(tajd.dat$location, "_", tajd.dat$ecotype )

# add in gene flow information

pop.file <- list.files(file.path(getwd(), "ev_prep_scripts"), pattern = "grouplist", full.names=TRUE)[1]
pop.file <- read.table(pop.file)
names(pop.file) <- c("pop", "location", "ecotype")



# plot tajd vs. pos

tajd.dat %>%
  ggplot(aes(x = pos1, y = tajd, color = pop))+
  geom_smooth()+
  facet_wrap(~lg)


tajd.dat %>%
  ggplot(aes(x = log(recomb_rate+1), y = tajd, color = pop))+
  geom_smooth(se = FALSE, size= 1)