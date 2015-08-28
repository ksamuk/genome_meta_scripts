#### process clustering script output


## libraries
library("dplyr")

## input files
cluster.folder <- "analysis_ready/clustering_tmp"
cluster.files <- list.files(cluster.folder, full.names = TRUE)
header.row <- names (read.table(file = cluster.files[1], header = TRUE, stringsAsFactors = FALSE))

# combine all cluster
read_files_chopped <- function(cluster.file){
  return(read.table(file = cluster.file, stringsAsFactors = FALSE, skip = 1))
}

file.chopped <- lapply(cluster.files, read_files_chopped)
cluster.df <- bind_rows(file.chopped)
names(cluster.df) <-  header.row

# relaxed groups
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short pop codes
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

# associate pops in coeff dat with regions in region dat
region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop1","reg1")
cluster.df$reg1 <- region.sub$reg1[match(cluster.df$pop1, region.sub$pop1)]

region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop2","reg2")
cluster.df$reg2 <- region.sub$reg2[match(cluster.df$pop2, region.sub$pop2)]

# make new geographic categories
cluster.df$geography2 <- ifelse(cluster.df$reg1==cluster.df$reg2, "para", "allo")

#make new groups :o
cluster.df$group2 <- paste0(cluster.df$geography2,"_",cluster.df$ecology)

# elaborating groupings
cluster.df <- cluster.df %>%
  mutate(group = paste0(geography, "_", ecology)) %>%
  mutate(nnd.diff = nnd.mean.null -  nnd.mean.emp)

#grouped 
grouped.df <- cluster.df %>%
  mutate(group = paste0(geography, "_", ecology)) %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01) %>%
  group_by(group2) %>%
  summarise(num.cluster = mean(nnd.sig, na.rm = TRUE)) 
  
# dem plots?
cluster.df %>%
  ggplot(aes(x = group2, y = log(nnd.diff+1))) +
  geom_boxplot()

cluster.df %>%
  mutate(nnd.sig = nnd.emp.pvalue < 0.01)%>%
  filter(nnd.sig == TRUE) %>%
    ggplot(aes(x = group, y = log(nnd.diff+1))) +
    geom_jitter()

cluster.df %>%
  ggplot(aes(x = group2, y = disp.out+1)) +
  geom_boxplot()





