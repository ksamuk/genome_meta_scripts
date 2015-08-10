# process a folder of stats files, compute clustering, perform permutation test
# ks august 2015

# libraries

library("dplyr")
library("IRanges")
library("data.table")

# folders

stats.folder <- "stats/snp_all"
stats.files <- list.files(stats.folder, full.names = TRUE)

# functions

# chrom to numeric function
chrom.to.num <- function(x){
  x <- gsub("group", "", x)
  chrom.rom <- as.character(as.roman(c(1:21)))
  return(match(x, chrom.rom))
}

######## ADD MAP DISTANCES TO A SNP FILE

# adds map distance to a stats file

add_map_distance <- function(snp.file){
  
  # read in snp file
  
  snp.file$gen.pos <- NA
  snp.file$pos <- as.numeric(snp.file$pos)
  
  snp.file <- snp.file %>%
    arrange(lg,pos)
  
  # read in genetic map data
  map.file <- list.files(file.path("ev_prep_scripts"), pattern="roesti_recomb_estimates.txt", full.names = TRUE)
  map.file <- read.table(map.file, header = TRUE)
  names(map.file)[1] <- "lg"
  
  
  # calculate genetic distance
  
  #place holder vector
  
  snp.map.position <- c()
  
  for (i in unique(snp.file$lg)){
    
    snp.file.lg <- snp.file %>%
      filter(lg == i)
    
    map.file.lg <- map.file %>%
      filter(lg == i)
    
    snp.map.position.lg <- snp.file.lg$gen.pos
    
    for (k in 1:length(snp.file.lg$pos)){
      
      snp.pos <- snp.file.lg$pos[k]
      
      #snp.pos <- 1190613 - 12500
      
      # find matching window for snp
      snp.pos.window <- map.file.lg[map.file.lg$pos1 < snp.pos & map.file.lg$pos2 >= snp.pos,]
      
      if (length(snp.pos.window$lg)>0){
        
        # how far (proportionally) snp is along window in bp
        snp.window.position <- 1 - (snp.pos.window$pos2 - snp.pos) / (snp.pos.window$pos2 - snp.pos.window$pos1)
        
        # convert to genetic distance
        snp.map.position.lg[k] <- (snp.window.position * (snp.pos.window$gen2 - snp.pos.window$gen1)) + snp.pos.window$gen1
      } else{
        snp.map.position.lg[k] <- NA
      }
      
      
    }
    
    snp.map.position<-c(snp.map.position, snp.map.position.lg)
    
  }
  
  snp.file$gen.pos <- snp.map.position
  return(snp.file)
}

######## END ADD MAP DISTANCES FUNCTION

######## FIND OUTLIERS FUNCTION

is.outlier <- function(x){
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  return(x >=x95)
}

######## END FIND OUTLIERS FUNCTION

######## DISPERSION FUNCTIONS


# calc for all loci
calculate_dispersion <- function (lg) {
  
  dispersion <- lg %>%
    select(gen.pos) %>% 
    arrange(gen.pos) %>%
    unlist %>% 
    diff %>% 
    (function(x) return(var(x,na.rm = TRUE)/mean(x,na.rm = TRUE)))
  
  return(dispersion)
  
}

# calc for outliers
calculate_dispersion_outliers <- function (lg) {
  
  start <- max(lg$gen.pos)
  end <- min(lg$gen.pos)
  
  outlier.positions <- lg %>% 
    filter(both.outlier == TRUE)%>%
    select(gen.pos) %>%
    unlist
  outlier.positions <- c(start, outlier.positions, end)
  
  outlier.dispersion <- outlier.positions %>%
    sort %>%
    diff %>% 
    (function(x) return(var(x,na.rm = TRUE)/mean(x,na.rm = TRUE)))
  
  return(outlier.dispersion)
  
}

######## END DISPERSION FUNCTIONS

######## NND FUNCTIONS

calculate_nndist_all_lg <- function (stats.file) {
  
  ## first, build null distributions of nndists for each linkage group:
  
  nnd.stats <- list()
  
  for (j in unique(stats.file$lg)){
    
    # subset for lg j
    stats.file.lg <- stats.file %>%
      filter(stats.file$lg == j)
   
    # the number of outliers on that lg
    num.outliers <- sum(as.numeric(stats.file.lg$both.outlier), na.rm = TRUE)
    
    if (num.outliers > 1){
    
      # draw 10000 samples of num.outliers random loci, take the mean, and return the ecdf and mean
      null.mean.nnds <- rep(NA, 10000)
      for (i in 1:10000){
        
        site.sample <- stats.file.lg %>%
        filter(!is.na(gen.pos)) %>%
        select(gen.pos) %>%
        sample_n(num.outliers) %>%
        arrange(gen.pos) %>%
        mutate(dist.1 = c(NA,diff(gen.pos))) %>%
        mutate(dist.2 = c(diff(sort(gen.pos)),NA))
        
        nn.dist <- rep(NA, length(site.sample$genpos))
        for (k in 1:length(site.sample$gen.pos)){
          
          if(!is.na(site.sample$dist.1[k]) & !is.na(site.sample$dist.2[k])){
            nn.dist[k] <- min(c(site.sample$dist.1[k],site.sample$dist.2[k]))
          }else if(is.na(site.sample$dist.1[k])){
            nn.dist[k] <- site.sample$dist.2[k]
          } else if(is.na(site.sample$dist.2[k])){
            nn.dist[k] <- site.sample$dist.1[k]
          }
        }
        null.mean.nnds[i] <- mean(nn.dist, na.rm = TRUE)
      }
      
    # calculate the estimate mean null nndist
    null.mean <- mean(null.mean.nnds, na.rm = TRUE)
    null.ecdf <- ecdf(null.mean.nnds)
    
    # calculate the empirical nndist for real outliers
    
    site.sample <- stats.file.lg %>%
      filter(!is.na(gen.pos)) %>%
      filter(stats.file.lg$both.outlier == TRUE) %>%
      select(gen.pos) %>%
      arrange(gen.pos) %>%
      mutate(dist.1 = c(NA,diff(gen.pos))) %>%
      mutate(dist.2 = c(diff(sort(gen.pos)),NA))
    
    nn.dist <- rep(NA, length(site.sample$genpos))
    for (k in 1:length(site.sample$gen.pos)){
      
      if(!is.na(site.sample$dist.1[k]) & !is.na(site.sample$dist.2[k])){
        nn.dist[k] <- min(c(site.sample$dist.1[k],site.sample$dist.2[k]))
      }else if(is.na(site.sample$dist.1[k])){
        nn.dist[k] <- site.sample$dist.2[k]
      } else if(is.na(site.sample$dist.2[k])){
        nn.dist[k] <- site.sample$dist.1[k]
      }
    }
    empirical.mean.nnd <- mean(nn.dist, na.rm = TRUE)
    
    #number of total loci
    n.sites <- stats.file.lg %>% filter(!is.na(gen.pos)) %>% select(gen.pos) %>% unlist %>% length
    
    nnd.stats[[j]] <- data.frame(lg = unique(stats.file.lg$lg),
                                 n.sites = n.sites,
                                 num.outliers = num.outliers,
                                 nnd.mean.null = null.mean, 
                                 nnd.sd.null = sd(null.mean.nnds, na.rm = TRUE),
                                 nnd.mean.emp = empirical.mean.nnd,
                                 nnd.emp.percentile = null.ecdf(empirical.mean.nnd),
                                 nnd.emp.zscore = (empirical.mean.nnd - null.mean)/sd(null.mean.nnds, na.rm = TRUE),
                                 nnd.emp.pvalue = 2*pnorm(-abs((empirical.mean.nnd - null.mean)/sd(null.mean.nnds, na.rm = TRUE))))
    }
  }
  return(do.call("rbind", nnd.stats))
}



################################################# MAIN FUNCTION

calculate_coeff_dispersion_stats_file <- function(stats.filename){
  
  #read in file
  
  print(paste0("Processing ",stats.filename,"..."))
  
  stats.file <- data.table(read.table(stats.filename, stringsAsFactors = FALSE, header=TRUE))
  
  setnames(stats.file, tolower(names(stats.file)))
  stats.file <- stats.file[,chrom.pos:=NULL]
  setnames(stats.file,1,"lg")
  stats.file$lg <- chrom.to.num(stats.file$lg)
  
  # call outliers
  stats.file <- stats.file %>%
    filter(!is.na(fst)) %>%
    filter(!is.infinite(fst)) %>%
    filter(fst > 0) %>%
    filter(!is.infinite(dxy)) %>%
    filter(!is.na(dxy)) %>%
    filter(dxy > 0) %>%
    mutate(fst.outlier = is.outlier(fst)) %>% 
    mutate(dxy.outlier = is.outlier(dxy)) %>% 
    mutate(both.outlier = dxy.outlier & fst.outlier)
  
  #add map distances
  stats.file <- add_map_distance(stats.file)
  
  ### calculate coefficients of dispersion for each lg
  
  dispersion.stats <- list()
  for (j in unique(stats.file$lg)){
    
    stats.file.lg <- stats.file %>%
      filter(lg == j)
    
    disp.lg.all <- calculate_dispersion(stats.file.lg )
    disp.lg.outl <- calculate_dispersion_outliers(stats.file.lg )
    
    dispersion.stats[[j]] <- data.frame(lg = j, disp.all = disp.lg.all, disp.out = disp.lg.outl)
  }
  
  # bind dispersion estimates into a df
  disp.df <- do.call("rbind", dispersion.stats)
  nnd.df <- calculate_nndist_all_lg(stats.file)
  
  #format the cluster df for output
  cluster.df <- left_join(nnd.df, disp.df, by = "lg")
  
  file.name.stripped <- sapply(strsplit(stats.filename, split = "/"), function(x)gsub(".txt","",x[length(x)]))
  file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
  cluster.df$pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
  cluster.df$ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
  cluster.df$pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
  cluster.df$ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
  cluster.df$geography <- file.name.split[3]
  cluster.df$ecology <- file.name.split[4]
  
  cluster.df <- cluster.df %>% select(pop1, ecotype1, pop2, ecotype2, geography, ecology, everything())
  
  return(cluster.df)
  
}

# burn and turn
cluster.master <- mclapply(stats.files, calculate_coeff_dispersion_stats_file, mc.cores = 12, mc.silent = FALSE) 
cluster.master <- do.call("rbind", cluster.master)

date.stamp <- paste("_", format(Sys.time(), "%Y-%m-%d"), sep="")
out.file.name <- file.path("analysis_ready", paste("clustering_master", date.stamp, ".txt", sep=""))
write.table(cluster.master,file = out.file.name, row.names = FALSE, quote = FALSE)


