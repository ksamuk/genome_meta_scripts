# process a folder of stats files, compute clustering, perform permutation test
# ks august 2015

# libraries

library("dplyr")
library("IRanges")

# folders

stats.folder <- "stats/stats_all"
stats.files <- list.files(stats.folder, full.names = TRUE)

# functions

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
    
    print(paste0("processing lg ",i,"..."))
    
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

################################################# MAIN FUNCTION

calculate_coeff_dispersion_stats_file <- function(stats.filename){
  
  #read in file
  
  stats.file <- data.table(read.table(stats.filename, stringsAsFactors = FALSE))
  
  # call outliers
  stats.file <- stats.file %>%
    mutate(fst.outlier = is.outlier(fst)) %>% 
    mutate(dxy.outlier = is.outlier(dxy)) %>% 
    mutate(both.outlier = dxy.outlier & fst.outlier)
  
  #add map distances
  
  stats.file <- add_map_distance(stats.file)
  
  #calc dispersion of all loci
  
  stats.file
  
  #calc dispersion of outliers
  
}



