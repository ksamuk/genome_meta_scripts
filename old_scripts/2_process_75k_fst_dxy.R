####Match EVs to 75k stats files, fit models, write model results to file

#washy washy
rm(list=ls())

library("IRanges")
library("ggplot2")
library("dplyr")
library("robustbase")
library("devtools")
#install_github("karthik/wesanderson")
library("wesanderson")

# trusty chrom to num function
chrom.to.num <- function(x){
  x <- gsub("group", "", x)
  chrom.rom <- as.character(as.roman(c(1:21)))
  return(match(x, chrom.rom))
}

# call outlier function
is.outlier <- function(x){
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  return(x >=x95)
}

#ev dir location and file list
ev.dir <- file.path(getwd(),"evs")
ev.files <- list.files(ev.dir, pattern = "txt",full.names = TRUE)

#stats files location and list
stats.dir <- file.path(getwd(),"stats/75k_all")
stats.files <- list.files(stats.dir,"*sliding*", full.names = TRUE)

######################## match_evs function

match_evs <- function(stats.file.name){
  
  print(paste0("processing ", stats.file.name,"..."))
  #initialize data to be matched 
  stats.file <- read.table(stats.file.name, header=TRUE)
  if (length(stats.file[,1]) < 1){
    return(NA)
  }
  names(stats.file) <- c("lg","pos1","pos2","midpos","sites","var.sites","dxy","fst","hexp1","hexp2")
  stats.file$lg <- chrom.to.num(stats.file$lg)
  stats.file$lg <- as.numeric(stats.file$lg)
  matched.all <- data.frame(stats.file)
  ###looping through all the ev files
  for (i in 1:length(ev.files)){
    
    #read ev file
    ev <- read.table(file = ev.files[i], header=TRUE)
    
    #prep empty ev dataframe
    matched.evs <- data.frame()
    length(matched.evs) <- 4
    names(matched.evs) <- c("lg", "pos1", "pos2", "ev")
    
    #sanitize file of NAs and duped rows
    ev <- ev[complete.cases(ev), ]
    ev <- arrange(ev, lg, pos1)
    ev <- unique(ev)
    
    ###loop through lgs, matching evs as we god
    for (j in 1:max(stats.file$lg)){
    #for (j in 1:2){
      
      #subset ev and stat by lg
      ev.chr <- subset(ev,ev$lg == j)
      stat.chr <- subset(stats.file, stats.file$lg == j)
      stat.chr <- stat.chr[,1:3]
      
      #build IRanges objects for overlap finding
      ev.range <- IRanges(start = ev.chr$pos1, end = ev.chr$pos2)
      #minimal set. note that this means overlapping evs end up returning the same value, rather than say, their average
      #(this shouldn't matter for most types of evs)
      ev.range <- reduce(ev.range, min.gapwidth = 0L)
      stat.range <- IRanges(start = stat.chr$pos1, end = stat.chr$pos2)
      
      #find ovelaps amd build an "overlap df"
      overlap <- findOverlaps(stat.range, ev.range, select = "all")
      overlap.df <- data.frame(lg = stat.chr[queryHits(overlap),]$lg,
                             pos1 = stat.chr[queryHits(overlap),]$pos1,
                             pos2 = stat.chr[queryHits(overlap),]$pos2,
                             ev.start = start(ev.range[subjectHits(overlap)]),
                             ev.end = end(ev.range[subjectHits(overlap)]))
      
      #truncate overlaps that start before window
      overlap.df$ev.start[overlap.df$ev.start < overlap.df$pos1] <- overlap.df$pos1[overlap.df$ev.start < overlap.df$pos1]
      
      #truncate overlaps that extend past end of window
      overlap.df$ev.end[overlap.df$ev.end > overlap.df$pos2] <- overlap.df$pos2[overlap.df$ev.end > overlap.df$pos2]
      
      #calc widths
      overlap.df$width <- (overlap.df$ev.end - overlap.df$ev.start) + 1
      
      overlap.df$ev <- ev.chr[subjectHits(overlap), 4]
      overlap.df <- unique(overlap.df)
      
      #calculate the proportion of total overlap
      prop.overlap <- overlap.df %>%
        group_by(pos1) %>%
        summarise(total.overlap = sum(width))
      
      #add it back into the overlap dataframe
      prop.overlap <- data.frame(prop.overlap)
      overlap.df$total.overlap <- prop.overlap$total.overlap[match(overlap.df$pos1, prop.overlap$pos1)]
      
      #the weighted contribution this to the total (deals with multiple ev hits to a single window)
      overlap.df$ev.prop <- overlap.df$ev * (overlap.df$width / overlap.df$total.overlap)
      
      #calculate the weighted evs for each window
      weighted.evs <- overlap.df %>%
        group_by(pos1) %>%
        summarise(ev = sum(ev.prop))
      weighted.evs <- data.frame(weighted.evs)
      #add back into df
      overlap.df$ev.avg <- weighted.evs$ev[match(overlap.df$pos1, weighted.evs$pos1)]
      
      ####DIVERGES FROM SNP BASED ANALYSIS
      #HACK TO ACCOMODATE WINDOWS re-read stats file
      stat.chr <- subset(stats.file, stats.file$lg==j)
      stat.chr$ev <- overlap.df$ev.avg[match(stat.chr$pos1, overlap.df$pos1)]
      matched.evs.chr <- stat.chr
      
      #rbind chromos
      matched.evs <- rbind(matched.evs, matched.evs.chr)
      
    }
    ###end lg loop
    
    #attach real name of ev and cbind to stats file
    names(matched.evs)[11] <- sapply(strsplit(ev.files[i], split="/"), function(x)gsub(".txt","",x[length(x)]))
    matched.all <- suppressMessages(left_join(matched.all, matched.evs))
    
  }

  # call outliers
  matched.all <- matched.all %>%
    mutate(fst.outlier = is.outlier(fst))%>%
    mutate(dxy.outlier = is.outlier(dxy))%>%
    mutate(both.outlier = dxy.outlier == TRUE & fst.outlier == TRUE)%>%
    mutate(hs = (hexp1+hexp2)/2)
  
  #filter matched.all
  
  # Recombination distances >25cM
  #matched.all$recomb_rate[matched.all$recomb_rate >= 25] <- NA
  # dS
  matched.all$ds[matched.all$ds >= 3] <- NA
  # gene_count
  matched.all$gene_count[matched.all$gene_count >= 15] <- NA
  # no lg 19
  matched.all <- matched.all %>%
    filter(lg!=19)
  #hs
  matched.all <- matched.all %>%
    mutate(hs = (hexp1+hexp2)/2)
  # fit model
  
  model.out <- glm(as.numeric(both.outlier) ~ recomb_rate + ds + gene_count,
                   na.action = "na.omit",
                   family = binomial,
                   data = matched.all)
  
  coeffs <- as.list(model.out$coefficients)
  names(coeffs)[1] <- "intercept"
  
  # parse file name
  file.name.stripped <- sapply(strsplit(stats.file.name, split = "/"), function(x)gsub(".txt","",x[length(x)]))
  file.name.split <- strsplit(file.name.stripped, split = "[.]") %>% unlist
  pop1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype1 <- file.name.split[1] %>% strsplit(split="_") %>% unlist %>% .[2]
  pop2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[1]
  ecotype2 <- file.name.split[2] %>% strsplit(split="_") %>% unlist %>% .[2]
  geography <- file.name.split[3]
  ecology <- file.name.split[4]
  
  row.out <- data.frame(pop1, ecotype1, pop2, ecotype2, geography, ecology, coeffs)
  
  return(row.out)
  
}

######################## end match_evs function

# apply the above function to the list of stats.files
coeff.df <- lapply(stats.files, match_evs)

#bind into a data frame
coeff.dat <- do.call("rbind", coeff.df)

#remove nas
coeff.dat <- coeff.dat[!is.na(coeff.dat$recomb_rate),]
coeff.dat <- coeff.dat %>%
  mutate(group = paste0(geography,"_",ecology))

## add in region data (for looser geography)
region.dat <- read.table(file = "meta_data/population_regions.txt", header = TRUE, stringsAsFactors = FALSE)

# make short pop codes
region.dat$pop.code <- strsplit(region.dat$pop.code, split = "_") %>% lapply(function(x)return(x[1])) %>% unlist

# associate pops in coeff dat with regions in region dat
region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop1","reg1")
coeff.dat$reg1 <- region.sub$reg1[match(coeff.dat$pop1, region.sub$pop1)]

region.sub <- region.dat[,c(1,4)]
names(region.sub) <- c("pop2","reg2")
coeff.dat$reg2 <- region.sub$reg2[match(coeff.dat$pop2, region.sub$pop2)]

# make new geographic categories
coeff.dat$geography2 <- ifelse(coeff.dat$reg1==coeff.dat$reg2, "para", "allo")

#make new groups :o
coeff.dat$group2 <- paste0(coeff.dat$geography2,"_",coeff.dat$ecology)

# write to file
write.table(coeff.dat, file = "analysis_ready/75k_stats_model_dxy_fits.txt", row.names = FALSE, quote = FALSE)


