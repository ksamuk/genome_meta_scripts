# calculate genetic distances for snp file

library(data.table)
library(dplyr)

# read in snp file
snp.file <- list.files(file.path("stats/snp_filtered"), pattern=".txt$", full.names = TRUE)
snp.file <- fread(snp.file)


snp.file$gen.pos <- NA

# read in genetic map data
map.file <- list.files(file.path("ev_prep_scripts"), pattern=".txt$", full.names = TRUE)
map.file <- read.table(map.file,header=TRUE)
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

snp.file.out <- list.files(file.path("stats/snp_filtered"),pattern=".txt$", full.names = TRUE)
write.table(snp.file, file = snp.file.out, row.names = FALSE)
