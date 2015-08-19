####create data table from (pre-processed) codeml output
####outputs gene ids and ds values for gacu ONLY
####KS mar-3-2015

library(ape)

#home dir set up
#home.dir<-"E:/Genome Meta Analysis/ev_prep_scripts/paml_analysis"
#home.dir<-"~/Documents/Science/Projects/Ph.D./Genome Meta Analysis/ev_prep_scripts/paml_analysis"

#contained pre-processed paml output (the dn and ds trees)
paml.output.dir <- file.path("ev_prep_scripts/paml/output")

#dat file list
file.list <- list.files(paml.output.dir, full.names = TRUE)

#set up empty output dataframe
gene.id <- vector(mode = "character", length = length(file.list))
ds <- vector(mode = "numeric", length = length(file.list))
gacu.ds <- data.frame(gene.id, ds, stringsAsFactors = FALSE)

for (i in 1:length(file.list)){
  #read in file
  file.con <- file(file.list[i])
  file.lines <- readLines(file.con)
  close(file.con)
  if (!is.na(file.lines[2])){  
    # make the tree
    ds.tree <- read.tree(text = file.lines[3])
    
    # distance matrix for branches
    dist.mat <- cophenetic.phylo(ds.tree)
    
    # the gacu row in the distance matrix
    gacu.row <- row.names(dist.mat) %>% grep("ENSGACP",.) %>% dist.mat[.,]
    
    # only retain estimates for pformosa, oniloticus, xmaculatus
    gacu.row <- gacu.row %>% names %>% grep("pformosa|oniloticus|xmaculatus",.) %>% gacu.row[.]
    
    # the ds estimate 
    gacu.ds$ds[i] <- mean(gacu.row, na.rm = TRUE)
    
    #get gene name and find it in the ds tree
    gene.name <- strsplit(file.list[i], split=".cml.txt") %>% gsub(paml.output.dir,"",.) %>% gsub("/","",.)
    
    #pull out ds value for gacu
    gacu.ds$gene.id[i] <- gene.name
    
    #tracker
    if (i %% 1000 == 0){
      print(i)
    }
    
  }
  if (is.na(file.lines[1])){
    gacu.ds$gene.id[i] <- NA
    gacu.ds$ds[i] <- NA
  }
  
}

# add in locations

gene.dat <- read.table(file = "evs/additional/gene_id.txt", stringsAsFactors = FALSE, header = TRUE)
gacu.ds$ensembl_gene_id <- gacu.ds$gene.id %>% gsub("P","G",.)

#output

test <- left_join(gacu.ds, gene.dat)


setwd(home.dir)
write.table(gacu.ds, file="ds_estimates_gacu.txt")

