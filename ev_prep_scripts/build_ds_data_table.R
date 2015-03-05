####create data table from (pre-processed) codeml output
####

library(ap)ea

home.dir<-"E:/Genome Meta Analysis/ev_prep_scripts/paml_analysis"
setwd(home.dir)

file.con<-file("ENSGACP00000014140.cml.txt")
file.lines<-readLines(file.con)

ds.tree<-read.tree(text=file.lines[2])

plot(ds.tree)

gacu.edge<-which.edge(ds.tree, c("ENSGACP00000014140"))

ds.tree$edge.length[gacu.edge]
ds.

branchlist <- 
