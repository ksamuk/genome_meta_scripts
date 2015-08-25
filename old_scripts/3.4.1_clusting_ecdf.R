#ecdf test


library(dplyr)
library(ggplot2)
##### PROOF OF CONCEPT
test.dat1 <- rnorm(100)
test.dat2 <- runif(100)
ecdf1 <- ecdf(test.dat1)
ecdf2 <- ecdf(test.dat2)
plot(ecdf2)
ks.test(ecdf1(test.dat1), ecdf2(test.dat2))
ks.test(test.dat1, test.dat2)
##### PROOF OF CONCEPT

snp.file <- list.files(file.path("stats/snp_filtered"),pattern="genpos", full.names = TRUE)
snp.file <- read.table(snp.file, header = TRUE, stringsAsFactors = FALSE)

snp.file <- list.files(file.path("stats/snp_filtered"),pattern="genpos", full.names = TRUE)
snp.file <- read.table(snp.file, header = TRUE, stringsAsFactors = FALSE)

snp.file[!is.na(snp.file$fst) & snp.file$fst<0, ]$fst <- NA

is.outlier <- function(x){
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  return(x >=x95)
}

outlier.quantile <- function(x){
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  return(x95)
}

snp.file<-snp.file%>%
  group_by(study_com)%>%
  mutate(fst.outlier = is.outlier(fst)) %>% 
  mutate(cutoff = outlier.quantile(fst)) %>% 
  mutate(dxy.outlier = is.outlier(dxy)) %>% 
  mutate(both.outlier = dxy.outlier & fst.outlier) %>% 
  ungroup

cutoffs<-snp.file%>%
  select(study_com, cutoff) %>%
  distinct

snp.file %>%
  ggplot(aes(x = fst))+
    geom_histogram()+
    facet_wrap(~study_com)
  

ecdf.dat.out <- snp.file %>%
  filter(both.outlier == TRUE) %>%
  group_by(study_com,lg) %>%
  do(outlier.ecdf = ecdf(.$gen.pos))

ecdf.dat.snp <- snp.file %>%
  group_by(study_com,lg) %>%
  do(outlier.ecdf = ecdf(.$gen.pos))

plot(ecdf.dat$outlier.ecdf[[4]])

plot(ecdf.dat.snp$outlier.ecdf[[4]])
