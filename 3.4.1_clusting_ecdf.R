#ecdf test

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

is.outlier <- function(x){
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  return(x >=x95)
}

snp.file<-snp.file%>%
  group_by(study_com)%>%
  mutate(fst.outlier = is.outlier(fst)) %>% 
  mutate(dxy.outlier = is.outlier(dxy)) %>% 
  mutate(both.outlier = dxy.outlier & fst.outlier) %>% 
  ungroup


ecdf.dat.out <- snp.file %>%
  filter(both.outlier == TRUE) %>%
  group_by(study_com,lg) %>%
  do(outlier.ecdf = ecdf(.$gen.pos))

ecdf.dat.snp <- snp.file %>%
  group_by(study_com,lg) %>%
  do(outlier.ecdf = ecdf(.$gen.pos))

plot(ecdf.dat$outlier.ecdf[[4]])

plot(ecdf.dat.snp$outlier.ecdf[[4]])
