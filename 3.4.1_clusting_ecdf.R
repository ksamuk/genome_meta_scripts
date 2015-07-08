#ecdf test


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


outlier.ecdf <- snp.file %>%
  filter(both.outlier == TRUE) %>%
  group_by(study_com,lg) %>%
  do(outlier.ecdf = ecdf(.$gen.pos))


