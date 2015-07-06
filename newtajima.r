rm(list=ls())

library(ggplot2)
library(dplyr)

taj.dat <- read.table(file = "analysis_ready/tajimasD.combined.D")

# classify windows by tajima's d value (for colour)

make_taj_class <- function(x){
  
  if (x == 0){
    return("zero")
  }else if(x > 0){
    return ("positive")
  }else if(x < 0){
    return("negative")
  }
  
}

# classify windows as outliers (pos or neg)
is_taj_outlier <- function(x){
  if (x >= quantile(x, na.rm = TRUE, probs = 0.95)[1]){
    return(TRUE)
  }else if(x <= quantile(x, na.rm = TRUE, probs = 0.05)[1]){
    return(TRUE)
  }
}

taj.dat <- taj.dat %>%
  mutate(id = paste0(study, location, ecotype)) %>%
  filter(taj.dat$n.sites != 0) %>%
  group_by(id) %>%
  mutate(outlier = is_taj_outlier(tajd)) %>%
  mutate(taj.class = make_taj_class(tajd))


taj.dat %>%
   %>%
  filter(study == "benlim") %>%
  filter(location == "Pax") %>%
  filter(ecotype == "B") %>%
  filter(outlier == TRUE) %>%
  ggplot(aes(y = tajd, x = pos1, color = taj.class)) +
  geom_point() +
  facet_wrap(~lg)