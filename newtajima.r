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
  x5 <- quantile(x, na.rm = TRUE, probs = 0.05)[1]
  x95 <- quantile(x, na.rm = TRUE, probs = 0.95)[1]
  
  return(x <= x5 | x >=x95)
}

taj.dat <- taj.dat %>%
  mutate(id = paste0(study, location, ecotype)) %>%
  filter(taj.dat$n.sites > 2) %>%
  group_by(id) %>%
  mutate(outlier = is_taj_outlier(tajd)) %>%
  ungroup() %>%
  mutate(taj.class =  sapply(tajd,make_taj_class))


taj.dat %>%
  filter(study == "benlim") %>%
  #filter(location == "Pax") %>%
  #filter(ecotype == "B") %>%
  filter(outlier == TRUE) %>%
  ggplot(aes(y = tajd, x = pos1, color = taj.class)) +
  geom_point() +
  facet_wrap(~lg)