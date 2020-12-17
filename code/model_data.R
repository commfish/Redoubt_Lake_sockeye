# Original code was from Steve Fleischman and Xinxian Zhang. code was modified
# and converted to JAGS by Sara Miller on 28 June 2017 and modified and
# maintained by Rich Brenner....lots of input from Ben Williams! Thanks Ben.
# to run model with only years 1999+, uncomment lines 17-20 and comment out lines 11-14
library(tidyverse)

# data----
brood <- read.csv("data/Redoubt_sockeye.csv", header = TRUE)
brood

# clean up data
#brood %>% 
#  dplyr::select(S=spawn, R=recruit) %>% 
#  mutate(lnRS = log(R/S), n()) -> sr

# clean up data and only run SR analysis from 1999 on for comparison
brood %>% 
  filter(year > 1998)%>%
  dplyr::select(S=spawn, R=recruit) %>% 
  mutate(lnRS = log(R/S), n()) -> sr

n <- nrow(sr) #calculates the number of years of data.
dat <- list(n = n, S = sr$S, lnRS = sr$lnRS)
