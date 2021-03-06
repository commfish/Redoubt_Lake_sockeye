# R CODE TO CREATE PROFILES DATA AND FIGURES (Horsetail plots, age comp, maturity at age,
# point estimate plots, yield plots)
# created by Ben Williams (Ben.Williams@alaska.gov);Nov 3, 2016; 2017-7-4
# Changes made by Sara Miller (Sara.Miller@alaska.gov); April 2017

LowerB <- 10000 #lower bound of recommended escapement goal range
UpperB <- 25000 #upper bound of recommended escapement goal range
SMSY <- 18174.978 #Lambert W version of SMSY from file: AR_quantiles_lambert
lnalpha.c <- 2.58690031604307
lnalpha <- 2.07841039608424
beta <- 0.0000449797444705782  


#load----
library(plyr)
library(tidyverse)
library(reshape2)
library(grid)
library(gsl)
library(scales)
library(backports)
devtools::install_github("commfish/fngr")
library(fngr)
library(FField)
#font_import() # only need run once
extrafont::font_import()
windowsFonts(Times=windowsFont("TT Times New Roman"))
theme_set(theme_report(base_size = 14))
source('code/functions.r')

if(!dir.exists(file.path("output", "base_model", "processed"))){dir.create(file.path("output", "base_model", "processed"))}

# data----
spawnrecruitdat <- read.csv("data/Redoubt_sockeye.csv") #Load Data File (make sure this file is updated)
coda <- read.csv("output/base_model/coda.csv") 

# data clean----
# profile spawnrecruitdat
coda %>% 
  mutate(S.eq.c = lnalpha.c/beta, 
         S.msy.c = (1-lambert_W0(exp(1-lnalpha.c)))/beta, #Lambert W
         R.msy.c = S.msy.c*exp(lnalpha.c-beta*S.msy.c), 
         MSY.c = R.msy.c-S.msy.c, 
         Rmax = exp(lnalpha)*(1/beta)*exp(-1)) -> coda
# analysis----
# create function for probability profiles and figures
profile(i=10, z=50, xa.start=0, xa.end=8000,lnalpha.c, beta) #can change i,z, xa.start, xa.end
QM <- read.csv("output/base_model/processed/QM.csv")
CI <- read.csv("output/base_model/processed/CI.csv")
num <- nrow(QM)
QM %>%
  dplyr::select(c(escapement)) -> x
coda %>%
  dplyr::select(c(lnalpha.c, beta)) %>%
  filter(row_number()==1:50) %>%
  slice(rep(1:n(), each = num)) %>%
  cbind(., x) %>% #lnalpha.c, beta, and S 
  mutate(recruitment = escapement*exp(lnalpha.c-beta*escapement),
         variable = rep(1:50,each=num),
         year = "") -> dataset

QM %>%
  dplyr::select(c(escapement)) %>%
  mutate (lnalpha.c = lnalpha.c,
          beta = beta,
          recruitment = escapement*exp(lnalpha.c-beta*escapement),
          variable = 51,
          year = "")%>%
  rbind(., dataset) -> dataset

spawnrecruitdat %>%
  mutate(escapement = spawn,
         lnalpha.c = NA,
         beta = NA,
         recruitment =recruit,
         variable = 52) %>%
  dplyr::select(c(escapement, lnalpha.c, beta, recruitment, variable, year)) %>%
  rbind(., dataset) %>%
  arrange(variable, escapement) -> dataset

ggplot(data=dataset, aes(x=escapement, y=recruitment, group=variable))+
  geom_line(data=subset(dataset,dataset$variable<52),linetype="solid", size=0.5, color="grey80")+
  scale_y_continuous(labels = comma,breaks = seq(0, 250000, 25000), limits = c(0, 250000))+
  scale_x_continuous(labels= comma,breaks = seq(0, 150000, 25000), limits = c(0, 150000))+
  ylab("Recruits (R)")+xlab("Spawners (S)")+
  geom_line(data=dataset, aes(x=escapement, y=escapement, group=1),linetype="solid", size=1)+#replacement line
  geom_line(data=subset(dataset,variable==51),colour = "black", lty=2, size=2)+
  geom_point(data=subset(dataset,variable==52),colour = "black", pch=16, size=1)+
  geom_text(size=3, data=dataset, aes(x=escapement, y=recruitment, group=52, label=year,family="Times", 
                                      hjust = -0.1, vjust= -0.4)) 
ggsave("output/base_model/processed/horsetail.png", dpi = 500, height = 6, width = 8, units = "in")

#Alternative horsetail plot
dataset %>%
  filter (variable %in% c(52)) %>%
  mutate(escapement1 = escapement) %>%
  dplyr::select(-c(lnalpha.c, beta, escapement)) %>%
  mutate(escapement = 'NA',
         Median = 'NA',
         q95 = 'NA',
         q90 ='NA',
         q10 ='NA',
         q5 = 'NA') -> dataset
CI %>%
  mutate(year = 'NA',
         variable = 51,
         recruitment = 'NA',
         escapement1 ='NA')%>%
  rbind(., dataset) %>%
  mutate_if(is.character, as.numeric) -> dataset1

x.fact <- 100/max(dataset1$escapement) 
y.fact <- 100/max(dataset1$recruitment1)
coords <- FFieldPtRep(coords = cbind(dataset1$escapement1 * x.fact, dataset1$recruitment * y.fact), rep.fact = 40)
x.t <- coords$x/x.fact 
y.t <- coords$y/y.fact 

ggplot(data=dataset1, aes(x=escapement, y=Median, group=variable)) + 
  geom_line(size=1, lty=2, group=51) +
  geom_ribbon(aes(ymin = q5, ymax = q95, group=51), alpha=.08) +
  geom_ribbon(aes(ymin = q10, ymax = q90, group=51), alpha=.08) +
  xlab('Spawners (S)') +
  ylab('Recruits (R)') +
  scale_y_continuous(labels = comma,breaks = seq(0, 250000, 25000), limits = c(0, 250000)) +
  scale_x_continuous(labels = comma,breaks = seq(0, 150000, 25000), limits = c(0, 150000)) +
  geom_line(aes(x=escapement, y=escapement, group=51),linetype="solid", size=1) +
  geom_point(data=dataset1, aes(x=x.t, y=y.t, group=52),pch=16, size=1) +
  geom_point(data=dataset1, aes(x=escapement1, y=recruitment, group=52),pch=16, size=1) +
  geom_text(size=3, data=dataset1, aes(x=escapement1, y=recruitment, group=52, label=year,family="Times", 
                                       hjust = -0.1, vjust= -0.4)) 
ggsave("output/base_model/processed/horsetail2.png", dpi = 500, height = 6, width = 8, units = "in")

