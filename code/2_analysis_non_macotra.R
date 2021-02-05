#### Analysis of 8 strains that are not in MACOTRA

### Explored these to look at linear relationship 

#################**************** (1) Libraries, data and code needed *******************###############
library(tidyverse) 
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(ggplot2)
library(patchwork) # for combining plots
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

#####*************************** READ IN DATA *******************###############
ddm_orig <- read.csv("output/cut_all_time_series_fit_params.csv")[,-1]

ddm <- ddm_orig %>% filter(source == "Other")

param <- read.csv("output/cut_all_model_fit_params.csv")[,-1]

strains <- unique(ddm$strain)

ggplot(param %>% filter(strain_name %in% strains), aes(x= inocl, y = timepeak, group = strain_name)) + 
  geom_point(aes(col = strain_name)) + 
  facet_wrap(~strain_name, nrow = 2) + theme(legend.position = "none") + 
  scale_x_continuous(breaks = seq(2:7), labels = function(x) parse(text=paste("10^",x)),"Inoculum", limits = c(1.5,6.5)) + 
  scale_y_continuous(expression(paste("Time to first peak (", italic(tmax),")")))
ggsave("plots/time_peak_vs_inoculum_strain_col.pdf", width = 12)


ggplot(param %>% filter(strain_name %in% strains), aes(x = inocl, y = timepeak, group = inocl)) + 
  geom_boxplot() + 
  scale_x_continuous(breaks = seq(2:7), labels = function(x) parse(text=paste("10^",x)),"Inoculum", limits = c(1.5,6.5)) + 
  scale_y_continuous(expression(paste("Time to first peak (", italic(tmax),")")))
ggsave("plots/time_peak_vs_inoculum.pdf")
