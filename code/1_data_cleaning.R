### Data cleaning
## Take output from calirometer and tidy up
# (1) remove contaminated runs
# (2) convert to correct units

## libraries needed
library(reshape2) # for data manipulation
library(ggplot2) # for plotting
library(plyr) # for data manipulation
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(MASS)
library(here)
library(tidyverse)
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

source("code/1_data_set1.R")
source("code/1_data_set2.R")
source("code/1_data_set3.R")
source("code/1_data_set4.R")
source("code/1_data_set5.R")
source("code/1_data_set6.R")
source("code/1_data_set7.R")
source("code/1_data_set8.R")
source("code/1_data_set9.R")
source("code/1_data_set10.R")
source("code/1_data_set11.R")
source("code/1_data_set12.R")
source("code/1_data_set13.R")
source("code/1_data_set14.R")

## outputs cleaned files into data
