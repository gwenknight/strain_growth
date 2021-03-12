### Analysis of growth data

###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
## libraries needed

library(reshape2) # for data manipulation
library(ggplot2) # for plotting
library(plyr) # for data manipulation
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(MASS)
library(tidyverse)
library(RColorBrewer)
library(here)
theme_set(theme_bw(base_size=6)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## Load in grofit functions if package no longer working
source("code/grofit_functions.R") # have copied those we use into this R script

## Load in functions for this analysis
source("code/functions_for_heat_curves.R")

## where is the data? These are the outputs from 1_data_set#.R: standardised all variable names etc in here
ddm1 <- as.data.table(read.csv("data/ddm_set1.csv")[,-1])
ddm2 <- as.data.table(read.csv("data/ddm_set2.csv")[,-1])
ddm3 <- as.data.table(read.csv("data/ddm_set3.csv")[,-1])
ddm4 <- as.data.table(read.csv("data/ddm_set4.csv")[,-1])
ddm5 <- as.data.table(read.csv("data/ddm_set5.csv")[,-1])
ddm6 <- as.data.table(read.csv("data/ddm_set6.csv")[,-1])
ddm7 <- as.data.table(read.csv("data/ddm_set7.csv")[,-1])
ddm8 <- as.data.table(read.csv("data/ddm_set8.csv")[,-1])
ddm9 <- as.data.table(read.csv("data/ddm_set9.csv")[,-1])
ddm10 <- as.data.table(read.csv("data/ddm_set10.csv")[,-1])
ddm11 <- as.data.table(read.csv("data/ddm_set11.csv")[,-1])
ddm12 <- as.data.table(read.csv("data/ddm_set12.csv")[,-1])
ddm13 <- as.data.table(read.csv("data/ddm_set13.csv")[,-1])
ddm14 <- as.data.table(read.csv("data/ddm_set14.csv")[,-1])
ddm <- as.data.frame(rbind(ddm1,ddm2,ddm3,ddm4,ddm5,ddm6,ddm7,ddm8,ddm9,ddm10,ddm11,ddm12,ddm13,ddm14) )

length(unique(ddm$strain)) # 98 in sets 1-13


#### NAME code: how to label output depending on strains being analysed
#name_code <- "1_2_6_7_" # old example
name_code <- "cut_"

###******* PLOTTING *************#######################################################################################################################################
# To look at raw data: commented out as only do once at the start

# w0 <- which(ddm$drytime == 0)
# w24 <- which(ddm$drytime == 24)
# w168 <- which(ddm$drytime == 168)
# 
# ggplot(ddm[w0,], aes(x=Time,y=value, group = interaction(rep, inoc), col = factor(inoc))) +
#   facet_wrap(~strain, scales="free") +
#   geom_line() + ggtitle("Raw data: t0")
# ggsave("output/all_t0.pdf")
# 
# ggplot(ddm[w24,], aes(x=Time,y=value_J, group = interaction(rep, inoc), col = factor(inoc))) +
#   facet_wrap(~strain) +
#   geom_line() + ggtitle("Raw data: t24")
# ggsave("output/all_t24.pdf")
# 
# ggplot(ddm[w168,], aes(x=Time,y=value_J, group = interaction(rep, inoc), col = factor(inoc))) +
#   facet_wrap(~strain) +
#   geom_line() + ggtitle("Raw data: t168")
# ggsave("output/all_t168.pdf")

###******** UNITS / DATA *************#######################################################################################################################################
ddm$value_J = ddm$value

###******* MODEL FITTING *************#######################################################################################################################################
## Fit growth curves to each set of data to determine the underlying parameters 

## What are the strains?
u <- as.character(unique(ddm$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(q)*3, 19); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
drying_times <- c(0,168)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] # which replicates for this strain (all have different names)
  for(ii in 1:length(r)){ # for each replicate
    for(kk in c(1,2)){ #each of the experimental conditions
      for(ll in 1:length(q)){ #each of the inocula
        
        strain <- u[jj];
        replicate <- r[ii]
        condition <- drying_times[kk]
        inocl <- q[ll]
        
        wi <- intersect(which(ddm$strain == strain),which(ddm$rep == replicate)) # if fit to each replicate
        wj <- intersect(wi, which(ddm$drytime == condition))
        w <- intersect(wj, which(ddm$inoc == as.numeric(inocl)))
        
        if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
          data1 <- ddm[w,] # Grab data
          
          print(c(strain, replicate, condition, inocl)) # output so can track how it is working
          p <- cut_extract(data1, "Time", "value_J", paste(strain, replicate, condition, inocl,sep="_")) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
          
          ## Required parameters
          
          param[index,] <- c(strain, replicate, condition, inocl, p$param)
          index <- index + 1 # counting for storing matrix - next row each iteration
        }
        
      }
    }
  }
}

# The original set 
param_orig <- param

## Fitted parameters
param <- as.data.frame(param)
colnames(param) <- c("strain_name","rep","drytime","inocl",
                     "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v",
                     "cut_exp", "timepeak", "valpeak")

w<-which(param$lag!=0); param <- param[w,] # remove 0 values

dim(param)

#### CHECK FOR ODD BEHAVIOUR
param$any_odd <- as.numeric(param$odd_peaks) + as.numeric(param$odd_width) + as.numeric(param$odd_shoulder)
param$odd_type <- "0"
param$odd_type_db <- "0"
param[which(param$odd_peaks > 0),"odd_type"] <- paste0(param[which(param$odd_peaks > 0),"odd_type"], "1")
param[which(param$odd_width > 0),"odd_type"] <- paste0(param[which(param$odd_width > 0),"odd_type"], "2")
param[which(param$odd_shoulder > 0),"odd_type"] <- paste0(param[which(param$odd_shoulder > 0),"odd_type"],"3")
param[which(param$odd_double > 0),"odd_type_db"] <- paste0(param[which(param$odd_double > 0),"odd_type"],"4")
unique(param$odd_type)
unique(param$odd_type_db)

##### Save output 
write.csv(param,paste0("output/",name_code,"all_model_fit_params.csv"))


######****** ODD behaviour & CUT POINT into main timeseries data ******#################
no_odd <- c()
ddm$odd_type <- "0"
ddm$odd_type_db <- "0"
ddm$shoulder_point_t <- 0
ddm$shoulder_point_v <- 0
ddm$shoulder_cut <- 0

#### store odd type in main timeseries dataframe
for(jj in 1:length(u)){ # for each strain
  
  pp <- param %>% dplyr::filter(strain_name == u[jj])
  ii <- which(ddm$strain == u[jj])
  print(u[jj])
  
  # Look at odd ones
  w<-which(pp$any_odd > 0)
  if(length(w > 0)){
    for(ww in 1:length(w)){
      w1 <- which(as.numeric(unlist(ddm[ii,"rep"])) == pp[w[ww],c("rep")])
      w2 <- intersect(w1,which(unlist(ddm[ii,"drytime"]) == pp[w[ww],"drytime"]))
      w3 <- intersect(w2,which(unlist(ddm[ii,'inoc']) == pp[w[ww],"inocl"]))
      
      oddv <- ii[w3]
      
      ddm[oddv,"odd_type"] <- pp[w[ww],"odd_type"] # label as odd across timeseries
      ddm[oddv,"odd_type_db"] <- pp[w[ww],"odd_type_db"] # label as odd across timeseries
      
    }
  }
  
  # Move over cut point or time to peak 
  for(ww in 1:length(pp[,1])){ # for every row
    w1 <- which(as.numeric(unlist(ddm[ii,"rep"])) == pp[ww,c("rep")])
    w2 <- intersect(w1,which(unlist(ddm[ii,"drytime"]) == pp[ww,"drytime"]))
    w3 <- intersect(w2,which(unlist(ddm[ii,'inoc']) == pp[ww,"inocl"]))
    
    oddv <- ii[w3] ## this 
    
    ### If have a shoulder then add in new point up to this cut 
    if(pp[ww,"shoulder_point_t"]>0){
      ddm[oddv,c("shoulder_point_t","shoulder_point_v","shoulder_cut")] <- c(pp[ww,c("timepeak","valpeak")],1) # label as odd across timeseries
    }else{
      ddm[oddv,c("shoulder_point_t","shoulder_point_v","shoulder_cut")] <- c(pp[ww,c("timepeak", "valpeak")],0) # add in time to peak as cut point
    }
    
    # If have any odd behaviour then do the extra cut fit to check no earlier peaks (e.g. 11214 11.2)
    if(pp[ww,"any_odd"] > 0){
      ddm[oddv,c("shoulder_cut")] <- 1 
    }
    
  }
  
  
  
  
}

## Making odd_type a factor 
ddm$odd_type <- factor(ddm$odd_type, levels = c("0","1","2","3","12","13","23","123"))
ddm$odd_type_db <- factor(ddm$odd_type_db, levels = c("0","1","2","3","4","12","13","14","134","34","23","24","123","124","234","1234"))

#### Label MACOTRA vs. not
ddm$source <- "Macotra"
w<-which(ddm$strain %in% c("Newman","RWW12","SA3297","SA2704","RWW146","SAC042W", "Mu50", "M116"))
ddm[w,"source"] <- "Other"

##### Save output 
write.csv(ddm,paste0("output/",name_code,"all_time_series_fit_params.csv"))

