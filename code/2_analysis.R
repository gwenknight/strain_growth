### Analysis of sets 1, 2 & 6

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

setwd(here:here())

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
ddm <- as.data.frame(rbind(ddm1,ddm2,ddm3,ddm4,ddm5,ddm6,ddm7,ddm8,ddm9,ddm10,ddm11,ddm12,ddm13) )

length(unique(ddm$strain)) # 98 in 1-13

# Add in cumulative 
#ddm <- ddm %>% group_by(rep, drytime, strain, inoc) %>% mutate(csum = cumsum(value)) %>% ungroup()
#ggplot(ddm, aes(x=Time, y = csum, group = interaction(inoc, strain, rep, drytime))) + geom_line()


#### 1_2_6_7: for just a subset of analysis
#ddm1 <- as.data.table(read.csv("data/ddm_set1.csv")[,-1])
#ddm2 <- as.data.table(read.csv("data/ddm_set2.csv")[,-1])
#ddm6 <- as.data.table(read.csv("data/ddm_set6.csv")[,-1])
#ddm7 <- as.data.table(read.csv("data/ddm_set7.csv")[,-1])

#ddm <- as.data.frame(rbind(ddm1,ddm2,ddm6,ddm7))
#ddm <- ddm6

#### NAME code: how to label output dependening on strains being analysed
#name_code <- "1_2_6_7_"
name_code <- "all_(1_13)_"

###******* PLOTTING *************#######################################################################################################################################
# To look at raw data

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
## Also those with odd behaviours are flagged here

## What are the strains?
u <- as.character(unique(ddm$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

## Run thru each experiment (columns in original data) for each strain
# Fit separately to each replicate as otherwise don't have enough data for later predictions
# still having up to 3 experimental drytimes despite set2 not having 24hr data
# Where the parameters for each strain are stored
param_n <- matrix(0, length(u)*length(r)*length(q)*3, 4); # number of strains x number of replicates x number of experimental conditions
param <- matrix(0, length(u)*length(r)*length(q)*3, 12); # number of strains x number of replicates x number of experimental conditions
index <- 1 # for counting 
max <- c() # for calibration
p_double_curves <- c() # for storing the double curve plot data
drying_times <- c(0,24,168)

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1]
  for(ii in 1:length(r)){ # for each replicate: fit to all the data, not just each replicate
    for(kk in c(1,3)){ #each of the three experimental conditions (0, 24, 168): most just 0 168 now
      for(ll in 1:length(q)){ #each of the inocula
        
        print(c(jj,ii,kk,ll))
        
        p <- fit_growth_curve(u[jj], r[ii], drying_times[kk],q[ll], ddm, 0,90) # parameters for growth curve
        
        ## Required parameters
        if(length(p$param_n) > 0){ # IF DATA then store
          param_n[index,] <- p$param_n
          param[index,] <- p$param
          max <- c(max,p$max_level)
          index <- index + 1 # counting for storing matrix - next row each iteration
          
          if(p$para_dbl[1] > 0){ # IF DOUBLE CURVE fitted then store output
            p_double_curves <- rbind(p_double_curves, cbind(p$plot_dbl, u[jj], r[ii], drying_times[kk],q[ll]))
          }
        }
        
      }
    }
  }
}

# The original set 
param_orig <- param

## Fitted parameters
param <- as.data.frame(param)
param_n <- as.data.frame(param_n)
colnames(param_n) <- c("strain_name","rep","drytime","inocl")
colnames(param) <- c("t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v")

w<-which(param_n$inocl!=0); param_n <- param_n[w,]
w<-which(param$lag!=0); param <- param[w,]

## put in names of strains
param$strain_name <- param_n$strain_name
param$rep <- param_n$rep
param$drytime <- param_n$drytime
param$inocl <- param_n$inocl

dim(param) # 523 now (sets 1&2&6)


#### CHECK FOR ODD BEHAVIOUR
param$any_odd <- param$odd_peaks + param$odd_width + param$odd_shoulder 
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

colnames(p_double_curves) <- c("time","value_J","pred","normal_curve1","normal_curve2",
                               "a","b","c","d","e","f",
                               "strain","replicate","condition","inocl")
write.csv(p_double_curves,paste0("output/",name_code,"p_double_curves.csv"))

######****** ODD behaviour ******#################

no_odd <- c()
ddm$odd_type <- "0"
ddm$odd_type_db <- "0"

#### store odd type in main dataframe
for(jj in 1:length(u)){ # for each strain
  
  pp <- param %>% filter(strain_name == u[jj])
  ii <- which(ddm$strain == u[jj])
  
  # Look at odd ones
  w<-which(pp$any_odd > 0)
  if(length(w > 0)){
    for(ww in 1:length(w)){
      w1 <- which(as.numeric(unlist(ddm[ii,"rep"])) == pp[w[ww],c("rep")])
      w2 <- intersect(w1,which(unlist(ddm[ii,"drytime"]) == pp[w[ww],"drytime"]))
      w3 <- intersect(w2,which(unlist(ddm[ii,'inoc']) == pp[w[ww],"inocl"]))
      
      oddv <- ii[w3]
      
      ddm[oddv,"odd_type"] <- pp[w[ww],"odd_type"] # label as odd across timeseries
      
    }
  }
  
  ddm$odd_type <- factor(ddm$odd_type, levels = c("0","01","02","03","012","013","023","0123"))
  ddm$odd_type_db <- factor(ddm$odd_type_db, levels = c("0","01","02","03","04","012","013","014","0134","034","023","0123","0124","01234"))
}


#### CUT POINT into main dataframe
ddm$shoulder_point_t <- 0
ddm$shoulder_point_v <- 0
ddm$shoulder_cut <- 0
for(jj in 1:length(u)){ # for each strain
  print(u[jj])
  pp <- param %>% filter(strain_name == u[jj])
  ii <- which(ddm$strain == u[jj])
  
  # Move over cut point or time to peak 
  for(ww in 1:length(pp[,1])){ # for every row
    w1 <- which(as.numeric(unlist(ddm[ii,"rep"])) == pp[ww,c("rep")])
    w2 <- intersect(w1,which(unlist(ddm[ii,"drytime"]) == pp[ww,"drytime"]))
    w3 <- intersect(w2,which(unlist(ddm[ii,'inoc']) == pp[ww,"inocl"]))
    
    oddv <- ii[w3] ## this 
    
    ### If have a shoulder then add in new point up to this cut 
    if(pp[ww,"shoulder_point_t"]>0){
      ddm[oddv,c("shoulder_point_t","shoulder_point_v","shoulder_cut")] <- c(pp[ww,c("shoulder_point_t","shoulder_point_v")],1) # label as odd across timeseries
    }else{
      ddm[oddv,c("shoulder_point_t","shoulder_point_v","shoulder_cut")] <- c(pp[ww,c("t_m_h_flow","v_m_h_flow")],0) # add in time to peak as cut point
    }
    
    # If have any odd behaviour then do the extra cut fit to check no earlier peaks (e.g. 11214 11.2)
    if(pp[ww,"any_odd"] > 0){
      ddm[oddv,c("shoulder_cut")] <- 1 
    }
    
  }
}

ddm1 <- ddm %>% filter(drytime %in% c(0,168))
param$rep <- as.numeric(as.character(param$rep))
param$drytime <- as.numeric(as.character(param$drytime))
#### Plot individual strain behaviour from model - odd ones highlighted
for(jj in 1:length(u)){ # for each strain
  
  dd <- ddm1 %>% filter(strain == u[jj])
  pp <- param %>% filter(strain_name == u[jj]) %>% filter(shoulder_point_v > 0)
  
  ggplot(dd, aes(x=Time, y = value_J)) + 
    geom_line(aes(group = inoc, col = factor(odd_type), linetype = factor(inoc))) + 
    facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
    scale_color_manual("Odd_type", breaks = c("0","01","02","03","012","013","023","0123"), 
                       labels = c("None","Peak","Width","Shoulder","Peak&Width","Peak&Shoulder",
                                  "Width&Shoulder","Peak Width&Shoulder"),
                       values = seq(1,8,1), drop = FALSE) + 
    scale_linetype_discrete("Inoc.") + 
    ggtitle(paste0(u[jj]," plotted:",Sys.Date()))
  ggsave(paste0("plots/",name_code,"odd_highlighted_",u[jj],".pdf")) # if any to highlight it is shown here
  
  
  #### If want to look at double curves... 
  # ggplot(dd, aes(x=Time, y = value_J)) + 
  #   geom_line(aes(group = inoc, col = odd_type_db, linetype = factor(inoc))) + 
  #   facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
  #   scale_color_manual("Odd_type", breaks = c("0","01","02","03","04","012","013","014","0134","034","023","0123","0124","01234"), 
  #                      labels = c("None","Peak","Width","Shoulder","Double","Peak&Width","Peak&Shoulder","Peak&Double",
  #                                 "Peak Shoulder&Double","Should&Double","Width&Shoulder","Peak Width&Shoulder","Peak Width&Double","All"),
  #                      values = mycolors, drop = FALSE) + 
  #   scale_linetype_discrete("Inoc.") #+ 
  # #geom_text(data = pp, aes(label = squared_dist, x = 10+as.numeric(inocl), y =as.numeric(inocl)*0.001, col = factor(inocl)),  size = 2)
  # ggsave(paste0("output/",name_code,"db_odd_highlighted_",u[jj],".pdf")) # if any to highlight it is shown here
  
  total_odd <- sum(pp$any_odd)
  
  if(total_odd < 3){print(paste0(u[jj]," has v few odd behaviours")); no_odd <- c(no_odd, u[jj])
  
  }
  
}

write.csv(ddm,paste0("output/",name_code,"_all_ddm.csv"))


#### Create cut data
## Up to time to peak or shoulder
## Remove first three hours 

ddm_cut <- ddm %>% 
  filter(Time > 3) %>% 
  group_by(strain, rep, drytime, inoc) %>% 
  mutate(cutpart = ifelse(Time < shoulder_point_t,1,0)) %>%
  filter(cutpart == 1)

write.csv(ddm_cut, "output/ddm_cut.csv")


