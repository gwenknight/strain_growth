### Analysis of sets 1, 2 & 6

###******* LOAD UP LIBRARIES AND DATA NEEDED *************#############################################################
## libraries needed

library(reshape2) # for data manipulation
library(ggplot2) # for plotting
library(gghighlight) # for highlighting lines
#library(grofit) # for fitting growth curves - no longer supported. 
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
ddm <- read_csv("output/ddm_cut.csv")[,-1]
length(unique(ddm$strain)) #98

name_code <- "cut"


###******** UNITS / DATA *************#######################################################################################################################################
dim(ddm)
ddm$value_J = ddm$value


###******* MODEL FITTING *************#######################################################################################################################################
## Fit growth curves to each set of data to determine the underlying parameters 
## Also those with odd behaviours are flagged here

## What are the strains?
u <- as.character(unique(ddm$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
#r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

## Run thru each experiment (columns in original data) for each strain
# Fit separately to each replicate as otherwise don't have enough data for later predictions
# still having up to 3 experimental drytimes despite set2 not having 24hr data
# Where the parameters for each strain are stored
drying_times <- c(0,24,168)

keep <- c()

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
### This asks if there is a shoulder? if yes then tidies up the cut point suggested by the prior analysis

for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1] %>% unlist() %>% as.character()
  for(ii in 1:length(r)){ # for each replicate: fit to all the data, not just each replicate
    for(kk in c(1,3)){ #each of the three experimental conditions (0, 24, 168): most just 0 168 now
      for(ll in 1:length(q)){ #each of the inoculums
        print(c(jj,ii,kk,ll))
        
        ### Info for this strain and condition set
        strain <- u[jj];
        replicate <- r[ii]
        condition <- drying_times[kk]
        inocl <- q[ll]
        data <- ddm # not ddm_cut?
        timepeak <- 0
        valpeak <- 0
        print(c(strain, replicate, condition, inocl))
        
        ## which rows of the data are this strain and replicate?
        wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
        wj <- intersect(wi, which(data$drytime == condition))
        w <- intersect(wj, which(data$inoc == as.numeric(inocl)))
        
        if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
          data1 <- data[w,] # %>% filter(Time > 5)# just get the data for this experiment (strain, time, value, drying time)
          
          if(data1$shoulder_cut[1] == 1){
            ## First time to peak 
            peaks_index = find_peaks(data1$value_J, m = 3)
            
            ## If there is a peak then cut at this peak time
            if(length(peaks_index) > 0){
              w <- which(data1[peaks_index,"Time"] < 3)
              if(length(w) > 0){peaks_index <- peaks_index[-w]}
              if(length(peaks_index) > 0){ # if any later than 3 yrs
                peaks_index = min(peaks_index)
                timepeak = data1[peaks_index,"Time"]
                valpeak = data1[peaks_index,"value_J"]
              }
              
            }else{ # if not peak
              valpeak = 0; timepeak = 0
              if(dim(data1)[1]>4){
                st <- c()
                for(i in 2:(-2 + dim(data1)[1])){ # fit a linear model to the data in segments
                  data11 <- data1[max(1,i-4):i,]
                  data12 <- data1[(i+1):dim(data1)[1],]
                  lm.1 <- lm(data11$value_J ~ data11$Time)
                  lm.2 <- lm(data12$value_J ~ data12$Time)
                  
                  st <- rbind(st, c(i,c(max(data11$Time),lm.1$coefficients[2],lm.2$coefficients[2])))
                }
                
                st <- as.data.frame(st)
                colnames(st) <- c("i","maxtime","f_ang","s_ang") # first angle, second angle
                st$d <- 10000
                st$d[1:(dim(st)[1]-1)] = abs(diff(st$s_ang)) # look at change in second angle: want to know when substantial change
                st_upper <- st%>% filter(maxtime > 0.50*max(st$maxtime)) # but want to be past halfway
                w1 <- which.min(st_upper$d)
                w2 <- which.min(st_upper$d[-w1])
                timepeak = min(st_upper[w1,"maxtime"], st_upper[-w1,"maxtime"][w2])
                valpeak = data1[which(data1$Time == timepeak),"value_J"]
              }else{valpeak = 0; timepeak = 0}
            }
            
            ## Check not too low: if cut point already good enough
            if(valpeak < 0.65*max(data1$value_J)){valpeak <- data1[dim(data1)[1],"value_J"]; timepeak =data1[dim(data1)[1],"Time"] }
            
            g1 <- ggplot(data1,aes(x=Time, y= value_J)) + geom_line() + 
              geom_point(data = data1[which(round(data1$Time,2) == as.numeric(round(timepeak,2))),c("Time","value_J")],col="red") + 
              ggtitle(paste0(u[jj],"_",r[ii], "_",drying_times[kk],"_",q[ll]))
            ggsave(paste0("plots/shoulder_curves/cutpoint_highlighted_",u[jj],"_",r[ii], "_",drying_times[kk],"_",q[ll],".pdf")) 
            
            data1topeak = data1[which(data1$Time < as.numeric(timepeak)),c("Time","value_J")]
            
            ## Growth Curve # (see fig 3 of vv33i07.pdf)
            ## This gives lag time and exponential growth rate cumulative
            gc_fit <- gcFitSpline(data1topeak$Time, data1topeak$value_J)
            # parameters from this fit
            s <- summary(gc_fit)
          }else{
            ## Growth Curve # (see fig 3 of vv33i07.pdf)
            timepeak = data1$shoulder_point_t[1]
            valpeak = data1$shoulder_point_v[1]
            ## This gives lag time and exponential growth rate cumulative
            gc_fit <- gcFitSpline(data1$Time, data1$value_J)
            # parameters from this fit
            s <- summary(gc_fit)
          }
          ###
          keep <- rbind(keep, as.numeric(c(strain,replicate,condition,inocl, s$mu.spline, timepeak, valpeak)))
        }
      }
    }
  }
}

keep <- as.data.frame(keep)
colnames(keep) <- c("strain","rep","drytime","inocl","cut_exp","cut_timepeak","cut_valpeak")
# keep$strain <- as.character(strain)
# keep$rep <- as.numeric(keep$rep)
# keep$drytime <- as.numeric(keep$drytime)
# keep$inocl <- as.numeric(keep$inocl)
# keep$cut_exp <- as.numeric(keep$cut_exp)
# keep$cut_timepeak <- as.numeric(keep$cut_timepeak)
# keep$cut_valpeak <- as.numeric(keep$cut_valpeak)

write.csv(keep, "output/cut_curves_fit.csv")

##### Get parameter input
param <- read.csv(paste0("output/","all_(1_13)_","all_model_fit_params.csv"))[,-1]
#ddm <- read.csv("output/all_(1_13)__all_ddm.csv")[,-1]


######****** ODD behaviour ******#################
keep$strain_name <- as.numeric(keep$strain)
keep$inoc <- as.numeric(keep$inocl)

# add on  cut 
# add in data on new exponential growth to cut
pk <- left_join(param, keep, by = c("strain_name","rep","drytime","inocl"))

dk <- left_join(ddm, keep[,c("strain","rep","drytime","inoc","cut_exp","cut_timepeak","cut_valpeak")], by = c("strain","rep","drytime","inoc"))

# Chop at this new time point for those with shoulders
dk_cut <- dk %>% 
  group_by(strain, rep, drytime, inoc) %>% 
  mutate(cutpart = ifelse(Time < cut_timepeak,1,0)) %>%
  filter(cutpart == 1)


dk_cut$odd_type <- as.character(dk_cut$odd_type)
pk$shoulder_point_t <- as.numeric(pk$shoulder_point_t)
pk$shoulder_point_v <- as.numeric(pk$shoulder_point_v)

dk_cut <- as.data.frame(dk_cut)
write.csv(dk_cut,"output/cut_all_ddm.csv")
write.csv(pk,"output/cut_all_param.csv")


