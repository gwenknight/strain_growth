### Proof of principle

### Are OD and calorimeter output comparable? 
#### Are heat flow and optical density data comparable? 

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
library(zoo)
theme_set(theme_bw(base_size=12)) # theme setting for plots: black and white (bw) and font size (24)
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(14)

setwd(here::here())

## Load in functions for this analysis
source("code/functions_for_heat_curves.R")
## Load in grofit functions if package no longer working
source("code/grofit_functions.R") # have copied those we use into this R script


## Data

data_od_orig <- read_csv("proof_of_principle/ddm_OD.csv")[,-1]
data_cs <- read_csv("proof_of_principle/ddm_CS.csv")[,-1]


### Look at data
g1 <- ggplot(data_cs, aes(x=Time, y = value_J, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Heat flow curve") +
  scale_x_continuous(lim = c(0,25)) + 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 
ggsave("proof_of_principle/CD_data.pdf")

ggplot(data_od_orig, aes(x=Time, y = value, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("OD")+ 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 
ggsave("proof_of_principle/OD_data.pdf")

### Smooth OD values
data_od <- data_od_orig %>% group_by(rep, exp, variable, strain, inoc) %>%
  mutate(ma_value = rollapply(value, 10, mean,fill = NA),
         differ = c(0,diff(ma_value))) %>%
  ungroup() %>%
  filter(Time > 1, Time < 22)

g2 <- ggplot(data_od, aes(x=Time, y = ma_value, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Optical density (600 nm)")+ 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 
ggsave("proof_of_principle/OD_data_smoothed.pdf")

g3 <- ggplot(data_od, aes(x=Time, y = differ, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Difference in optical density (600 nm) per time step")+ 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 
ggsave("proof_of_principle/OD_data_difference.pdf")

### Cumulate heat flow values
g4 <- ggplot(data_cs, aes(x=Time, y = csum, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Heat flow curve") + 
  scale_x_continuous(lim = c(0,25)) + 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 

(g1 + g4) / (g3 + g2) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
ggsave("proof_of_principle/OD_vs_CS.pdf",width = 10, height = 10)




## Join together
data_cs$value_comp <- data_cs$value_J
data_od$value_comp <- data_od$differ
data_all <- rbind(data_od %>% dplyr::select(c(Time, rep, exp, variable, strain, inoc, value_comp)), 
                  data_cs %>% dplyr::select(c(Time, rep, exp, variable, strain, inoc, value_comp))) %>% 
  filter(!is.na(value_comp))



###******* MODEL FITTING *************#######################################################################################################################################
## Fit growth curves to each set of data to determine the underlying parameters 

## What are the strains?
u <- as.character(unique(data_all$strain,stringsasFactors = FALSE)) # strains
## How many replicates? 
r <- unique(data_all$rep) # replicates
# What are the inoculums? 
q <- unique(data_all$inoc)
# What are the data types? 
da <- c("OD","Heat flow")

## Run thru each experiment (columns in original data) for each strain
# Fit separately to each replicate as otherwise don't have enough data for later predictions
# still having up to 3 experimental drytimes despite set2 not having 24hr data
# Where the parameters for each strain are stored
param <- matrix(0, length(u)*length(r)*length(q)*3*2, 19); # number of strains x number of replicates x number of experimental conditions x number of experiments
index <- 1 # for counting 
max <- c() # for calibration

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(data_all %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1])
  for(ii in 1:length(r)){ # for each replicate: fit to all the data, not just each replicate
    for(ll in 1:length(q)){ #each of the inocula
      for(exp in 1:2){ #each data type
        
        strain <- u[jj];
        replicate <- as.character(r[ii])
        inocl <- q[ll]
        exp <- da[exp]
        
        wi <- intersect(which(data_all$strain == strain),which(data_all$rep == replicate)) # if fit to each replicate
        wj <- intersect(wi, which(data_all$exp == exp)) # which type of data
        w <- intersect(wj, which(data_all$inoc == as.numeric(inocl))) # which inoculum
        
        if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
          print(c(strain, replicate, inocl, exp))
          data1 <- data_all[w,] # Grab data
          
          p <- cut_extract(data1, "Time", "value_comp", paste(strain, replicate, exp, inocl, sep="_"), early_cut = 0) ### NEW function: runs on any timeseries: gives baseline and cut parameter analysis
          
          ## Required parameters
          
          param[index,] <- unlist(c(strain, replicate, exp, inocl, p$param))
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
colnames(param) <- c("strain_name","rep","experiment","inocl",
                     "t_m_h_flow", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v",
                     "cut_exp", "timepeak", "valpeak")

w<-which(param$lag!=0); param <- param[w,] # remove 0 values
dim(param)

param$valpeak <- as.numeric(param$valpeak)
param$timepeak <- as.numeric(param$timepeak)





## STORE CUT POINT IN MAIN TIMESERIS 
data_all$shoulder_point_t <- 0
data_all$shoulder_point_v <- 0
data_all$shoulder_cut <- 0

#### store odd type in main timeseries dataframe
for(jj in 1:length(u)){ # for each strain
  
  pp <- param %>% filter(strain_name == u[jj])
  ii <- which(data_all$strain == u[jj])
  print(u[jj])
  
  # Move over cut point or time to peak 
  for(ww in 1:length(pp[,1])){ # for every row
    w1 <- which(unlist(data_all[ii,"rep"]) == pp[ww,c("rep")])
    w2 <- intersect(w1,which(unlist(data_all[ii,"exp"]) == pp[ww,"experiment"]))
    w3 <- intersect(w2,which(unlist(data_all[ii,'inoc']) == as.numeric(pp[ww,"inocl"])))
    
    oddv <- ii[w3] ## this 
    
    ### If have a shoulder then add in new point up to this cut 
    if(pp[ww,"shoulder_point_t"]>0){
      data_all[oddv,c("shoulder_point_t")] <- c(as.numeric(unlist(pp[ww,c("shoulder_point_t","shoulder_point_v")])),1)[1]
      data_all[oddv,c("shoulder_point_v")] <- c(as.numeric(unlist(pp[ww,c("shoulder_point_t","shoulder_point_v")])),1)[2]
      data_all[oddv,c("shoulder_cut")] <- c(as.numeric(unlist(pp[ww,c("shoulder_point_t","shoulder_point_v")])),1)[3]# label as odd across timeseries
    }else{
      data_all[oddv,c("shoulder_point_t")] <- c(pp[ww,c("timepeak", "valpeak")],0)[1]
      data_all[oddv,c("shoulder_point_v")] <- c(pp[ww,c("timepeak", "valpeak")],0)[2]
      data_all[oddv,c("shoulder_cut")] <- c(pp[ww,c("timepeak", "valpeak")],0)[3] # add in time to peak as cut point
    }
    
    
  }
  
  
}

### Output
write_csv(data_all, "proof_of_principle/data_all.csv")
write_csv(param, "proof_of_principle/parameters_output.csv")


### PLOTS
ggplot(data_all, aes(x=Time, y = value_comp, group = interaction(rep, inoc, exp))) + 
  geom_line(col = "grey", alpha = 0.2) + 
  facet_grid(exp~strain, scales = "free") + 
  geom_line(data = data_all %>% group_by(strain, rep, exp, variable) %>% filter(Time <= shoulder_point_t),aes(col = factor(inoc))) + 
  geom_point(data = data_all, aes(x=shoulder_point_t, y = shoulder_point_v)) + 
  scale_x_continuous(lim = c(0,25)) + 
  scale_y_continuous("Value for comparison (OD (600nm) or heat flow value)") + 
  scale_color_discrete("Inoculum", breaks = c(2,3,4,5,6), labels = c(str2expression("10^2"),str2expression("10^3"),
                                                                     str2expression("10^4"),str2expression("10^5"),
                                                                     str2expression("10^6"))) 
ggsave("proof_of_principle/od_vs_cs_comparison.pdf", width = 10, height = 10)  

param %>% ungroup() %>% dplyr::select(strain_name, inocl, rep, experiment, timepeak) %>% pivot_wider(names_from = experiment, values_from = timepeak) %>%
  ggplot(aes(x=OD, y = CS)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1) + 
  scale_y_continuous(lim = c(0,25))
ggsave("proof_of_principle/od_vs_cs_linear.pdf")

## Correlation in time to first peak
corr <- c()
for(i in unique(param$strain_name)){
  corr <- rbind(corr, 
                c(i, round(cor(param %>% filter(strain_name == i, experiment == "OD") %>% dplyr::select(timepeak),
            param %>% filter(strain_name == i, experiment == "CS") %>% dplyr::select(timepeak)),2)))
}

corr
write.csv(corr, "proof_of_principle/corr_output.csv")

