##### Data formatting for set13

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
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

## Load in grofit functions if package no longer working
source("code/grofit_functions.R")

## Load in functions for this analysis
source("code/functions_for_heat_curves.R")

### TO do this for multiple files: 
## put the experiment/file names in here: 
names_experiments<- c(29, 44, 45, 46, 30, 31, 33, 47)

dd <- c()

## loop thru experiment names
for (i in 1:length(names_experiments)){ 
  da <- paste0("data/exp",names_experiments[i],"a.txt")
  db <- paste0("data/exp",names_experiments[i],"b.txt")
  dc <- paste0("data/exp",names_experiments[i],"c.txt")
  
  if(file.exists(da)){dda <- read.table(da,header = TRUE)} # read in the data 
  if(file.exists(db)){ddb <- read.table(db,header = TRUE)} # read in the data 
  if(file.exists(dc)){ddc <- read.table(dc,header = TRUE)} # read in the data 
  
  # which experiment? hours of drying
  dda$exp <- 0
  ddb$exp <- 24
  ddc$exp <- 168
  
  # which replicate?
  dda$rep <- names_experiments[i]
  ddb$rep <- names_experiments[i]
  ddc$rep <- names_experiments[i]
  
  dd <- rbind(rbind(rbind(dd, dda), ddb), ddc)
}


## convert seconds to hours
dd$Time <- dd$Time / 3600

## pull down experiments into a column (instead of a separate column for each)
ddm <- dd %>% pivot_longer(cols = c(B2:D6,E2:E6),names_to = "variable")

## Cumulative heat curve in joules
# Heat flow in muW. 1 W = 1 Joule / Second. 
# # Assume that heat curve output is x muW over that interval of 0.25hrs
# intervals = 0.25/60/60 # time interval between readings in seconds
# ddm$value_J <-  ddm$value * intervals # convert to Joules (W = J/s => J = W*s)
# # Want to normalise - i.e. remove baseline value so that starts at 0 but tricky as some negative? not done here
# ddm <- ddply(ddm,.(name, rep,exp),transform,csum=cumsum(value_J))
# 

#Data in mJ, convert to J
ddm$value <- (ddm$value/1000)

## Plot the experiments
## Raw data
ggplot(subset(ddm,exp==0), aes(x=Time,y=value,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
  facet_wrap(~variable, nrow = 2) + ggtitle("Reference") + scale_y_continuous("Heat (J)")

ggplot(subset(ddm,exp==168), aes(x=Time,y=value,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
  facet_wrap(~variable, nrow = 2) + ggtitle("168hr drying") + scale_y_continuous("Heat (J)")



# Remove the contaminated data: MSSAs
# REFERENCE
# 24hr drying
w1 <- intersect(intersect(which(ddm$variable == "D5"), which(ddm$rep == 46)),which(ddm$exp == 24))
w2 <- intersect(intersect(which(ddm$variable == "E5"), which(ddm$rep == 46)),which(ddm$exp == 24))
ddm <-ddm[-c(w1,w2),]

# 168hr drying
w1 <- intersect(intersect(which(ddm$variable == "E3"), which(ddm$rep == 29)),which(ddm$exp == 168))
ddm <-ddm[-c(w1),]

# Remove the contaminated data: MRSAs
# REFERENCE
w1 <- intersect(which(ddm$variable == "E5"), which(ddm$rep == 33))
ddm <-ddm[-c(w1),]

# 24hr drying
w1 <- intersect(intersect(which(ddm$variable == "B6"), which(ddm$rep == 33)),which(ddm$exp == 24))
ddm <-ddm[-c(w1),]

# 168hr drying
w1 <- intersect(intersect(which(ddm$variable == "E2"), which(ddm$rep == 33)),which(ddm$exp == 168))
ddm <-ddm[-c(w1),]

### Replace 
ddm$strain <- ""; ddm$inoc <- 0

ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

# strain names
w<- intersect(which(ddm$rep %in% c(29, 44, 45, 46)),which(ddm$strain_label == "B"))
ddm[w,"strain"] = "Newman"
w<- intersect(which(ddm$rep %in% c(29, 44, 45, 46)),which(ddm$strain_label == "C"))
ddm[w,"strain"] = "RWW12"
w<- intersect(which(ddm$rep %in% c(29, 44, 45, 46)),which(ddm$strain_label == "D"))
ddm[w,"strain"] = "SA3297"
w<- intersect(which(ddm$rep %in% c(29, 44, 45, 46)),which(ddm$strain_label == "E"))
ddm[w,"strain"] = "SA2704"

# strain names
w<- intersect(which(ddm$rep %in% c(30, 31, 33, 47)),which(ddm$strain_label == "B"))
ddm[w,"strain"] = "RWW146"
w<- intersect(which(ddm$rep %in% c(30, 31, 33, 47)),which(ddm$strain_label == "C"))
ddm[w,"strain"] = "SAC042W"
w<- intersect(which(ddm$rep %in% c(30, 31, 33, 47)),which(ddm$strain_label == "D"))
ddm[w,"strain"] = "Mu50"
w<- intersect(which(ddm$rep %in% c(30, 31, 33, 47)),which(ddm$strain_label == "E"))
ddm[w,"strain"] = "M116"

# inoculum size
w<- which(ddm$inoc_name == "2")
ddm[w,"inoc"] = 6
w<- which(ddm$inoc_name == "3")
ddm[w,"inoc"] = 5
w<- which(ddm$inoc_name == "4")
ddm[w,"inoc"] = 4
w<- which(ddm$inoc_name == "5")
ddm[w,"inoc"] = 3
w<- which(ddm$inoc_name == "6")
ddm[w,"inoc"] = 2



## Plot tidy data 
ggplot(subset(ddm,exp==0), aes(x=Time,y=value,group = interaction(inoc,factor(rep)),colour=factor(inoc))) + geom_line(lwd = 1.5) + 
  facet_wrap(~strain, nrow = 2) + ggtitle("Reference") + scale_y_continuous("Heat (J)")


# OUTPUT
ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label
ddm$drytime <- ddm$exp
ddm <- ddm[,c("Time","rep","exp","variable","value","strain","inoc","drytime")]
write.csv(ddm, "data/ddm_set14.csv")
