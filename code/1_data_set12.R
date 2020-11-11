##### Data formatting for set12

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
names_experiments<- c(12.1,12.2,12.3)

## loop thru experiment names to get data
dd <- c()

for (i in 1:length(names_experiments)){ 
  
  name <- names_experiments[i]
  
  da <- paste0("data/set",name,"t0.txt") ##change name data storage folder
  db <- paste0("data/set",name,"t7.txt") ##change name data storage folder
  
  if(file.exists(da)){dda <- read.table(da,header = TRUE)} # read in the data 
  if(file.exists(db)){ddb <- read.table(db,header = TRUE)} # read in the data 
  
  # which experiment?
  dda$exp <- "t0"
  ddb$exp <- "t7"
  
  # which replicate?
  dda$rep <- name
  ddb$rep <- name
  
  dd <- rbind(rbind(dd, dda), ddb)
  
}


## convert seconds to hours
dd$Time <- dd$Time / 3600

## pull down experiments into a column (instead of a separate column for each)
ddm <-reshape2::melt(dd, id.vars=c("Time","rep","exp"))

## Cumulative heat curve in joules
# Heat flow in muW. 1 W = 1 Joule / Second. 
# Assume that heat curve output is x muW over that interval of 0.25hrs
#intervals = 0.25/60/60 # time interval between readings in seconds
#ddm$value_J <-  ddm$value * intervals # convert to Joules (W = J/s => J = W*s)
# Want to normalise - i.e. remove baseline value so that starts at 0 but tricky as some negative? not done here
#ddm <- ddply(ddm,.(variable, rep,exp),transform,csum=cumsum(value_J))

#Data in mJ, convert to J
ddm$value <- (ddm$value/1000)

## Plot the experiments
## Raw data
ggplot(subset(ddm,exp=="t0"), aes(x=Time,y=value,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
  facet_wrap(~variable, nrow = 4) + ggtitle("Reference") + scale_y_continuous("Heat (J)")
#ggsave("output/Raw_data_Reference_set12.tiff", width = 15, height = 10)
ggplot(subset(ddm,exp=="t7"), aes(x=Time,y=value,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
  facet_wrap(~variable, nrow = 4) + ggtitle("168hr drying") + scale_y_continuous("Heat (J)")
#ggsave("output/Raw_data_168hr drying_set12.tiff", width = 15, height = 10)


# Remove the contaminated data
# REFERENCE
#w1 <- intersect(which(ddm$variable == "B4"), which(ddm$rep == 1.1))
#w2 <- intersect(which(ddm$variable == "B6"), which(ddm$rep == 1.1))
#w3 <- intersect(which(ddm$variable == "C6"), which(ddm$rep == 26))
#w4 <- intersect(which(grepl("E",ddm$variable)), which(ddm$rep == 23))

#ddm <-ddm[-c(w1,w2,w3,w4,w5,w7,w8),]

# 1wk drying
#w1 <- intersect(intersect(which(ddm$variable == "D5"), which(ddm$rep == 1.1)),which(ddm$exp == "b"))
#w2 <- intersect(intersect(which(ddm$variable == "E5"), which(ddm$rep == 1.1)),which(ddm$exp == "b"))
#w3 <- intersect(intersect(which(ddm$variable == "D5"), which(ddm$rep == 26)),which(ddm$exp == "b"))
#w4 <- intersect(intersect(which(ddm$variable == "B5"), which(ddm$rep == 26)),which(ddm$exp == "b"))

#ddm <-ddm[-c(w1,w2),]

## Plot tidy data and save 
ggplot(subset(ddm,exp=="t0"), aes(x=Time,y=value,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
  facet_wrap(~variable, nrow = 4) + ggtitle("Reference") + scale_y_continuous("Heat (J)")
#ggsave("output/Raw_data_Reference_set12_tidy.tiff", width = 15, height = 10)
ggplot(subset(ddm,exp=="t7"), aes(x=Time,y=value,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
  facet_wrap(~variable, nrow = 4) + ggtitle("168hr drying") + scale_y_continuous("Heat (J)")
#ggsave("output/Raw_data_168hr drying_set12_tidy.tiff", width = 15, height = 10)


## Cumulative heat output separate
#ggplot(subset(ddm,exp=="a"), aes(x=Time,y=csum,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
#  facet_wrap(~variable) + ggtitle("Cumulative, Reference")+ scale_y_continuous("Heat in J")
#ggsave("output/Heat curve_Reference.pdf")
#ggsave("output/Heat curve_Reference.tiff")
#ggplot(subset(ddm,exp=="b"), aes(x=Time,y=csum,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
#  facet_wrap(~variable) + ggtitle("Cumulative, 24hr drying")+ scale_y_continuous("Heat in J")
#ggsave("output/Heat curve_24hr drying.pdf")
#ggsave("output/Heat curve_24hr drying.tiff")
#ggplot(subset(ddm,exp=="c"), aes(x=Time,y=csum,group = factor(rep),colour=factor(rep))) + geom_line(lwd = 1.5) + 
#  facet_wrap(~variable) + ggtitle("Cumulative, 168hr drying")+ scale_y_continuous("Heat in J")
#ggsave("output/Heat curve_168hr drying.pdf")
#ggsave("output/Heat curve_168hr drying.tiff")


## Cumulative heat output together
#ggplot(ddm, aes(x=Time,y=csum, group = variable, col = factor(rep))) + facet_wrap(~exp) + geom_line(aes(col=variable)) + ggtitle("Raw data")
#ggsave("output/all_toge_csum_tog.pdf")

## To plot a single piece of data use:
#gd <- dd[,c("Time","C5")] ## where "C5" is replaced by the experiment you want
#ggplot(gd, aes(x=Time,y=C5)) + geom_line() + ggtitle("Raw data")
#ggsave(paste0("output/",name,"_C5.pdf"))

### Replace - CHECK VALERIE! 
ddm$strain <- 0; ddm$inoc <- 0

ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

# strain names
w<- which(ddm$variable == "B2"|ddm$variable == "B3"|ddm$variable == "B4")
ddm[w,"strain"] = "11237"
w<- which(ddm$variable == "C2"|ddm$variable == "C3"|ddm$variable == "C4")
ddm[w,"strain"] = "11239"
w<- which(ddm$variable == "D2"|ddm$variable == "D3"|ddm$variable == "D4")
ddm[w,"strain"] = "11245"
w<- which(ddm$variable == "E2"|ddm$variable == "E3"|ddm$variable == "E4")
ddm[w,"strain"] = "11250"
w<- which(ddm$variable == "B5"|ddm$variable == "B6"|ddm$variable == "B7")
ddm[w,"strain"] = "11257"
w<- which(ddm$variable == "C5"|ddm$variable == "C6"|ddm$variable == "C7")
ddm[w,"strain"] = "11259"
w<- which(ddm$variable == "D5"|ddm$variable == "D6"|ddm$variable == "D7")
ddm[w,"strain"] = "11263"
w<- which(ddm$variable == "E5"|ddm$variable == "E6"|ddm$variable == "E7")
ddm[w,"strain"] = "11265"

# inoculum size
w<- which(ddm$variable == "B2"|ddm$variable == "C2"|ddm$variable == "D2"|ddm$variable == "E2"|
            ddm$variable == "B5"|ddm$variable == "C5"|ddm$variable == "D5"|ddm$variable == "E5")
ddm[w,"inoc"] = 5
w<- which(ddm$variable == "B3"|ddm$variable == "C3"|ddm$variable == "D3"|ddm$variable == "E3"|
            ddm$variable == "B6"|ddm$variable == "C6"|ddm$variable == "D6"|ddm$variable == "E6")
ddm[w,"inoc"] = 4
w<- which(ddm$variable == "B4"|ddm$variable == "C4"|ddm$variable == "D4"|ddm$variable == "E4"|
            ddm$variable == "B7"|ddm$variable == "C7"|ddm$variable == "D7"|ddm$variable == "E7")
ddm[w,"inoc"] = 3

### Experiment change to hours drying
wa<-which(ddm$exp == "t0")
wb<-which(ddm$exp == "t7")
ddm[wa,"drytime"] <- 0
ddm[wb,"drytime"] <- 168

### CHECK ALL OK? plots
ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(group = inoc)) + facet_wrap(drytime~strain)

#### OUTPUT
ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label

write.csv(ddm, "data/ddm_set12.csv")
