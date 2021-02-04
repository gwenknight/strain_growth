##Growth curve fitting using grofit curves##
##based on Gwen Knight's survival model script##

##Libraries##

library(reshape2) # for data manipulation
library(ggplot2) # for plotting
library(grofit) # for fitting growth curves
library(plyr) # for data manipulation
library(data.table)
library(reshape2)
library(magrittr)
library(dplyr)
library(MASS)
library(grid)
library(gridExtra)
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

##Home
home <- "D:/ErasmusMC/2020 OD concentratiereeks/Analyse/Growth curve analysis"
setwd(home)
getwd()

## Load in grofit functions if package no longer working
source("code/grofit_functions.R")

## Load in functions for this analysis
source("20_functions_for_heat_curves.R")

##Organize CS heat flow data - in muW ####
dd <- read.table("proof_of_principle/CS.txt",header = TRUE)
dd$Time <- dd$Time / 3600
ddm <-reshape2::melt(dd, id.vars=c("Time","rep","exp"))

# Cumulative heat curve in joules
# Heat flow in muW. 1 W = 1 Joule / Second. 
# Assume that heat curve output is x muW over that interval of 1/6hrs
intervals = (1/6)/60/60 # time interval between readings in seconds (every 10 min)
ddm$value_J <-  ddm$value * intervals # convert to Joules (W = J/s => J = W*s)
ddm <- ddply(ddm,.(variable, rep,exp),transform,csum=cumsum(value_J))


ddm$strain <- 0; ddm$inoc <- 0
ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

# strain names
w<- which(ddm$strain_label == "A")
ddm[w,"strain"] = "SAC042W"
w<- which(ddm$strain_label == "B")
ddm[w,"strain"] = "M116"
w<- which(ddm$strain_label == "C")
ddm[w,"strain"] = "SA2704"
w<- which(ddm$strain_label == "D")
ddm[w,"strain"] = "SAMUP15a"
# inoculum size
w<- which(ddm$inoc_name == "1")
ddm[w,"inoc"] = 6
w<- which(ddm$inoc_name == "2")
ddm[w,"inoc"] = 5
w<- which(ddm$inoc_name == "3")
ddm[w,"inoc"] = 4
w<- which(ddm$inoc_name == "4")
ddm[w,"inoc"] = 3
w<- which(ddm$inoc_name == "5")
ddm[w,"inoc"] = 2

# Check plots
ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(group = interaction(rep,inoc), col = factor(inoc))) + facet_wrap(~strain)
ggplot(ddm, aes(x=Time,y = value_J)) + geom_line(aes(group = interaction(rep,inoc), col = factor(inoc))) + facet_wrap(~strain)
#ggplot(ddm, aes(x=Time,y = csum)) + geom_line(aes(group = inoc)) + facet_wrap(~strain)
#ggsave("output/CS_curves.tiff")
#ggsave("output/CS_curves.png")

# OUTPUT
ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label
write.csv(ddm, "data/ddm_CS.csv")

# ##Organize CS heat data - in mJ ####
# dd <- read.table("data/CSHeat.txt",header = TRUE)
# dd$Time <- dd$Time / 3600
# ddm <-reshape2::melt(dd, id.vars=c("Time","rep","exp"))
# 
# ddm$strain <- 0; ddm$inoc <- 0
# ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))
# 
# # strain names
# w<- which(ddm$strain_label == "A")
# ddm[w,"strain"] = "SAC042W"
# w<- which(ddm$strain_label == "B")
# ddm[w,"strain"] = "M116"
# w<- which(ddm$strain_label == "C")
# ddm[w,"strain"] = "SA2704"
# w<- which(ddm$strain_label == "D")
# ddm[w,"strain"] = "SAMUP15a"
# # inoculum size
# w<- which(ddm$inoc_name == "1")
# ddm[w,"inoc"] = 6
# w<- which(ddm$inoc_name == "2")
# ddm[w,"inoc"] = 5
# w<- which(ddm$inoc_name == "3")
# ddm[w,"inoc"] = 4
# w<- which(ddm$inoc_name == "4")
# ddm[w,"inoc"] = 3
# w<- which(ddm$inoc_name == "5")
# ddm[w,"inoc"] = 2
# 
# # Check plots
# ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(group = inoc)) + facet_wrap(~strain)
# #ddms <- subset(ddm, strain == "SAC042W")
# #heatplot <- ggplot(ddms, aes(x=Time,y = value)) + geom_line(aes(group = inoc))
# #heatplot
# #ggsave("output/CSHeat_curves.png")
# 
# # OUTPUT
# ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label
# write.csv(ddm, "data/ddm_CSHeat.csv")

##Organize OD data ####
dd <- read.table("proof_of_principle//OD.txt",header = TRUE)
dd$Time <- dd$Time / 3600
ddm <-reshape2::melt(dd, id.vars=c("Time","rep","exp"))

ddm$strain <- 0; ddm$inoc <- 0
ddm[,c("strain_label","inoc_name")]<-colsplit(ddm$variable, "(?<=\\p{L})(?=[\\d+$])", c("char", "digit"))

# strain names
w<- which(ddm$strain_label == "A")
ddm[w,"strain"] = "SAC042W"
w<- which(ddm$strain_label == "B")
ddm[w,"strain"] = "M116"
w<- which(ddm$strain_label == "C")
ddm[w,"strain"] = "SA2704"
w<- which(ddm$strain_label == "D")
ddm[w,"strain"] = "SAMUP15a"
# inoculum size
w<- which(ddm$inoc_name == "1")
ddm[w,"inoc"] = 6
w<- which(ddm$inoc_name == "2")
ddm[w,"inoc"] = 5
w<- which(ddm$inoc_name == "3")
ddm[w,"inoc"] = 4
w<- which(ddm$inoc_name == "4")
ddm[w,"inoc"] = 3
w<- which(ddm$inoc_name == "5")
ddm[w,"inoc"] = 2

# Check plots
ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(group = interaction(rep,inoc), col = factor(inoc))) + facet_wrap(~strain)
#ggsave("output/OD_curves.png")
#ddms <- subset(ddm, strain == "SAC042W")
#odplot <- ggplot(ddms, aes(x=Time,y = value)) + geom_line(aes(group = inoc))
#odplot
#sacplot <- grid.arrange(odplot, heatplot, nrow = 1, top="SAC042W")
#ggsave("output/sacplot.png", sacplot, width = 12, height = 5)

# OUTPUT
ddm <- subset(ddm, select = -c(strain_label,inoc_name)) # remove label
write.csv(ddm, "data/ddm_OD.csv")

##Model fitting CS heat flow data ####

#Grab CS data
ddm <- as.data.table(read.csv("data/ddm_CS.csv")[,-1])

#Create matrices
u <- as.character(unique(ddm$strain)) #unique strains
r <- unique(ddm$rep) #unique replicates
q <- unique(ddm$inoc) #unique inocula

param_n <- matrix(0, length(u)*length(r)*length(q)*1, 3); # Strain x rep x inoc x number of experimental conditions x number of calculated variables
param <- matrix(0, length(u)*length(r)*length(q)*1, 5); #Strain x rep x inoc x number of experimental conditions x number of calculated variables
index <- 1

# Calculate parameters and fill matrices

for (jj in 1:length(u)){ #for each strain
  for (ii in 1:length(r)){# for each replicate
    for (ll in 1:length(q)){# for each inoculum 
      print(c(jj, ii, ll))
      
      #select data for unique combination of u, r and q
      w <- intersect(which(ddm$inoc == q[ll]),(intersect(which(ddm$strain == u[jj]), which(ddm$rep == r[ii]))))
      #w <- intersect(which(ddm$inoc == 6),(intersect(which(ddm$strain == "SAC042W"), which(ddm$rep == "A"))))
      
      if(length(w) > 0){
        data1 <- ddm[w,]
        
        #Fit growth curve
        gc_fit <- gcFitSpline(data1$Time, data1$csum)
        #Calculate growth parameters
        wmax <- which.max(data1$value_J)
        time_max_heat_flow <- as.numeric(data1[wmax, "Time"])
        value_max_heat_flow <- as.numeric(data1[wmax, "value"])
        print(c(u[jj], time_max_heat_flow, value_max_heat_flow))
        
        #Plot data1, csum and fit
        #data1
        datam <- reshape2::melt(data1[,c("Time", "value_J", "csum", "strain")], id.vars = c("Time", "strain"))
        gg <- ggplot(datam, aes(x=Time, y=value)) + geom_point() + facet_wrap(~strain, scales = "free")
        #model fit
        gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
        colnames(gc_df) <- c("Time", "csum")
        #gc_df$value_J <- c(gc_df[1, "csum"], diff(gc_df$csum))
        gc_dfm <- reshape2::melt(gc_df, id.vars = "Time")
        #add fit to plot
        gg <- gg + geom_line(data = gc_dfm, aes(x=Time, y=value), col = "red")
        gg
        ##save
        #ggsave(paste0("output/",u[jj],"_rep_", r[ii],"_inoc_", q[ll], "_model_fit_CS.pdf"))
        
        #Plot heat flow, reference line time peak
        ggplot(data1, aes(x=Time, y=value)) + geom_point() +
          geom_vline(xintercept = time_max_heat_flow, col = "black") + geom_text(aes(time_max_heat_flow, 0, label = round(time_max_heat_flow, digits = 2))) +
          scale_x_continuous("Time (h)", limits = c(0, 25), minor_breaks = seq(0, 25, 5)) + 
          scale_y_continuous(expression(paste("Heatflow (", mu, "W)")), limits = c(0,100), breaks = seq(0, 100, 20)) +
          ggtitle(paste("Strain",u[jj],"rep",r[ii]," inoc",q[ll]))
        #ggsave(paste0("output/",u[jj],"_rep_", r[ii],"_inoc_", q[ll], "_heat flow_timemax_CS.png"))
        
        #Save parameters
        s <- summary(gc_fit)
        #u_d <- levels(droplevels(u))
        #param_n[index,] <- c(unlist(strsplit(u_d[jj],"")), r[ii], q[ll])
        param_n[index,] <- c(u[jj], r[ii], q[ll])
        param[index,] <- c(time_max_heat_flow, value_max_heat_flow, s$mu.spline, s$lambda.spline,s$integral.spline)
        #param[index,] <- c(time_max_heat_flow, value_max_heat_flow, gc_fit$parameters$mu, gc_fit$parameters$lambda, gc_fit$parameters$integral)
        index <- index + 1
        
      }
    }
  }
}


param <- as.data.frame(param)
param_n <- as.data.frame(param_n)
colnames(param_n) <- c("strain","rep","inoc")
colnames(param) <- c("tmax", "ymax", "expgr","lag","auc")

param_orig_cs <- param
param$strain <- param_n$strain
param$inoc <- param_n$inoc
param$rep <- param_n$rep

#w<-which(param[,1] == 0)
#if(length(w) > 0){param <- param[-w,]} # remove 0 rows only if there are some

dim(param) 

#Output
write.csv(param, "output/param_CS.csv")

# ##Model fitting CS heat data ####
# 
# #Grab CSHeat data
# ddm <- as.data.table(read.csv("data/ddm_CSHeat.csv")[,-1])
# 
# #Create matrices
# u <- as.character(unique(ddm$strain)) #unique strains
# r <- unique(ddm$rep) #unique replicates
# q <- unique(ddm$inoc) #unique inocula
# 
# param_n <- matrix(0, length(u)*length(r)*length(q)*1, 3); # Strain x rep x inoc x number of experimental conditions x number of calculated variables
# param <- matrix(0, length(u)*length(r)*length(q)*1, 5); #Strain x rep x inoc x number of experimental conditions x number of calculated variables
# index <- 1
# 
# max(ddm$value)
# min(ddm$value)
# # Calculate parameters and fill matrices
# 
# for (jj in 1:length(u)){ #for each strain
#   for (ii in 1:length(r)){# for each replicate
#     for (ll in 1:length(q)){# for each inoculum 
#       print(c(jj, ii, ll))
#       
#       #select data for unique combination of u, r and q
#       w <- intersect(which(ddm$inoc == q[ll]),(intersect(which(ddm$strain == u[jj]), which(ddm$rep == r[ii]))))
#       #w <- intersect(which(ddm$inoc == 6),(intersect(which(ddm$strain == "SAC042W"), which(ddm$rep == "A"))))
#       
#       if(length(w) > 0){
#         data1 <- ddm[w,]
#         
#         #Fit growth curve
#         gc_fit <- gcFitSpline(data1$Time, data1$value)
#         #plot.gcFitSpline(gc_fit)
#         #dev.copy(png, "output/_modelspline_OD.png")
#         
#         #Calculate growth parameters
#         y.spl     <- smooth.spline(data1$Time, data1$value)
#         dydt.spl   <- predict(y.spl, data1$Time, deriv = 1)
#         inde      <- which.max(dydt.spl$y)         
#         t.max      <- dydt.spl$x[inde]
#         dydt.max   <- max(dydt.spl$y)
#         y.max      <- y.spl$y[inde]
#         mu.spl     <- dydt.max;
#         b.spl      <- y.max-dydt.max*t.max
#         
#         #Plot data1, splines, tangent and ablines
#         gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
#         colnames(gc_df) <- c("Time", "value")
#         
#         ggplot(data1, aes(x=Time, y=value)) + geom_point() +
#           geom_line(data = gc_df, aes(x=Time, y=value), col = "red") +
#           geom_abline(intercept = b.spl, slope = mu.spl, col = "red") +
#           geom_vline(xintercept = t.max, col = "black") + geom_text(aes(t.max, 0, label = round(t.max, digits = 2))) +
#           scale_x_continuous("Time (h)", limits = c(0, 25), minor_breaks = seq(0, 25, 5)) + 
#           scale_y_continuous("Heat (mJ)", limits = c(-5, 2750), breaks = seq(0, 2500, 500)) +
#           ggtitle(paste("Strain",u[jj],"rep",r[ii]," inoc",q[ll])) +
#           theme(plot.title=element_text(size=20, hjust=0.5, vjust=-1))
#         #ggsave(paste0("output/",u[jj],"_rep_", r[ii],"_inoc_", q[ll], "_model_fit_CSHeat.png"))
#         
#         #Save parameters
#         s <- summary(gc_fit)
#         t.max <- as.numeric(t.max)
#         y.max <- as.numeric(y.max)
#         param_n[index,] <- c(u[jj], r[ii], q[ll])
#         param[index,] <- c(t.max, y.max, s$mu.spline, s$lambda.spline,s$integral.spline)
#         #param[index,] <- c(s$mu.spline, s$lambda.spline,s$integral.spline)
#         index <- index + 1
#         
#       }
#     }
#   }
# }
# 
# 
# param <- as.data.frame(param)
# param_n <- as.data.frame(param_n)
# colnames(param_n) <- c("strain","rep","inoc")
# colnames(param) <- c("tmax", "ymax", "expgr","lag","auc")
# 
# param_orig_csheat <- param
# param$strain <- param_n$strain
# param$inoc <- param_n$inoc
# param$rep <- param_n$rep
# 
# dim(param)  
# 
# #Output
# write.csv(param, "output/param_CSHeat.csv")


##Model fitting OD data ####

#Grab OD data
ddm <- as.data.table(read.csv("data/ddm_OD.csv")[,-1])
#ggplot(ddm, aes(x=Time,y = value)) + geom_line(aes(group = inoc)) + facet_wrap(~strain)

#Create matrices
u <- as.character(unique(ddm$strain)) #unique strains
r <- unique(ddm$rep) #unique replicates
q <- unique(ddm$inoc) #unique inocula

param_n <- matrix(0, length(u)*length(r)*length(q)*1, 3); # Strain x rep x inoc x number of experimental conditions x number of calculated variables
param <- matrix(0, length(u)*length(r)*length(q)*1, 5); #Strain x rep x inoc x number of experimental conditions x number of calculated variables
index <- 1

# Calculate parameters and fill matrices

for (jj in 1:length(u)){ #for each strain
  for (ii in 1:length(r)){# for each replicate
    for (ll in 1:length(q)){# for each inoculum 
      print(c(jj, ii, ll))
      
      #select data for unique combination of u, r and q
      w <- intersect(which(ddm$inoc == q[ll]),(intersect(which(ddm$strain == u[jj]), which(ddm$rep == r[ii]))))
      #w <- intersect(which(ddm$inoc == 6),(intersect(which(ddm$strain == "SAC042W"), which(ddm$rep == "A"))))
      
      if(length(w) > 0){
        data1 <- ddm[w,]
        
        #Fit growth curve
        gc_fit <- gcFitSpline(data1$Time, data1$value)
        #plot.gcFitSpline(gc_fit)
        #dev.copy(png, "output/_modelspline_OD.png")
        
        #Calculate growth parameters
        y.spl     <- smooth.spline(data1$Time, data1$value)
        dydt.spl   <- predict(y.spl, data1$Time, deriv = 1)
        inde      <- which.max(dydt.spl$y)         
        t.max      <- dydt.spl$x[inde]
        dydt.max   <- max(dydt.spl$y)
        y.max      <- y.spl$y[inde]
        mu.spl     <- dydt.max;
        b.spl      <- y.max-dydt.max*t.max
        
        #Plot data1, splines, tangent and ablines
        gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
        colnames(gc_df) <- c("Time", "value")
        
        ggplot(data1, aes(x=Time, y=value)) + geom_point() +
          geom_line(data = gc_df, aes(x=Time, y=value), col = "red") +
          geom_abline(intercept = b.spl, slope = mu.spl, col = "red") +
          geom_vline(xintercept = t.max, col = "black") + geom_text(aes(t.max, 0, label = round(t.max, digits = 2))) +
          scale_x_continuous("Time (h)", limits = c(0, 25), minor_breaks = seq(0, 25, 5)) + 
          scale_y_continuous("OD at 600nm", limits = c(0, 1), breaks = seq(0.25, 1, 0.25)) +
          ggtitle(paste("Strain",u[jj],"rep",r[ii]," inoc",q[ll]))
        #ggsave(paste0("output/",u[jj],"_rep_", r[ii],"_inoc_", q[ll], "_model_fit_OD.png"))
        
        #Save parameters
        s <- summary(gc_fit)
        t.max <- as.numeric(t.max)
        y.max <- as.numeric(y.max)
        param_n[index,] <- c(u[jj], r[ii], q[ll])
        param[index,] <- c(t.max, y.max, s$mu.spline, s$lambda.spline,s$integral.spline)
        #param[index,] <- c(s$mu.spline, s$lambda.spline,s$integral.spline)
        index <- index + 1
        
      }
    }
  }
}


param <- as.data.frame(param)
param_n <- as.data.frame(param_n)
colnames(param_n) <- c("strain","rep","inoc")
colnames(param) <- c("tmax", "ymax", "expgr","lag","auc")

param_orig_od <- param
param$strain <- param_n$strain
param$inoc <- param_n$inoc
param$rep <- param_n$rep

dim(param) 

#Output
write.csv(param, "output/param_OD.csv")

##Correlation between CSHeat and OD growth parameters ####

#Only useful for time derived parameters, maybe exp gr? Not for y-max/value peak as different units. AUC is depended on selected data

#data manipulation
pod <- as.data.table(read.csv("output/param_OD.csv")[,-1])
pcs <- as.data.table(read.csv("output/param_CS.csv")[,-1])
pcsh <- as.data.table(read.csv("output/param_CSHeat.csv")[,-1])

#cor(pcs$time, pod$time)
#cor(pcs$ymax, pod$value) #different units
#cor(pcs$expgr, pod$exp_gr)
#cor(pcs$lag, pod$lag)

#CS heat flow vs OD
#SAC042W tmax
cor.test(pcs$tmax[pcs$strain == "SAC042W"| pcs$rep == "A"], pod$tmax[pod$strain == "SAC042W" | pod$rep == "A"])
cor.test(pcs$tmax[pcs$strain == "SAC042W"| pcs$rep == "B"], pod$tmax[pod$strain == "SAC042W" | pod$rep == "B"])
cor.test(pcs$tmax[pcs$strain == "SAC042W"| pcs$rep == "C"], pod$tmax[pod$strain == "SAC042W" | pod$rep == "C"])
cor.test(pcs$tmax[pcs$strain == "SAC042W"], pod$tmax[pod$strain == "SAC042W"])

#M116 tmax
cor.test(pcs$tmax[pcs$strain == "M116"| pcs$rep == "A"], pod$tmax[pod$strain == "M116" | pod$rep == "A"])
cor.test(pcs$tmax[pcs$strain == "M116"| pcs$rep == "B"], pod$tmax[pod$strain == "M116" | pod$rep == "B"])
cor.test(pcs$tmax[pcs$strain == "M116"| pcs$rep == "C"], pod$tmax[pod$strain == "M116" | pod$rep == "C"])
cor.test(pcs$tmax[pcs$strain == "M116"], pod$tmax[pod$strain == "M116"])

#SA2704 tmax
cor.test(pcs$tmax[pcs$strain == "SA2704"| pcs$rep == "A"], pod$tmax[pod$strain == "SA2704" | pod$rep == "A"])
cor.test(pcs$tmax[pcs$strain == "SA2704"| pcs$rep == "B"], pod$tmax[pod$strain == "SA2704" | pod$rep == "B"])
cor.test(pcs$tmax[pcs$strain == "SA2704"| pcs$rep == "C"], pod$tmax[pod$strain == "SA2704" | pod$rep == "C"])
cor.test(pcs$tmax[pcs$strain == "SA2704"], pod$tmax[pod$strain == "SA2704"])

#SAMUP15a tmax
cor.test(pcs$tmax[pcs$strain == "SAMUP15a"| pcs$rep == "A"], pod$tmax[pod$strain == "SAMUP15a" | pod$rep == "A"])
cor.test(pcs$tmax[pcs$strain == "SAMUP15a"| pcs$rep == "B"], pod$tmax[pod$strain == "SAMUP15a" | pod$rep == "B"])
cor.test(pcs$tmax[pcs$strain == "SAMUP15a"| pcs$rep == "C"], pod$tmax[pod$strain == "SAMUP15a" | pod$rep == "C"])
cor.test(pcs$tmax[pcs$strain == "SAMUP15a"], pod$tmax[pod$strain == "SAMUP15a"])

#SAC042W lag
cor.test(pcs$lag[pcs$strain == "SAC042W"| pcs$rep == "A"], pod$lag[pod$strain == "SAC042W" | pod$rep == "A"])
cor.test(pcs$lag[pcs$strain == "SAC042W"| pcs$rep == "B"], pod$lag[pod$strain == "SAC042W" | pod$rep == "B"])
cor.test(pcs$lag[pcs$strain == "SAC042W"| pcs$rep == "C"], pod$lag[pod$strain == "SAC042W" | pod$rep == "C"])
cor.test(pcs$lag[pcs$strain == "SAC042W"], pod$lag[pod$strain == "SAC042W"])

#M116 lag
cor.test(pcs$lag[pcs$strain == "M116"| pcs$rep == "A"], pod$lag[pod$strain == "M116" | pod$rep == "A"])
cor.test(pcs$lag[pcs$strain == "M116"| pcs$rep == "B"], pod$lag[pod$strain == "M116" | pod$rep == "B"])
cor.test(pcs$lag[pcs$strain == "M116"| pcs$rep == "C"], pod$lag[pod$strain == "M116" | pod$rep == "C"])
cor.test(pcs$lag[pcs$strain == "M116"], pod$lag[pod$strain == "M116"])

#SA2704 lag
cor.test(pcs$lag[pcs$strain == "SA2704"| pcs$rep == "A"], pod$lag[pod$strain == "SA2704" | pod$rep == "A"])
cor.test(pcs$lag[pcs$strain == "SA2704"| pcs$rep == "B"], pod$lag[pod$strain == "SA2704" | pod$rep == "B"])
cor.test(pcs$lag[pcs$strain == "SA2704"| pcs$rep == "C"], pod$lag[pod$strain == "SA2704" | pod$rep == "C"])
cor.test(pcs$lag[pcs$strain == "SA2704"], pod$lag[pod$strain == "SA2704"])

#SAMUP15a lag
cor.test(pcs$lag[pcs$strain == "SAMUP15a"| pcs$rep == "A"], pod$lag[pod$strain == "SAMUP15a" | pod$rep == "A"])
cor.test(pcs$lag[pcs$strain == "SAMUP15a"| pcs$rep == "B"], pod$lag[pod$strain == "SAMUP15a" | pod$rep == "B"])
cor.test(pcs$lag[pcs$strain == "SAMUP15a"| pcs$rep == "C"], pod$lag[pod$strain == "SAMUP15a" | pod$rep == "C"])
cor.test(pcs$lag[pcs$strain == "SAMUP15a"], pod$lag[pod$strain == "SAMUP15a"])

#CS heat flow vs OD
#SAC042W expgr
cor.test(pcs$expgr[pcs$strain == "SAC042W"| pcs$rep == "A"], pod$expgr[pod$strain == "SAC042W" | pod$rep == "A"])
cor.test(pcs$expgr[pcs$strain == "SAC042W"| pcs$rep == "B"], pod$expgr[pod$strain == "SAC042W" | pod$rep == "B"])
cor.test(pcs$expgr[pcs$strain == "SAC042W"| pcs$rep == "C"], pod$expgr[pod$strain == "SAC042W" | pod$rep == "C"])
cor.test(pcs$expgr[pcs$strain == "SAC042W"], pod$expgr[pod$strain == "SAC042W"])

#M116 expgr
cor.test(pcs$expgr[pcs$strain == "M116"| pcs$rep == "A"], pod$expgr[pod$strain == "M116" | pod$rep == "A"])
cor.test(pcs$expgr[pcs$strain == "M116"| pcs$rep == "B"], pod$expgr[pod$strain == "M116" | pod$rep == "B"])
cor.test(pcs$expgr[pcs$strain == "M116"| pcs$rep == "C"], pod$expgr[pod$strain == "M116" | pod$rep == "C"])
cor.test(pcs$expgr[pcs$strain == "M116"], pod$expgr[pod$strain == "M116"])

#SA2704 expgr
cor.test(pcs$expgr[pcs$strain == "SA2704"| pcs$rep == "A"], pod$expgr[pod$strain == "SA2704" | pod$rep == "A"])
cor.test(pcs$expgr[pcs$strain == "SA2704"| pcs$rep == "B"], pod$expgr[pod$strain == "SA2704" | pod$rep == "B"])
cor.test(pcs$expgr[pcs$strain == "SA2704"| pcs$rep == "C"], pod$expgr[pod$strain == "SA2704" | pod$rep == "C"])
cor.test(pcs$expgr[pcs$strain == "SA2704"], pod$expgr[pod$strain == "SA2704"])

#SAMUP15a expgr
cor.test(pcs$expgr[pcs$strain == "SAMUP15a"| pcs$rep == "A"], pod$expgr[pod$strain == "SAMUP15a" | pod$rep == "A"])
cor.test(pcs$expgr[pcs$strain == "SAMUP15a"| pcs$rep == "B"], pod$expgr[pod$strain == "SAMUP15a" | pod$rep == "B"])
cor.test(pcs$expgr[pcs$strain == "SAMUP15a"| pcs$rep == "C"], pod$expgr[pod$strain == "SAMUP15a" | pod$rep == "C"])
cor.test(pcs$expgr[pcs$strain == "SAMUP15a"], pod$expgr[pod$strain == "SAMUP15a"])

#CS heat flow vs OD
#SAC042W ymax
cor.test(pcs$ymax[pcs$strain == "SAC042W"| pcs$rep == "A"], pod$ymax[pod$strain == "SAC042W" | pod$rep == "A"])
cor.test(pcs$ymax[pcs$strain == "SAC042W"| pcs$rep == "B"], pod$ymax[pod$strain == "SAC042W" | pod$rep == "B"])
cor.test(pcs$ymax[pcs$strain == "SAC042W"| pcs$rep == "C"], pod$ymax[pod$strain == "SAC042W" | pod$rep == "C"])
cor.test(pcs$ymax[pcs$strain == "SAC042W"], pod$ymax[pod$strain == "SAC042W"])

#M116 ymax
cor.test(pcs$ymax[pcs$strain == "M116"| pcs$rep == "A"], pod$ymax[pod$strain == "M116" | pod$rep == "A"])
cor.test(pcs$ymax[pcs$strain == "M116"| pcs$rep == "B"], pod$ymax[pod$strain == "M116" | pod$rep == "B"])
cor.test(pcs$ymax[pcs$strain == "M116"| pcs$rep == "C"], pod$ymax[pod$strain == "M116" | pod$rep == "C"])
cor.test(pcs$ymax[pcs$strain == "M116"], pod$ymax[pod$strain == "M116"])

#SA2704 ymax
cor.test(pcs$ymax[pcs$strain == "SA2704"| pcs$rep == "A"], pod$ymax[pod$strain == "SA2704" | pod$rep == "A"])
cor.test(pcs$ymax[pcs$strain == "SA2704"| pcs$rep == "B"], pod$ymax[pod$strain == "SA2704" | pod$rep == "B"])
cor.test(pcs$ymax[pcs$strain == "SA2704"| pcs$rep == "C"], pod$ymax[pod$strain == "SA2704" | pod$rep == "C"])
cor.test(pcs$ymax[pcs$strain == "SA2704"], pod$ymax[pod$strain == "SA2704"])

#SAMUP15a ymax
cor.test(pcs$ymax[pcs$strain == "SAMUP15a"| pcs$rep == "A"], pod$ymax[pod$strain == "SAMUP15a" | pod$rep == "A"])
cor.test(pcs$ymax[pcs$strain == "SAMUP15a"| pcs$rep == "B"], pod$ymax[pod$strain == "SAMUP15a" | pod$rep == "B"])
cor.test(pcs$ymax[pcs$strain == "SAMUP15a"| pcs$rep == "C"], pod$ymax[pod$strain == "SAMUP15a" | pod$rep == "C"])
cor.test(pcs$ymax[pcs$strain == "SAMUP15a"], pod$ymax[pod$strain == "SAMUP15a"])

#CS heat vs OD
#cor.test(pcsh$tmax[pcsh$strain == "M116"| pcsh$rep == "A"], pod$tmax[pod$strain == "M116" | pod$rep == "A"])
#cor.test(pcsh$tmax[pcsh$strain == "M116"| pcsh$rep == "B"], pod$tmax[pod$strain == "M116" | pod$rep == "B"])
#cor.test(pcsh$tmax[pcsh$strain == "M116"| pcsh$rep == "C"], pod$tmax[pod$strain == "M116" | pod$rep == "C"])
#cor.test(pcsh$tmax[pcsh$strain == "M116"], pod$tmax[pod$strain == "M116"])
#cor(pcsh$value, pod$value) #verschillende units
#cor(pcsh$exp_gr, pod$exp_gr)
#cor(pcsh$lag, pod$lag)


