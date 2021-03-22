#FIGURE FOR HEAT FLOW VARIATION

library(ggplot2)
library(patchwork)
theme_set(theme_bw(base_size=14)) # theme setting for plots: black and white (bw) and font size (24)

# DATASET 1 TYPICAL 
strain <- "Newman" 
replicate <- 44 
condition <- 0
inocl <- 5

####### DATA
data <- ddm_orig
wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
w <- intersect(wi, which(data$inoc == as.numeric(inocl)))

data_orig <- data[w,]  # just get the data for this experiment (time, value, )
data1 <- subset(data_orig, drytime == condition) 

h1 <- ggplot(data1, aes(x = Time, y = value)) + 
  geom_line() +
  scale_x_continuous("Time (h)", minor_breaks = seq(0, 25, 1), limits = c(0,25)) + 
  scale_y_continuous(expression(paste("Heatflow (mW)")) , limits = c(-0.005, 0.09), breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  labs(title = expression(paste("Newman (",10^5, ") before dehydration"))) 

# DATASET 2 TYPICAL 
strain <- "SAC042W" 
replicate <- 31 
condition <- 0
inocl <- 5

####### DATA
data <- ddm_orig
wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
w <- intersect(wi, which(data$inoc == as.numeric(inocl)))

data_orig <- data[w,]  # just get the data for this experiment (time, value, )
data2 <- subset(data_orig, drytime == condition) 

h2 <- ggplot(data2, aes(x = Time, y = value)) + 
  geom_line() +
  scale_x_continuous("Time (h)", minor_breaks = seq(0, 25, 1), limits = c(0,25)) + 
  scale_y_continuous(expression(paste("Heatflow (mW)")) , limits = c(-0.005, 0.09), breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  labs(title = expression(paste("SAC042W (",10^5, ") before dehydration"))) 

# DATASET 3 DOUBLE PEAK 
strain <- "11016" 
replicate <- 4.3 
condition <- 0
inocl <- 5

####### DATA
data <- ddm_orig
wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
w <- intersect(wi, which(data$inoc == as.numeric(inocl)))

data_orig <- data[w,]  # just get the data for this experiment (time, value, )
data3 <- subset(data_orig, drytime == condition) 

h3 <- ggplot(data3, aes(x = Time, y = value)) + 
  geom_line() +
  scale_x_continuous("Time (h)", minor_breaks = seq(0, 25, 1), limits = c(0,25)) + 
  scale_y_continuous(expression(paste("Heatflow (mW)")) , limits = c(-0.005, 0.09), breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  labs(title = expression(paste("11016 (",10^5, ") before dehydration"))) 

# DATASET 4 SHOULDER PEAK 
strain <- "11283" 
replicate <- 2.1 
condition <- 0
inocl <- 5

####### DATA
data <- ddm_orig
wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
w <- intersect(wi, which(data$inoc == as.numeric(inocl)))

data_orig <- data[w,]  # just get the data for this experiment (time, value, )
data4 <- subset(data_orig, drytime == condition) 

h4 <- ggplot(data4, aes(x = Time, y = value)) + 
  geom_line() +
  scale_x_continuous("Time (h)", minor_breaks = seq(0, 25, 1), limits = c(0,25)) + 
  scale_y_continuous(expression(paste("Heatflow (mW)")) , limits = c(-0.005, 0.09), breaks = c(0, 0.02, 0.04, 0.06, 0.08)) +
  labs(title = expression(paste("11283 (",10^5, ") before dehydration"))) 

#ALL PLOTS
(h1 + h2) / (h3 + h4) + plot_layout(guides = 'collect', widths = c(1,2)) + plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom')
ggsave("plots/final/figure3.pdf", width = 20, height = 15)
