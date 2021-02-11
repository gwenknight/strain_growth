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
ggsave("plots/final/time_peak_vs_inoculum_strain_col.pdf", width = 12)


ga <- ggplot(param %>% filter(strain_name %in% strains), aes(x = inocl, y = timepeak, group = inocl)) + 
  geom_boxplot() + 
  scale_x_continuous(breaks = seq(2:7), labels = function(x) parse(text=paste("10^",x)),"Inoculum", limits = c(1.5,6.5)) + 
  scale_y_continuous(expression(paste("Time to first peak (", italic(tmax),")")))
ggsave("plots/final/time_peak_vs_inoculum.pdf")


## Data
data_od <- read_csv("proof_of_principle/ddm_OD.csv")[,-1]
data_cs <- read_csv("proof_of_principle/ddm_CS.csv")[,-1]


### Look at data
g1 <- ggplot(data_cs %>% filter(strain == "SA2704"), aes(x=Time, y = value_J, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Heat curve") +
  scale_x_continuous(lim = c(0,25)) + 
  scale_color_discrete("Inoculum")

### Smooth OD values
data_od <- data_od %>% group_by(rep, exp, variable, strain, inoc) %>%
  mutate(ma_value = rollapply(value, 10, mean,fill = NA),
         differ = c(0,diff(ma_value))) %>%
  ungroup()

g2 <- ggplot(data_od %>% filter(strain == "SA2704"), aes(x=Time, y = ma_value, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("OD")+ 
  scale_color_discrete("Inoculum")

g3 <- ggplot(data_od %>% filter(strain == "SA2704"), aes(x=Time, y = differ, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("OD")+ 
  scale_color_discrete("Inoculum")
ggsave("proof_of_principle/OD_data_difference.pdf")

### Cumulate heat flow values
g4 <- ggplot(data_cs %>% filter(strain == "SA2704"), aes(x=Time, y = csum, group = interaction(rep, inoc))) + 
  geom_line(aes(col = factor(inoc))) + 
  facet_wrap(~strain) + 
  scale_y_continuous("Heat curve") + 
  scale_x_continuous(lim = c(0,25)) + 
  scale_color_discrete("Inoculum") 

gb <- ((g1 + g4) / (g3 + g2)) + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A') & theme(legend.position='bottom')
gb |  ga  
ggsave("plots/final/linear_OD_vs_CS.pdf",width = 15, height = 10)