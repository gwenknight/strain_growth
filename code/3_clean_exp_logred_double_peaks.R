#### Next steps
## Where do we go from the parameter set? 
# (1) Get the parameter sets from 2_analysis.R 

# (2) Clean datasets: 
## Remove those with exponential growth outside the set range (see 3_exponential_growth_variation.R)

# (3) Predict inoculum reduction 


#################**************** (1) Libraries, data and code needed *******************###############
library(tidyverse) 
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(ggplot2)
library(patchwork) # for combining plots
library(RColorBrewer)
library(here)
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here::here())

#####*************************** READ IN DATA *******************###############
ddm_orig <- read_csv("data_paper2/output/cut_all_time_series_fit_params.csv")[,-1]
ddm <- ddm_orig %>% filter(source == "Macotra")

param_orig <- read_csv("data_paper2/output/cut_all_model_fit_params.csv")[,-1]
param <- param_orig %>% filter(strain_name %in% unique(ddm$strain))

# Change from inoc to scalar
# scalar = 1 at 10^6, 1/10th of each next dilution
param$scalar <- 0

w<- which(param$inocl == "2")
param[w,"scalar"] = 100
w<- which(param$inocl == "3")
param[w,"scalar"] = 1000
w<- which(param$inocl == "4")
param[w,"scalar"] = 10000
w<- which(param$inocl == "5")
param[w,"scalar"] = 100000
w<- which(param$inocl == "6")
param[w,"scalar"] = 1000000

param$scalar2 <- (param$scalar)^2
strains <- unique(param$strain_name)
reps <- unique(param$rep)
drytimes <- unique(param$drytime)

# Change replicates from within experiment given names to replicate 1 to 3 for all 
po <- param %>% group_by(strain_name) %>% dplyr::mutate(maxx = max(rep), minn = min(rep), 
                                                        ones = ifelse(rep == minn, 1, 0), threes = ifelse(rep == maxx,1,0), twos = ifelse(ones == 0, ifelse(threes == 0,1,0),0)) %>%
  mutate(rep_st = case_when((ones == 1) ~ 1,
                            (threes == 1)  ~ 3,
                            (twos == 1) ~ 2)) %>% # tries pmax etc didn't work # Labels reps as 1 2 3
  dplyr::select(-c(maxx,minn, ones, threes, twos))

write.csv(po, "data_paper2/output/param_labelled_repst.csv")


#####*************************** FILTERED plot - only the clean data *******************###############
all_strains = unique(param$strain_name)

colourCount = length(unique(ddm$odd_type)) - 1
getPalette = colorRampPalette(brewer.pal(9, "Set3"))
cols = c("#000000",getPalette(colourCount))

dir.create(file.path(here(), "data_paper2/plots/final_data_split_highlighted/"),showWarnings = FALSE)

for(jj in 1:length(all_strains)){ # for each strain
  
  # Want to keep the unique combination of replicate and dataset that are clean
  clean = param %>% filter(strain_name == all_strains[jj])
  
  ddm_strain <- ddm %>% filter(strain == all_strains[jj])
  ddm_orig_s <- ddm_orig %>% filter(strain == all_strains[jj])
  
  wc <- c()
  for(i in 1:dim(clean)[1]){
    w1 <- intersect(which(ddm_strain$rep == as.numeric(clean[i,"rep"])),which(ddm_strain$inoc == as.numeric(clean[i,"inocl"])))
    wc<-c(wc,intersect(w1,which(ddm_strain$drytime == as.numeric(clean[i,"drytime"]))))
  }
  
  dd <- ddm_strain[wc,] %>% group_by(strain, inoc, rep) %>% filter(Time > 3.5)
  dd$odd_type <- as.character(dd$odd_type)
  ddm_orig_s$odd_type <- as.character(ddm_orig_s$odd_type)
  #dd$odd_type_db <- as.character(dd$odd_type_db)
  #ddm_orig_s$odd_type_db <- as.character(ddm_orig_s$odd_type_db)
  
  # 1 = peak 
  # 2 = width 
  # 3 = shoulder_before
  # 4 = shoulder_after
  # 5 = minor_peak
  
  ggplot(dd, aes(x=Time, y = value)) + 
    geom_line(aes(group = inoc, col = odd_type, linetype = factor(inoc)), lwd = 2) + 
    facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
    scale_color_manual("Odd_type", 
                       breaks = c("0", "0125", "045", "04", "02", "024", "025", 
                                  "01245", "03", 
                                  "0145", "035", "01235", 
                                  "0135", "023", "015", "0235", "0245", 
                                  "012345", "0345", "02345", "05", "034"),
                       labels = c("None","peak_width_minorpeak","shoulderafter_minorpeak","shoulderafter","width","width_shoulderafter","width_minorpeak",
                                  "peak_width_shouldafter_minorpeak","shoulderbefore",
                                  "peak_shoulderafter_minorpeak","shoulderbefore_minorpeak","peak_width_shoulderbefore_minorpeak",
                                  "peak_shoulderbefore_minorpeak","width_shoulderbefore","peak_minorpeak","width_shoulderbefore_minorpeak","width_shoulderafter_minorpeak",
                                  "all","shoulderbefore_shoulderafter_minorpeak","all_but_peak","minorpeak","shoudlerbefore_shoulderafter"),
                       values = cols, drop = FALSE) + 
    scale_linetype_discrete("Inoculum size 10^x") + 
    geom_line(data =  ddm_orig_s, aes(group = inoc, col = odd_type, linetype = factor(inoc)), alpha = 0.2, lwd = 2) + 
    geom_point(data = dd %>% filter(shoulder_point_t > 0), aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
    geom_point(data = dd %>% filter(shoulder_point_t > 0), aes(x=shoulder_point_t, y = shoulder_point_v), col = "red") + 
    geom_point(data = dd %>% filter(shoulder_point_past_t > 0), aes(x=shoulder_point_past_t, y = shoulder_point_past_v), col = "blue") + 
    geom_point(data = dd %>% filter(mp_t1 > 0), aes(x=mp_t1, y = mp_h1), col = "green") +
    geom_point(data = dd %>% filter(mp_t2 > 0), aes(x=mp_t2, y = mp_h2), col = "green") +
    geom_point(data = dd %>% filter(mp_t3 > 0), aes(x=mp_t3, y = mp_h3), col = "green") +
    geom_point(data = dd %>% filter(mp_t4 > 0), aes(x=mp_t4, y = mp_h4), col = "green") +
    labs(y = "Heat flow (mW)") +
    ggtitle(all_strains[jj]) + 
    theme(legend.key.width=unit(2,"cm"))

  
  ggsave(paste0("data_paper2/plots/final_data_split_highlighted/",all_strains[jj],"_filtered.pdf"), height = 10, width = 15) # if any to highlight it is shown here
  ggsave(paste0("data_paper2/plots/final_data_split_highlighted/",all_strains[jj],"_filtered.png"), height = 10, width = 15) # if any to highlight it is shown here
}




