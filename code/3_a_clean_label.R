#### Next steps
## Where do we go from the parameter set? 
# (1) Get the parameter sets from 2_analysis.R then 2_analysis_cut.R

# (2) Clean datasets: 
### How many odd datasets in a drytime for a single rep? 
##### If < 50% of datasets in a rep at a drytime are odd then remove just this dataset
##### If > 50% of datasets are odd in one drytime then label the rep as odd
### How many odd replicates are there for a strain?
##### If more than one replicate is odd then strain is odd
##### If only one then remove the replicate from the strain
### > NOW just label for "removal" but as only use up to peak OK to keep in

# (3) Check exponential growth 
### Should be same across replicate

#################**************** (1) Libraries, data and code needed *******************###############
library(tidyverse) # good library for filtering / summarising data
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
library(ggplot2)
library(patchwork) # for combining plots
theme_set(theme_bw(base_size=24)) # theme setting for plots: black and white (bw) and font size (24)

setwd(here:here())

##### READ IN DATA
ddm <- read.csv("output/cut_all_ddm.csv")[,-1]
param <- read.csv("output/cut_all_param.csv")[,-1]

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

### Code for linear model 
source("code/function_linear_model.R")

#################**************** (2) LABEL ODD STRAINS *******************###############

####### How many odd? 
w_odd <- which(param$any_odd>=1) # any_odd is the sum of all odd indicators. 1= just one indicator, 3 = three indicators etc
length(w_odd)
100*length(w_odd) / dim(param)[1] # 16% of datasets
# e.g. 1_2_6_ has 19% of datasets are odd
100*length(unique(param[w_odd,"strain_name"])) / length(unique(param$strain_name)) # 54% of strains have an odd dataset

### By looking at the distribution of "odd" behaviour in strains can we see those that are non-typical? Not clear
pdf("plots/distribution_oddness_in_strains.pdf")
hist(table(param[w_odd,"strain_name"]), breaks = seq(0,15,1))
dev.off()

t <- table(param[w_odd,"strain_name"])
odd_by_freq <- names(t[order(t)][(dim(t)-10):dim(t)]) # 10 with the greatest number of odd? 

#### (1) LABEL replicate if > 50% odd ( at least 3 replicates per strain). For removal?
#### (2) LABEL dataset if < 50% odd (usually at least 3 datasets per replicate and drying time). For removal?
remove_reps <- c()
remove_dataset <- c()

for(i in strains){ # for all strains
  for(j in reps){ # and replicates
    print(c(i,j))
    # Which data is for this strain and rep? 
    wthis <- intersect(which(param$strain_name == i),which(param$rep == j))
    
    # Within a dry time (e.g. 168hrs) what is odd? 
    for(k in drytimes){
      wk <- intersect(wthis, which(param$drytime == k))
      if(length(wk)>0){
        # how many odd datasets? 
        n_odd <- nnzero(param[wk,"any_odd"])
        # how many datasets?
        n_measures <- length(wk)
        # proportion odd
        p_odd_m <- n_odd / n_measures 
        if(p_odd_m >= 0.5){remove_reps <- rbind(remove_reps, c(i, j))} # if more than 50% write it down to remove this rep
        if(n_odd > 0 & p_odd_m < 0.5){ # if there are odd ones but proportion < 50% just remove this dataset
          where_odd <- which(param[wk,"any_odd"] > 0) # which are odd
          for(iww in where_odd){
            remove_dataset <- rbind(remove_dataset, c(i,j,k, as.numeric(param[wk[iww],"inocl"])))} # store this dataset to remove
        }
      }
    }
  }
}

# (1) Check if any reps to remove but keep strain typical
tr <- table(remove_reps[,1]) # How many replicates to remove? 
w1 <- which(tr == 1) # Only remove a replicate if there is only one to remove
# Cycle through those to remove
if(length(w1)>0){
  strain_remove_one <- names(tr)[w1] # which have only one to remove? 
  rep_to_remove <- c()
  # For those strains that have one replicate to remove
  for(i in strain_remove_one){
    wi <- which(remove_reps[,1] == i) # which is the replicate?
    rep_to_remove <- rbind(rep_to_remove,
                           cbind(rep(i,length(wi)), remove_reps[which(remove_reps[,1] == i),2]))
  }
  wremove <- c()
  # Index of where this rep is to remove
  for(i in 1:length(rep_to_remove[,1])){
    wremove <- c(wremove, 
                 intersect(which(param$strain_name == rep_to_remove[i,1]),
                           which(param$rep == rep_to_remove[i,2])))
  }
}

# LABEL those reps that should be removed 
param$removed_rep <- 0 # Note this
param[wremove,"removed_rep"] <- 1


# (2) Check if any datasets to remove but keep strain typical
t <- table(remove_dataset[,1]) # How many datasets to remove? 

## Label all these datasets
param$removed_dataset <- 0
for(i in 1:length(remove_dataset[,1])){
  wst <- which(param$strain_name == remove_dataset[i,1])
  wstr <- intersect(wst, which(param$rep == remove_dataset[i,2]))
  wstrt <- intersect(wstr, which(param$drytime == remove_dataset[i,3]))
  wstrti <- intersect(wstrt, which(param$inocl == remove_dataset[i,4]))
  
  param[wstrti,"removed_dataset"] <- 1
}


# which have more than one replicate that should be removed? ODD STRAINS
w1 <- which(tr > 1) # Which have more than one odd replicate? 
odd_strains_names <- names(tr[w1])

param <- param %>% mutate(odd_strains = ifelse(strain_name %in% odd_strains_names, 1, 0))

typical_strains <- as.numeric(unlist(param %>% filter(odd_strains == 0) %>% summarise(unique(strain_name))))

### STORE 
write.csv(param, paste0("output/param_labelled.csv"))

#### Move plots of odd strains into one folder
dir.create(paste0("plots/odd"), showWarnings = FALSE) # warnings = FALSE otherwise get warning if folder already exists
dir.create(paste0("plots/typical"), showWarnings = FALSE)


name_code <- "all_(1_13)_"
for(i in odd_strains){
  file.copy(paste0("plots/",name_code,"odd_highlighted_",i,".pdf"), paste0("plots/odd"))
}
for(i in typical_strains){
  file.copy(paste0("plots/",name_code,"odd_highlighted_",i,".pdf"), paste0("plots/typical"))
}

ddm$odd_type <- as.character(ddm$odd_type)
write.csv(ddm, "output/data_to_cut.csv")

#### FILTERED plot - only the clean data
### ALL STRAINS 
all_strains = unique(param$strain)

ddm_orig <- read.csv(paste0("output/",name_code,"_all_ddm.csv"))[,-1] # Whole time course data

cols = c(1,brewer.pal(n = 8, name = "Set1"))

for(jj in 1:length(all_strains)){ # for each strain
  
  # Wnat to keep the unique combinatino of replicate and dataset that are clean
  clean = param #%>% filter(removed_rep == 0, removed_dataset == 0, strain_name == all_strains[jj])
  
  ddm_strain <- ddm %>% filter(strain == all_strains[jj])
  ddm_orig_s <- ddm_orig %>% filter(strain == all_strains[jj])
  
  wc <- c()
  for(i in 1:length(clean[,1])){
    w1 <- intersect(which(ddm_strain$rep == clean[i,"rep"]),which(ddm_strain$inoc == clean[i,"inocl"]))
    wc<-c(wc,intersect(w1,which(ddm_strain$drytime == clean[i,"drytime"])))
  }
  
  dd <- ddm_strain[wc,]
  ddm_orig_s$odd_type <- as.character(ddm_orig_s$odd_type)
  
  ggplot(dd, aes(x=Time, y = value_J)) + 
    geom_line(aes(group = inoc, col = odd_type, linetype = factor(inoc)), lwd = 1) + 
    facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
    scale_color_manual("Odd_type", breaks = c("0","1","2","3","12","13","23","123"), 
                       labels = c("None","Peak","Width","Shoulder","Peak&Width","Peak&Shoulder",
                                  "Width&Shoulder","Peak Width&Shoulder"),
                       values = cols, drop = FALSE) + 
    scale_linetype_discrete("Inoc.") + 
    geom_line(data =  ddm_orig_s, aes(group = inoc, col = odd_type, linetype = factor(inoc)), alpha = 0.2, size = 1) + 
    geom_point(aes(x=cut_timepeak, y = cut_valpeak), col = "red") + 
    ggtitle(all_strains[jj])
    
  ggsave(paste0("plots/final_data_split_highlighted/",all_strains[jj],"_filtered.pdf")) # if any to highlight it is shown here
}


#################**************** (3) CHECK EXPONENTIAL GROWTH*******************###############
param_clean <- param #%>% filter(removed_rep == 0, removed_dataset == 0) ### COULD FILTER ON BEHAVIOUR PAST PEAK 

# Look at the distribution of exponential growth for typical strains
for(i in unique(param_clean$strain_name)){
  pa <- param_clean %>% filter(strain_name == i)
  g1 <- ggplot(pa, aes(x=rep, y = cut_exp, group = rep)) + geom_boxplot() + ggtitle(paste0(i," all"))+ scale_x_continuous("Replicate") + scale_y_continuous("Exponential growth") 
  
  g2 <- ggplot(pa, aes(x=rep, y = cut_exp, group = interaction(rep,drytime))) + geom_boxplot(aes(col = factor(drytime))) + ggtitle(paste0("Split by drying time")) + 
    scale_color_discrete("Dry time") + scale_x_continuous("Replicate") + scale_y_continuous("Exponential growth") + 
    theme(legend.position="bottom")
  
  g1 + g2 + plot_layout(widths = c(1, 2)) 
  if(pa$odd_strains == 0){ggsave(paste0("plots/exp_growth/typical_",i,".pdf"))}else{ggsave(paste0("plots/exp_growth/odd_",i,".pdf"))}
}

po <- param_clean %>% group_by(strain_name) %>% dplyr::mutate(maxx = max(rep), minn = min(rep), 
                                                       ones = ifelse(rep == minn, 1, 0), threes = ifelse(rep == maxx,1,0), twos = ifelse(ones == 0, ifelse(threes == 0,1,0),0)) %>%
  mutate(rep_st = case_when((ones == 1) ~ 1,
                            (threes == 1)  ~ 3,
                            (twos == 1) ~ 2)) # tries pmax etc didn't work

write.csv(po, "output/param_labelled_repst.csv")

##### Diagnostic / exploring plots
ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_replicates.pdf")

ggplot(po, aes(x=inocl, y = cut_exp, group = interaction(inocl, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_inocl_bp.pdf")

ggplot(po, aes(x=inocl, y = cut_exp, group = interaction(inocl, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_inocl_bpgrp.pdf")

ggplot(po, aes(x=inocl, y = cut_exp,aes(group = drytime))) + geom_point() +  
  geom_smooth(method = "loess") +  #, #, formula = y ~ a * x + b,method.args = list(start = list(a = 0.1, b = 0.1))) + 
  facet_wrap(drytime~odd_strains) + 
  scale_y_continuous(lim = c(0,0.1),"Exponential growth") + 
  scale_x_continuous("Inoculum") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/exp_growth_across_inocl_pointsline.pdf")

ggplot(po, aes(x=rep_st, y = t_m_h_flow, group = interaction(rep_st, drytime, strain_name))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous("Time to peak") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/time_to_peak_across_replicates.pdf")

ggplot(po, aes(x = cut_exp, y = t_m_h_flow,group = interaction(rep_st, drytime, strain_name))) + 
  geom_point(aes(col = factor(drytime))) + 
  geom_line() +
  facet_wrap(~odd_strains)


ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_jitter(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous(lim = c(0,0.1),"exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/across_replicates_zoom.pdf")

ggplot(po, aes(x=interaction(rep_st, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
  geom_jitter(width = 0.1,aes(col = factor(drytime))) + 
  facet_wrap(~strain_name, scales = "free_x") + 
  scale_y_continuous("Exponential growth") + 
  scale_x_discrete("Replicate/Strain") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/exp_growth/across_replicates_zoom_strain_rep.pdf",width = 30, height = 30)

ggplot(po, aes(x=interaction(rep_st, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
  geom_jitter(width = 0.1,aes(col = factor(inocl), pch = factor(drytime))) + 
  facet_wrap(~strain_name, scales = "free_x") + 
  scale_y_continuous("Exponential growth") + 
  scale_x_discrete("Replicate/Strain") +   
  scale_shape_discrete("Dry time") + 
  scale_color_manual("Inoculum", values = c("red","black","turquoise")) + 
  theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/exp_growth/across_replicates_zoom_strain_rep_inocl.pdf",width = 30, height = 30)

ggplot(po, aes(x=interaction(inocl, strain_name), y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + 
  geom_jitter(width = 0.1,aes(col = factor(inocl), pch = factor(drytime))) + 
  facet_wrap(~strain_name, scales = "free_x") + 
  scale_y_continuous("Exponential growth") + 
  scale_x_discrete("Inoculum/Strain") +   
  scale_shape_discrete("Dry time") + 
  scale_color_manual("Inoculum", values = c("red","black","turquoise")) + 
  theme(legend.position="bottom") + theme(axis.text.x = element_text(angle = 90))
ggsave("plots/exp_growth/across_replicates_zoom_strain_rep_inoclgrp.pdf",width = 30, height = 30)


ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime, strain_name))) + geom_jitter(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous("exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/across_replicates_all.pdf")

ggplot(po, aes(x=rep_st, y = cut_exp, group = interaction(rep_st, drytime))) + geom_boxplot(aes(col = factor(drytime))) + 
  facet_wrap(~odd_strains) + 
  scale_y_continuous(lim = c(0,0.1),"exponential growth") + 
  scale_x_continuous("Replicate") +   scale_color_discrete("Dry time") + 
  theme(legend.position="bottom")
ggsave("plots/exp_growth/group_replicates.pdf")

ggplot(po %>% filter(drytime == 0), aes(x=inocl, y = cut_exp, group = interaction(rep,strain))) + 
  geom_line(aes(col = strain))

#### Normalised start
ddm <-ddm %>% group_by(strain,rep, inoc, drytime) %>% dplyr::mutate(initial = first(value_J), value_J_norm = value_J - initial)
ggplot(ddm, aes(x = Time, y = value_J_norm, group= interaction(rep, inoc, drytime))) + 
  geom_line(aes(col = factor(drytime))) + facet_wrap(~strain)
ggsave("plots/exp_growth/normalised_start.pdf", width = 20, height = 20)


#### Example 
# 11272
d2 <- ddm %>% filter(strain == "11277", rep == 1.1) %>% ungroup() %>% select(Time,value_J)
s <- summary(gcFitSpline(d2$Time, d2$value_J))

# ##### Check exp growth OK across inocula of the same strain at set drying times
# perc <- 0.2
# param_exp_gr_lab <- param %>%
#   group_by(strain_name, rep) %>%  # Mean over strain_name and rep - want to be same over dry times
#   mutate(mean_peak_exp_gr = mean(exp_gr),
#          mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
#          mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
#          outside = ifelse(exp_gr < mean_peak_exp_gr_m10 | exp_gr > mean_peak_exp_gr_p10, 1, 0)) %>%
#   dplyr::mutate(sum_outside_s_r = sum(outside)) # number of odd datasets in this strain and replicate
# 
# #### check for exp_gr - remove if only one is "wrong"
# param_typical <- param_exp_gr_lab %>%
#   filter(sum_outside_s_r < 2) # some have only one in the strain and replicate
# 
# length(unique(param_typical$strain_name))
# # For 1_2_6_ 9 strains still have typical behaviour
# 
# strains_typical = unique(param_typical$strain_name) # PERFECT strains
# 
# dim(param_typical_notexp_gr)[1] - dim(param_typical)[1] # 79 datasets removed because of exp_gr - one of those from the one strain that is removed ("RN6390B")
# 
# ggplot(param_exp_gr_lab[1:100,], aes(x=inocl, y = exp_gr)) + geom_point(aes(colour = factor(outside))) +
#   scale_color_manual("In the limits?", values = c("black","red")) +
#   facet_wrap(strain_name ~ rep + drytime) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_m10)) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr_p10)) +
#   geom_hline(aes(yintercept = mean_peak_exp_gr),lty = "dashed") +
#   ggtitle(paste0(100*perc," percent from mean exponential growth"))
# ggsave(paste0("plots/exp_growth/",name_code,100*perc,"perc_from_mean_exponential_growth.pdf"),width = 20, height = 20)
# 
# param %>% filter(strain_name == "11277", rep == 1.1, drytime == 0)
