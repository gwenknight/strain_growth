#### Next steps
## Where do we go from the parameter set? 
# (1) Get the parameter sets from 2_analysis_set1&2&6.R
# (2) For those that are “odd”, check if odd behaviour the same across different drying times within the same replicate. 
# If they are the same then use the data. BRANCH 1 (see below)

library(tidyverse) # good library for filtering / summarising data
library(matrixStats) # for row sd calculations
library(Matrix) # for nnzero function
theme_set(theme_bw(base_size=6)) # theme setting for plots: black and white (bw) and font size (24)

##### READ IN DATA
### CHOOSE IN LINE WITH DATA CLEANING IN 2_analysis.R
#name_code <- "1_2_6_7_"
name_code <- "all_(1_13)_"

param <- read.csv(paste0("output/",name_code,"all_model_fit_params.csv"), stringsAsFactors = FALSE)[,-1]
param <- as.data.frame(param)

# Change from inoc to scalar
# scalar = 1 at 10^6, 1/10th of each next dilution
param$scalar <- 0

w<- which(param$inocl == "2")
param[w,"scalar"] = 10
w<- which(param$inocl == "3")
param[w,"scalar"] = 100
w<- which(param$inocl == "4")
param[w,"scalar"] = 1000
w<- which(param$inocl == "5")
param[w,"scalar"] = 10000
w<- which(param$inocl == "6")
param[w,"scalar"] = 100000

param$scalar2 <- (param$scalar)^2
strains <- unique(param$strain_name)
reps <- unique(param$rep)
drytimes <- unique(param$drytime)

### Code for linear model 
source("20_function_linear_model.R")

###################************** CHECK EXPONENTIAL GROWTH OK 
##### Check exp growth OK across inocula of the same strain at set drying times
perc <- 0.2
param_exp_gr_lab <- param %>%
  group_by(strain_name, rep) %>%  # Mean over strain_name and rep - want to be same over dry times
  mutate(mean_peak_exp_gr = mean(exp_gr),
         mean_peak_exp_gr_p10 = mean_peak_exp_gr + perc*mean_peak_exp_gr,
         mean_peak_exp_gr_m10 = mean_peak_exp_gr - perc*mean_peak_exp_gr,
         outside = ifelse(exp_gr < mean_peak_exp_gr_m10 | exp_gr > mean_peak_exp_gr_p10, 1, 0)) %>%
  dplyr::mutate(sum_outside_s_r = sum(outside)) # number of odd datasets in this strain and replicate

#### check for exp_gr - remove if only one is "wrong" 
param_typical <- param_exp_gr_lab %>%
  filter(sum_outside_s_r < 2) # some have only one in the strain and replicate

length(unique(param_typical$strain_name)) 
# For 1_2_6_ 9 strains still have typical behaviour

strains_typical = unique(param_typical$strain_name) # PERFECT strains

dim(param_typical_notexp_gr)[1] - dim(param_typical)[1] # 79 datasets removed because of exp_gr - one of those from the one strain that is removed ("RN6390B")

ggplot(param_exp_gr_lab[1:100,], aes(x=inocl, y = exp_gr)) + geom_point(aes(colour = factor(outside))) +
  scale_color_manual("In the limits?", values = c("black","red")) +
  facet_wrap(strain_name ~ rep + drytime, scales = "free") +
  geom_hline(aes(yintercept = mean_peak_exp_gr_m10)) +
  geom_hline(aes(yintercept = mean_peak_exp_gr_p10)) +
  geom_hline(aes(yintercept = mean_peak_exp_gr),lty = "dashed") +
  ggtitle(paste0(100*perc," percent from mean exponential growth"))
ggsave(paste0("plots/",name_code,100*perc,"perc_from_mean_exponential_growth.pdf"),width = 20, height = 20)










################################## **** ODD TYPICAL ASSIGNATION 

####### How many odd? 
w_odd <- which(param$any_odd>=1) # any_odd is the sum of all odd indicators. 1= just one indicator, 3 = three indicators etc
length(w_odd) # 145 datasets for 1_2_6_7
100*length(w_odd) / dim(param)[1] 
# e.g. 1_2_6_ has 19% of datasets are odd
length(unique(param[w_odd,"strain_name"])) 
# e.g. 1_2_6_ has 17 odd datasets
### By looking at the distribution of "odd" behaviour in strains can we see those that are non-typical? 
pdf("plots/distribution_oddness_in_strains.pdf")
hist(table(param[w_odd,"strain_name"]), breaks = seq(0,15,1))
dev.off()
t <- table(param[w_odd,"strain_name"])
odd_by_freq <- names(t[order(t)][50:76]) # 5 and greater

#### (1) REMOVE replicate if > 50% odd ( at least 3 replicates per strain)
#### (2) REMOVE dataset if < 50% odd (usually at least 3 datasets per replicate and drying time)
remove_reps <- c()
remove_dataset <- c()

for(i in strains){
  for(j in reps){
    print(c(i,j))
    # Which data is for this strain and rep? 
    wthis <- intersect(which(param$strain_name == i),which(param$rep == j))
    
    # Within a dry time (e.g. 24hrs) what is odd? 
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

# REMOVE those reps that should be removed 
param$removed_rep <- 0 # Note this
param[wremove,"removed_rep"] <- 1
param_no_odd <- param[-wremove,] # parameter sets minus the reps to remove


# (2) Check if any datasets to remove but keep strain typical
t <- table(remove_dataset[,1]) # How many datasets to remove? 

## Remove all these datasets
for(i in 1:length(remove_dataset[,1])){
  wst <- which(param_no_odd$strain_name == remove_dataset[i,1])
  wstr <- intersect(wst, which(param_no_odd$rep == as.numeric(remove_dataset[i,2])))
  wstrt <- intersect(wstr, which(param_no_odd$drytime == remove_dataset[i,3]))
  wstrti <- intersect(wstrt, which(param_no_odd$inocl == remove_dataset[i,4]))
  if(length(wstrti) > 0){ # some replicates may already have been removed
    param_no_odd <- param_no_odd[-wstrti,]
  }
}

## Label all these datasets
param$removed_dataset <- 0
for(i in 1:length(remove_dataset[,1])){
  wst <- which(param$strain_name == remove_dataset[i,1])
  wstr <- intersect(wst, which(param$rep == remove_dataset[i,2]))
  wstrt <- intersect(wstr, which(param$drytime == remove_dataset[i,3]))
  wstrti <- intersect(wstrt, which(param$inocl == remove_dataset[i,4]))
  
  param[wstrti,"removed_dataset"] <- 1
}


### param_no_odd is all parameter data minus
# (1) datasets that could be removed
# (2) reps that could be removed 

# which have more than one replicate that should be removed? ODD STRAINS
w1 <- which(tr > 1) # Only remove a replicate if there is only one to remove
odd_strains <- names(tr[w1])

param_no_odd <- param_no_odd %>%
  filter(!strain_name %in% odd_strains)

typical_strains <- unique(param_no_odd$strain_name)

length(typical_strains)
length(odd_strains)

intersect(odd_strains, odd_by_freq)
length(intersect(odd_strains, odd_by_freq)) # all in odd by freq

### STORE 

write.csv(param_no_odd, paste0("plots/",name_code,"param_typical.csv"))

param_odd <- param %>%
  filter(strain_name %in% odd_strains)

write.csv(param_odd, paste0("plots/",name_code,"param_non_typical.csv"))

# Store with new information on what was removed
param <- mutate(param, typical = ifelse(strain_name%in% odd_strains, 0, 1))
write.csv(param, paste0("plots/",name_code,"information_all_model_fit_params.csv"))

#### Move plots of odd strains into one folder
dir.create(paste0("plots/",name_code,"odd"))
dir.create(paste0("plots/",name_code,"typical"))

for(i in odd_strains){
  file.copy(paste0("plots/",name_code,"odd_highlighted_",i,".pdf"), paste0("plots/",name_code,"odd"))
}
for(i in typical_strains){
  file.copy(paste0("plots/",name_code,"odd_highlighted_",i,".pdf"), paste0("plots/",name_code,"typical"))
}

#### FILTERED plot - only the clean data
for(jj in 1:length(typical_strains)){ # for each strain
  
  # Wnat to keep the unique combinatino of replicate and dataset that are clean
  clean = subset(param_no_odd, strain_name == typical_strains[jj])
  
  ddm_strain <- ddm %>% filter(strain == typical_strains[jj])
  #reps_clean = unique(subset(param_no_odd, strain_name == typical_strains[jj])$rep)
  #inoc_clean = unique(subset(param_no_odd, strain_name == typical_strains[jj])$inocl)
  wc <- c()
  for(i in 1:length(clean[,1])){
    w1 <- intersect(which(ddm_strain$rep == clean[i,"rep"]),which(ddm_strain$inoc == clean[i,"inocl"]))
    wc<-c(wc,intersect(w1,which(ddm_strain$drytime == clean[i,"drytime"])))
  }
  
  dd <- ddm_strain[wc,]
  
  ggplot(dd, aes(x=Time, y = value_J)) + 
    geom_line(aes(group = inoc, col = odd_type, linetype = factor(inoc))) + 
    facet_wrap(drytime~rep, nrow = length(unique(dd$drytime))) + 
    scale_color_manual("Odd_type", breaks = c("0","01","02","03","012","013","023","0123"), 
                       labels = c("None","Peak","Width","Shoulder","Peak&Width","Peak&Shoulder",
                                  "Width&Shoulder","Peak Width&Shoulder"),
                       values = seq(1,8,1), drop = FALSE) + 
    scale_linetype_discrete("Inoc.") #+ 
  #geom_text(data = pp, aes(label = squared_dist, x = 10+as.numeric(inocl), y =as.numeric(inocl)*0.001, col = factor(inocl)),  size = 2)
  ggsave(paste0("plots/",name_code,"typical/",name_code,"odd_highlighted_",typical_strains[jj],"_filtered.pdf")) # if any to highlight it is shown here
}
