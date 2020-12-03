#### Explore / tweak final shoulder curves
### 11050 / 11280 / 11283 issues

# EXPLORE BEHAVIOUR FOR ONE STRAIN 
# DATASET

u <- c("11280", "11283","11050")
## How many replicates? 
r <- unique(ddm$rep) # replicates
# What are the inoculums? 
q <- unique(ddm$inoc)

drying_times <- c(0,24,168)
ddm <- ddm %>% ungroup()
parah <- c()

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unlist(unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1])
  for(ii in 1:length(r)){ # for each replicate: fit to all the data, not just each replicate
    for(kk in c(1,3)){ #each of the three experimental conditions (0, 24, 168): most just 0 168 now
      for(ll in 1:length(q)){ #each of the inocula
        
        print(c(jj,ii,kk,ll))
        
        strain <- u[jj];
        replicate <- r[ii]
        condition <- drying_times[kk]
        inocl <- q[ll]
        data <- ddm
        plot = 0
        thresh_wide = 80
        
        ## which rows of the data are this strain and replicate?
        wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
        wj <- intersect(wi, which(data$drytime == condition))
        w <- intersect(wj, which(data$inoc == as.numeric(inocl)))
        
        ## is this strain replicate ODD? Set all ODD indicators to zero initially
        odd_peak <- 0 # 0 = one peak
        odd_width <- 0 # 0 = not a wide peak
        odd_shoulder <- 0 # 0 = no shoulder
        max_level <- 0 # height of peak
        odd_double <- 0 # double curve fits this data
        
        if(length(w) > 0){ # if this replicate exists for this strain (i.e. there is data)
          data1 <- data[w,]  # just get the data for this experiment (strain, time, value, drying time)
          
          ## Growth Curve # (see fig 3 of vv33i07.pdf)
          ## This gives lag time and exponential growth rate cumulative
          gc_fit <- gcFitSpline(data1$Time, data1$value_J)
          # parameters from this fit
          s <- summary(gc_fit)
          

          
          if(length(close_peaks_i) > 1){odd_peak <- 1} # if multiple close time and height peaks
          
          # Only keep those peaks that are far apart
          time_peaks <- time_peaks[c(1,keep_time_far_apart)]
          peaks_index <- peaks_index[c(1,keep_time_far_apart)]
          
          
          ###### SHOULDER
          # (3) Is there a shoulder? 
          # GIVES distance from line to curve post peak - if far from line then there is a "shoulder"
          # First, draw straight line from peak
          # endline = bottom of line near time zero
          time_endline <- -40 #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
          value_endline <- 0
          ## Start of line is where? (top point)
          time_startline <- time_peaks[1] # highest peak
          value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
          
          # Draw straight line, assuming peak time = time 0. 
          grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
          # Only want exponential growth line
          exp_start <- which.min(abs(data1$Time-s$lambda.spline+1))
          times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
          
          
          # Predicted straight line
          pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
          
          ### Visualise where peaks are if needed
          plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(-40,25))
          points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
          # and can plot line against this if needed too
          lines(times_line, pred_points_fit, col= "blue")
          
          ### Run through lines. 
          odd_shoulder<- 0
          shoulder_point <- time_peaks[1] # shoulder before this
          ### Visualise where peaks are if needed
          plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(0,25))
          points(unlist(data1[peaks_index,"Time"]), unlist(data1[peaks_index,"value_J"]), col = 'black', pch = 19)
          
          for(i in -40:5){
            time_endline <- i #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
            value_endline <- 0
            #where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
            #value_endline <- data1[where.end, "value_J"] # Find value at this time point
            ## Start of line is where? (top point)
            time_startline <- time_peaks[1] # highest peak
            value_startline <- as.numeric(data1[peaks_index[1],"value_J"])
            
            # Draw straight line, assuming peak time = time 0. 
            grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
            times_line <- as.numeric(unlist(data1[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
            
            # Predicted straight line
            pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
            
            
            ## and can plot line against this if needed too 
            lines(times_line, pred_points_fit, col= "blue")
            
            #### How far is the predicted straight line from the data? 
            #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
            dist <- as.numeric(unlist(c(pred_points_fit - data1[peaks_index[1]:exp_start,"value_J"])))
            
            # Crossing points
            if(length(which(abs(diff(sign(dist)))==2)) > 2){# if cross more than twice then shoulder
              if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
                if(max(dist) > 0.001){
                  #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
                  fpsq <- find_peaks(as.numeric(unlist(dist)), m = 3) # Are there more than 1 peaks?
                  if(length(fpsq) != 0){ # of no preaks then not a clear shoulder
                    shoulder_point1 = min(time_peaks[1]-1,max(times_line[fpsq]))
                    ws <- which(round(data1$Time,5) == round(shoulder_point1,5)); shoulder_point_v1 <- data1[ws,"value_J"]
                    height <- shoulder_point_v1 / data1[peaks_index[1],"value_J"]
                    if(height > 0.5 && height < 0.96){ # if the shoulder is greater than half way up but not too close
                      shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
                      odd_shoulder <- 1;# print(c("new",i)) # Then odd shoulder  
                    }
                  }
                }
              }
            }
            
          }
          if(odd_shoulder == 0){shoulder_point <- 0; shoulder_point_v <- 0 # no shoulder_point if no shoulder!
          }else{ws <- which(round(data1$Time,5) == round(shoulder_point,5)); shoulder_point_v <- data1[ws,"value_J"]}
          
          #### If no shoulder but multiple peaks, want to grab time of first peak
          if(length(time_peaks_diff) > 1 & odd_shoulder == 0){ # if multiple peaks but no shoulder
            shoulder_point1 = min(time_peaks)
            shoulder_point_v1 = data1[which(round(data1$Time,5) == round(shoulder_point1,5)),"value_J"]
            if((shoulder_point_v1 / data1[peaks_index[1],"value_J"]) > 0.5){ # shoulder needs to be high still! 
              shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
            }
          }
          
          parah   <- rbind(parah,c(as.numeric(strain), replicate, condition, inocl, time_max_heat_flow, value_max_heat_flow, 
                         s$mu.spline, s$lambda.spline,s$integral.spline, 
                         odd_peak, odd_width, max_level, odd_shoulder, odd_double, shoulder_point, shoulder_point_v))
        }
      }
    }
  }
}
parah <- as.data.frame(parah)
colnames(parah) <- c("strain", "rep","drytime","inocl",
                     "t_m_h_flowc", "v_m_h_flow", "exp_gr","lag","auc",
                     "odd_peaks","odd_width","width_peak","odd_shoulder","odd_double","shoulder_point_t","shoulder_point_v")
parah$shoulder_point_t <- as.numeric(unlist(parah$shoulder_point_t))
parah$shoulder_point_v <- as.numeric(unlist(parah$shoulder_point_v))

dd <- ddm %>% filter(strain %in% u)
pp <- parah %>% filter(strain %in% u) %>% filter(shoulder_point_v > 0)

ggplot(dd, aes(x=Time, y = value_J, group = interaction(rep, strain, inoc, drytime))) + 
  geom_line(aes(col = odd_type, linetype = factor(inoc))) + 
  facet_wrap(drytime~strain + rep) + 
  scale_color_manual("Odd_type", breaks = c("0","1","2","3","12","13","23","123"), 
                     labels = c("None","Peak","Width","Shoulder","Peak&Width","Peak&Shoulder",
                                "Width&Shoulder","Peak Width&Shoulder"),
                     values = seq(1,8,1), drop = FALSE) + 
  scale_linetype_discrete("Inoc.") + 
  ggtitle(paste0(u[jj]," plotted:",Sys.Date())) + 
  geom_point(data = pp, aes(x=shoulder_point_t, y =shoulder_point_v), col = "red")

