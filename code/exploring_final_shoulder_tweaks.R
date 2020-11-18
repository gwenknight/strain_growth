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

## Run thru each strain/rep/drying time/ inoculum: fit a growth curve to each
for(jj in 1:length(u)){ # for each strain
  r <- unique(ddm %>% filter(strain == u[jj]) %>% dplyr::select(rep))[,1]
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
          
          ## What is the maximum heat flow and when?
          wmax <- which.max( data1$value_J[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
          time_max_heat_flow <- as.numeric(data1[wmax,"Time"])
          value_max_heat_flow <- as.numeric(data1[wmax,"value_J"])
          
          ## ODD 
          # (1) Is the peak broad? 
          interval_peak <- time_max_heat_flow + c(-2.5, 2.5) # time interval around peak time
          interval_value <- as.numeric(unlist(data1[c( which(round(data1$Time,4) == round(interval_peak[1],4)), which(round(data1$Time,4) == round(interval_peak[2],4))),"value_J"]))
          
          max_level <- 0;
          max_level <- max(max_level,round(100*interval_value / value_max_heat_flow,0)) # asymmetric peak then take highest
          
          # ODD? 
          if(max_level >= thresh_wide){odd_width <- 1}
          
          # (2) Are there multiple peaks? 
          # GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
          peaks_index = find_peaks(data1$value_J, m = 5)
          if(length(peaks_index) > 1){
            # REMOVE - early ones
            we<-which(data1[peaks_index,"Time"]>3) 
            # REMOVE - late ones
            w<-intersect(we,which(data1[peaks_index,"Time"]< 0.90*max(data1$Time))) # remove > 95% of time
            wl<-intersect(w,which(data1[peaks_index,"Time"] > 0.6*max(data1$Time))) # which in the odd 70-95% of the time range
            if(length(wl)>0){ # if a late peak check it is big 
              for(gg in 1:length(wl)){
                if(data1[peaks_index[wl[gg]],"value_J"] < 0.45*max(data1[peaks_index,"value_J"])){ # if not bigger than 40% of peak
                  w <-setdiff(w,wl[gg]) }}}   # then remove
            peaks_index <- peaks_index[w] # Keep the ones not at the beginning / end
            # Sort by height - want to compare to and keep the tallest (first now in peaks_index)
            o <- order(data1[peaks_index,"value_J"], decreasing = "TRUE")
            peaks_index <- peaks_index[o]
          }
          
          # When are the peaks? 
          time_peaks <- as.numeric(unlist(data1[peaks_index,"Time"]))
          time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
          # If multiple far apart then issue: double peaks
          keep_time_far_apart <- which(abs(time_peaks_diff) > 5) 
          
          if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
          if(any(data1[peaks_index,"value_J"] < 0.001)){odd_peak <- 1} # or if peak low
          
          # If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
          close_peaks <- which(abs(time_peaks_diff) <= 5)
          close_peaks_i <- 1
          if(length(close_peaks)>1){
            for(i in 2:length(close_peaks)){
              ifelse(data1[peaks_index[close_peaks[i]],"value_J"]/data1[peaks_index[close_peaks[1]],"value_J"]> 0.9,
                     close_peaks_i <- c(close_peaks_i,i),"")
            }}
          
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
          
          
        }
      }
    }
  }
}

