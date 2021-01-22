#### Function

### Input: times series data
### Output: time series data up to cut point: time of end of exponential growth
### Output: lag time / exponential growth / time of end of exponential growth 

cut_extract <- function(ts, strain, replicate, condition, inocl, thresh_wide = 80, plot = 0){
  ## ts = timeseries (time and value_J)
  ## strain, replicate, condition, inocl = labels for output
  ## plot = 0:  don't plot. 1 to plot
  ## thresh_wide = 80: % what is a wide peak? 
  
  
  
  ## is this strain replicate ODD? Set all ODD indicators to zero initially
  odd_peak <- 0 # 0 = one peak
  odd_width <- 0 # 0 = not a wide peak
  odd_shoulder <- 0 # 0 = no shoulder
  max_level <- 0 # height of peak
  odd_double <- 0 # double curve fits this data
  
  ## This gives lag time and exponential growth rate cumulative
  gc_fit <- gcFitSpline(ts$Time, ts$value_J)
  # parameters from this fit
  s <- summary(gc_fit)
  
  ## What is the maximum heat flow and when?
  wmax <- which.max( ts$value_J[-c(1:6)]) + 6 # take off funny initial behaviour < 2hr (take off in which.max and then add on 6 for position as taken 6 off)
  time_max_heat_flow <- as.numeric(ts[wmax,"Time"])
  value_max_heat_flow <- as.numeric(ts[wmax,"value_J"])
  
  ## ODD 
  # (1) Is the peak broad? 
  interval_peak <- time_max_heat_flow + c(-2.5, 2.5) # time interval around peak time
  interval_value <- as.numeric(unlist(ts[c( which(round(ts$Time,4) == round(interval_peak[1],4)), which(round(ts$Time,4) == round(interval_peak[2],4))),"value_J"]))
  
  max_level <- 0;
  max_level <- max(max_level,round(100*interval_value / value_max_heat_flow,0)) # asymmetric peak then take highest
  
  # ODD? 
  if(max_level >= thresh_wide){odd_width <- 1}
  
  # (2) Are there multiple peaks? 
  # GIVES WHERE ALL PEAKS ARE (greater or equal to (m) 5 points around them)
  peaks_index = find_peaks(ts$value_J, m = 5)
  if(length(peaks_index) > 1){
    # REMOVE - early ones
    we<-which(ts[peaks_index,"Time"]>3) 
    # REMOVE - late ones
    w<-intersect(we,which(ts[peaks_index,"Time"]< 0.90*max(ts$Time))) # remove > 95% of time
    wl<-intersect(w,which(ts[peaks_index,"Time"] > 0.6*max(ts$Time))) # which in the odd 70-95% of the time range
    if(length(wl)>0){ # if a late peak check it is big 
      for(gg in 1:length(wl)){
        if(ts[peaks_index[wl[gg]],"value_J"] < 0.45*max(ts[peaks_index,"value_J"])){ # if not bigger than 40% of peak
          w <-setdiff(w,wl[gg]) }}}   # then remove
    peaks_index <- peaks_index[w] # Keep the ones not at the beginning / end
    # Sort by height - want to compare to and keep the tallest (first now in peaks_index)
    o <- order(ts[peaks_index,"value_J"], decreasing = "TRUE")
    peaks_index <- peaks_index[o]
  }
  
  # When are the peaks? 
  time_peaks <- as.numeric(unlist(ts[peaks_index,"Time"]))
  time_peaks_diff <- time_peaks - time_peaks[1] # how far apart are they?
  # If multiple far apart then issue: double peaks
  keep_time_far_apart <- which(abs(time_peaks_diff) > 5) 
  
  if(length(keep_time_far_apart) >= 1){odd_peak <- 1} # if multiple peaks
  if(any(ts[peaks_index,"value_J"] < 0.001)){odd_peak <- 1} # or if peak low
  
  # If close and same height (90% of tallest) then odd (peak decline plateau decline OK)
  close_peaks <- which(abs(time_peaks_diff) <= 5)
  close_peaks_i <- 1
  if(length(close_peaks)>1){
    for(i in 2:length(close_peaks)){
      ifelse(ts[peaks_index[close_peaks[i]],"value_J"]/ts[peaks_index[close_peaks[1]],"value_J"]> 0.9,
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
  time_endline <- -40 #max(min(ts$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
  value_endline <- 0
  ## Start of line is where? (top point)
  time_startline <- time_peaks[1] # highest peak
  value_startline <- as.numeric(ts[peaks_index[1],"value_J"])
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  # Only want exponential growth line
  exp_start <- which.min(abs(ts$Time-s$lambda.spline+1))
  times_line <- as.numeric(unlist(ts[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
  
  
  # Predicted straight line
  pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
  
  ### Visualise where peaks are if needed
  # plot(ts$Time,ts$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(-40,25))
  #  points(unlist(ts[peaks_index,"Time"]), unlist(ts[peaks_index,"value_J"]), col = 'black', pch = 19)
  ## and can plot line against this if needed too 
  #  lines(times_line, pred_points_fit, col= "blue")
  
  ### Run through lines. 
  odd_shoulder<- 0
  shoulder_point <- time_peaks[1] # shoulder before this
  ### Visualise where peaks are if needed
  # plot(ts$Time,ts$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(0,25))
  # points(unlist(ts[peaks_index,"Time"]), unlist(ts[peaks_index,"value_J"]), col = 'black', pch = 19)
  
  for(i in -40:5){
    time_endline <- i #max(min(ts$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
    value_endline <- 0
    #where.end <- which.min(abs(ts$Time - time_endline)) # What time exactly is this? 
    #value_endline <- ts[where.end, "value_J"] # Find value at this time point
    ## Start of line is where? (top point)
    time_startline <- time_peaks[1] # highest peak
    value_startline <- as.numeric(ts[peaks_index[1],"value_J"])
    
    # Draw straight line, assuming peak time = time 0. 
    grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
    times_line <- as.numeric(unlist(ts[peaks_index[1]:exp_start,"Time"])) # The times for the line (x values)
    
    # Predicted straight line
    pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
    
    
    ## and can plot line against this if needed too 
    # lines(times_line, pred_points_fit, col= "blue")
    
    #### How far is the predicted straight line from the data? 
    #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
    dist <- c(pred_points_fit - as.numeric(unlist(ts[peaks_index[1]:exp_start,"value_J"])))
    
    # Crossing points
    if(length(which(abs(diff(sign(dist)))==2)) > 2){# if cross more than twice then shoulder
      if(max(-dist) > 0.001){ # far enough away to be a proper shoulder
        if(max(dist) > 0.001){
          #if(times_line[which.max(-dist)] < (time_peaks[1] - 2)){ # and not just another peak
          fpsq <- find_peaks(as.numeric(unlist(dist)), m = 3) # Are there more than 1 peaks?
          if(length(fpsq) != 0){ # of no preaks then not a clear shoulder
            shoulder_point1 = min(time_peaks[1]-1,max(times_line[fpsq]))
            ws <- which(round(ts$Time,5) == round(shoulder_point1,5)); shoulder_point_v1 <- ts[ws,"value_J"]
            height <- shoulder_point_v1 / ts[peaks_index[1],"value_J"]
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
  }else{ws <- which(round(ts$Time,5) == round(shoulder_point,5)); shoulder_point_v <- ts[ws,"value_J"]}
  
  #### If no shoulder but multiple peaks, want to grab time of first peak
  if(length(time_peaks_diff) > 1 & odd_shoulder == 0){ # if multiple peaks but no shoulder
    shoulder_point1 = min(time_peaks)
    shoulder_point_v1 = ts[which(round(ts$Time,5) == round(shoulder_point1,5)),"value_J"]
    if((shoulder_point_v1 / ts[peaks_index[1],"value_J"]) > 0.5){ # shoulder needs to be high still! 
      shoulder_point = shoulder_point1; shoulder_point_v = shoulder_point_v1 # Update shoulder values to ones that are > half way
    }
  }
  
  ##### Double curves? 
  x <- ts$Time
  y <- ts$value_J
  
  # Guesses for the parameters in the curve
  startl=list(a=unlist(ts[peaks_index,"value_J"])[1]*1.5, 
              b=diff(interval_peak)/3, 
              c=unlist(ts[peaks_index,"Time"])[1],
              d=unlist(ts[peaks_index,"value_J"])[1]*1.5, 
              e=diff(interval_peak), 
              f=unlist(ts[peaks_index,"Time"])[1])
  # Set indicator/output parameters to zero
  fit1 <- 0
  plot_p <- 0
  p <- 0
  # Try to fit two normal curves using non-linear lear squares methods 
  try(fit1 <- nls(y~(a/b)*exp(-(x-c)^2/(2*b^2))+(d/e)*exp(-(x-f)^2/(2*e^2)),start = startl), silent = TRUE)
  
  # If manage to fit
  if(length(fit1) > 1){
    # print that can 
    print(paste("Double curve fit", strain," ", replicate, condition, inocl, sep = " "))
    odd_double <- 1
    # Save the predictions of the fit 
    pred_p <- predict(fit1)
    # Save the coefficients of the fit
    p <- coef(fit1)
    # Pull out the actual curves using the parameters in p
    g1p <- (p["a"]/p["b"])*exp(-(x-p["c"])^2/(2*p["b"]^2)) 
    g2p <- (p["d"]/p["e"])*exp(-(x-p["f"])^2/(2*p["e"]^2))
    # Store the curves
    plot_p <- as.data.frame(cbind(x,y,pred_p, g1p, g2p))
    colnames(plot_p) <- c("time","value_J","fit","normal_curve1","normal_curve2")
    # Store the parameters
    plot_p$a <- p["a"]; plot_p$b <- p["b"]; plot_p$c <- p["c"]
    plot_p$d <- p["d"]; plot_p$e <- p["e"]; plot_p$f <- p["f"]
  }
  
  ## Plot ts, cumulative and fit
  # if functions takes in a command to plot
  if(plot == 1){
    #### FIT TO Csum? 
    # datam <- melt(ts[,c("Time","value_J","csum")], id.vars = "Time")
    # gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() + facet_wrap(~strain, scales = "free") 
    # ## model fit:
    # gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
    # colnames(gc_df) <- c("Time","csum")
    # gc_df$value_J <- c(gc_df[1,"csum"],diff(gc_df$csum))
    # gc_dfm <- melt(gc_df, id.vars = "Time")
    # ## add fit to data plot
    # gg <- gg + geom_line(data = gc_dfm, aes(x=Time,y=value), col= "red")
    # 
    
    ### Fit to value_J : could add cumulative plot in here but not done
    datam <- reshape2::melt(ts[,c("Time","value_J")], id.vars = "Time")
    gg <- ggplot(datam, aes(x=Time,y=value)) + geom_point() #+ facet_wrap(~strain, scales = "free") 
    ## model fit:
    gc_df <- as.data.frame(cbind(gc_fit$fit.time, gc_fit$fit.data))
    colnames(gc_df) <- c("Time","value")
    #gc_df$value_J <- c(gc_df[1,"csum"],diff(gc_df$csum)) # CHANGE TO CUMULATIVE? 
    gc_dfm <- reshape2::melt(gc_df, id.vars = "Time")
    ## add fit to data plot
    gg <- gg + geom_line(data = gc_df, aes(x=Time,y=value), col= "red")
    
    ## save
    ggsave(paste0("plots/",strain,"_",replicate, "_",condition,"_", inocl, "_model_fit.pdf"))
  } 
  
  ## Build vectors of required parameters to output
  param_o   <- c(time_max_heat_flow, value_max_heat_flow, 
                 s$mu.spline, s$lambda.spline,s$integral.spline, 
                 odd_peak, odd_width, max_level, odd_shoulder, odd_double, shoulder_point, shoulder_point_v)
  
  
  ####******************************************* This is the first analysis. *******************************************
  # Now cut at shoulder or first peak
  if(shoulder_point > 0){
    cut_point_t <- shoulder_point; cut_point_v <- shoulder_point_v}else{
      cut_point_t <- t_m_h_flow; cut_point_v <- v_m_h_flow
    }
  
  # NEW TS: cut up to first peak or shoulder 
  ts <- ts %>% 
    filter(Time > 3) %>% # cut off first 3hrs 
    mutate(cutpart = ifelse(Time <= cut_point_t,1,0)) %>%
    filter(cutpart == 1) # Trim off extra parts
  
  
  ### Look at this cut data
  ## Currently the time to use is: peak growth 
  timepeak = cut_point_t
  valpeak = cut_point_v[1]
  
  ## What if there is a peak in this?
  peaks_index = find_peaks(ts$value_J, m = 3)
  
  ## If there is a peak then reassign peak 
  if(length(peaks_index) > 0){
    w_early <- which(ts[peaks_index,"Time"] < 5) # Remove early ones
    if(length(w_early) > 0){peaks_index <- peaks_index[-w_early]}
    w_high <- which(ts[peaks_index,"value_J"] < 0.6*max(ts$value_J)) # Remove low ones
    if(length(w_high) > 0){peaks_index <- peaks_index[-w_high]}
    if(length(peaks_index) > 0){ # if any later than 3 
      peaks_index = min(peaks_index)
      timepeak = ts[peaks_index,"Time"]
      valpeak = ts[peaks_index,"value_J"]
    }
  }
  
  # Reassign up to peak 
  ts = ts[which(ts$Time <= as.numeric(timepeak)),]
  
  ### Check if the slope changes substantially in this period: may be a shoulder or a plateau near peak 
  if(dim(ts)[1]>4){ # If enough data
    st <- c()
    
    #plot(ts$Time, ts$value_J, "l")
    #points(ts$Time, ts$value_J)
    
    for(i in 2:(-1 + dim(ts)[1])){ # fit a linear model to the data in segments
      ts1 <- ts[max(1,i-6):i,]
      ts2 <- ts[(i):dim(ts)[1],]
      lm.1 <- lm(ts1$value_J ~ ts1$Time)
      lm.2 <- lm(ts2$value_J ~ ts2$Time)
      
      #lines(ts1$Time, lm.1$coefficients[1] + lm.1$coefficients[2]*ts1$Time, col = "blue")
      #lines(ts2$Time, lm.2$coefficients[1] + lm.2$coefficients[2]*ts2$Time, col = "red")
      
      st <- rbind(st, c(i,c(max(ts1$value_J), max(ts1$Time),lm.1$coefficients[2],lm.2$coefficients[2])))
    }
    
    st <- as.data.frame(st)
    colnames(st) <- c("i","maxval","maxtime","f_ang","s_ang") # first angle, second angle
    st$d <- 0; st$da <- 1000
    st$d[1:(dim(st)[1]-1)] = diff(st$s_ang) # look at change in second angle: want to know when substantial change
    st$da[1:(dim(st)[1]-1)] = abs(diff(st$s_ang)) # look at change in second angle: want to know when substantial change
    
    #plot(st$maxtime, st$d) 
    #plot(ts$Time, ts$value_J)
    
    ## Check if shoulder
    st_upper <- st%>% filter(maxtime > 0.50*max(st$maxtime)) # but want to be past halfway
    w1 <- which.min(st_upper$da) # as absolute this detects a plateau 
    w2 <- which.min(st_upper$da[-w1])
    timepeak_s = min(st_upper[w1,"maxtime"], st_upper[-w1,"maxtime"][w2]) # earliest 
    valpeak_s = as.numeric(ts[which(ts$Time == timepeak_s),"value_J"])
    
    ## Check not too low: if cut point already good enough. 
    # If OK then cut at shoulder value
    if(valpeak_s > 0.65*max(ts$value_J)){
      valpeak <- valpeak_s; timepeak = timepeak_s}else{
        
        ## Check if end peak should be moved forward at all (i.e. a peak: want end of exponential growth )
        st_upper <- st %>% filter(maxval > 0.80*max(st$maxval)) # want to be near the end 
        #if(dim(st_upper)[1] > 5){ # have at least 4 d values to look at
        w1 <- which.max(st_upper$d)
        w2 <- which.max(st_upper$d[-w1])
        if(min(w1,w2) != 1){ # if its not just a slope down - if there is a plateau near the top? i.e. index now just the earliest point
          timepeak = st_upper[min(w1,w2),"maxtime"]
          valpeak = as.numeric(ts[which(ts$Time == timepeak),"value_J"])}
      }
  }
  
  
  g1 <- ggplot(ts,aes(x=Time, y= value_J)) + geom_line() + 
    geom_point(data = ts[which(round(ts$Time,2) == as.numeric(round(timepeak,2))),c("Time","value_J")],col="red") + 
    geom_point(data= ts[1,], aes(x=shoulder_point_t, y = shoulder_point_v), col = "black") + 
    ggtitle(paste0(u[jj],"_",r[ii], "_",drying_times[kk],"_",q[ll]))
  ggsave(paste0("plots/shoulder_curves/cutpoint_highlighted_",strain,"_",replicate, "_",condition,"_",inocl,".pdf")) 
  
  tstopeak = ts[which(ts$Time <= as.numeric(timepeak)),c("Time","value_J")]
  
  ## Growth Curve # (see fig 3 of vv33i07.pdf)
  ## This gives lag time and exponential growth rate cumulative
  gc_fit <- gcFitSpline(tstopeak$Time, tstopeak$value_J)
  # parameters from this fit
  s <- summary(gc_fit)
  
  
  ## Build vectors of required parameters to output
  param_o   <- c(param_o, s$mu.spline, timepeak, valpeak)
  
  return(list(param = param_o))
  
}