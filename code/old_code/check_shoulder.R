#### Check for extra shoulder

#### THOUGHT ABUOT PUTTING THIS IN 2_analysis_cut to check for additional shoulder - driven by 11016 but not sensitive enough? 

###### SHOULDER
# (3) Is there a shoulder? 
# GIVES distance from line to curve post peak - if far from line then there is a "shoulder"
# First, draw straight line from peak
# endline = bottom of line near time zero
time_endline <- -40 #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
value_endline <- 0
## Start of line is where? (top point)
time_startline <- as.numeric(timepeak)
value_startline <- as.numeric(valpeak)

# Draw straight line, assuming peak time = time 0. 
grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line

times_line <- seq(time_startline,time_endline,-0.25) # The times for the line (x values)


# Predicted straight line
pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)

### Visualise where peaks are if needed
plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(-40,25))
points(timepeak, valpeak, col = 'black', pch = 19)
# and can plot line against this if needed too
lines(times_line, pred_points_fit, col= "blue")

### Run through lines. 
odd_shoulder<- 0
#shoulder_point <- time_peaks[1] # shoulder before this
### Visualise where peaks are if needed
plot(data1$Time,data1$value_J, type = 'l', xlab = "Time", ylab = "Value_J", xlim = c(0,25))
points(timepeak, valpeak, col = 'black', pch = 19)

for(i in -40:5){
  time_endline <- i #max(min(data1$Time)+0.5, time_peaks[1] - 10) # 10 hrs from first peak or at least 30min in to recording data
  value_endline <- 0
  #where.end <- which.min(abs(data1$Time - time_endline)) # What time exactly is this? 
  #value_endline <- data1[where.end, "value_J"] # Find value at this time point
  ## Start of line is where? (top point)
  time_startline <- as.numeric(timepeak) # highest peak
  value_startline <- as.numeric(valpeak) #as.numeric(data1[peaks_index[1],"value_J"])
  
  # Draw straight line, assuming peak time = time 0. 
  grad = as.numeric((value_endline - value_startline)/(time_endline - time_startline)) # Gradient of line
  times_line <- seq(time_startline,i,-0.25)  # The times for the line (x values)
  
  # Predicted straight line
  pred_points_fit <- grad*(times_line - time_startline) + value_startline # remove times_startline as taking this to be time zero (so can use value_startline as y intercept)
  
  
  ## and can plot line against this if needed too 
  lines(times_line, pred_points_fit, col= "blue")
  
  #### How far is the predicted straight line from the data? 
  #if(length(pred_points_fit) > 28){leng = 28}else{leng = length(pred_points_fit)}
  #peaks_index[1]:exp_start
  top <- which(data1$Time == as.numeric(timepeak))
  tail <- max(which(round(data1$Time,0) == 3))
  
  dist <- as.numeric(unlist(c(pred_points_fit[1:(top-tail+1)] - data1[top:tail,"value_J"])))
  
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

odd_shoulder

if(odd_shoulder == 0){shoulder_point <- 0; shoulder_point_v <- 0 # no shoulder_point if no shoulder!
}else{ws <- which(round(data1$Time,5) == round(shoulder_point,5)); shoulder_point_v <- data1[ws,"value_J"]}