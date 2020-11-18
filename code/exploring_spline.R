#### Exploring spline

strain <- "11277"
replicate <- 1.1
condition <- 0
inocl <- 3
data <- ddm

## which rows of the data are this strain and replicate?
wi <- intersect(which(data$strain == strain),which(data$rep == replicate)) # if fit to each replicate
wj <- intersect(wi, which(data$drytime == condition))
w <- intersect(wj, which(data$inoc == as.numeric(inocl)))

data1 <- data[w,] 

y.spl <- smooth.spline(data1$Time,data1$value_J,spar=grofit.control()$smooth.gc)

dydt.spl   <- predict(y.spl, data1$Time, deriv = 1) # derivative of smoothed spline
index      <- which.max(dydt.spl$y)          #index of maximum derivative
t.max      <- dydt.spl$x[index] # time at maximum derivative
dydt.max   <- max(dydt.spl$y) # Maximum derivative
y.max      <- y.spl$y[index] # y value at max deriv
mu.spl     <- dydt.max; # slope at maximum derivative = dydt.max
b.spl      <- y.max-dydt.max*t.max           #intercept: y value at end of lag time (y = mx + b => y - mx = b) 
lambda.spl <- -b.spl/mu.spl # 
integral   <- low.integrate(y.spl$x,y.spl$y) #
