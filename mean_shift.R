## Implementation of Mean-Shift Algorithm for finding modes 
## of the density for a 1D sample drawn from
## some distribution p

# Iteration step: with sp as starting point, we iterate by the
# weighted average of densities from all data points
# h hyperparameter for tuning the kernel regression
next_step <- function(data, sp, h){
  # divide by total density contribution at sp
  den <- sum(dnorm((sp - data)/h)) 
  # numerator given by sum of data point*density contribution at sp
  numerator <- sum(data*dnorm((sp - data)/h))
  return(numerator/den)
}


# iterate n times using previous function
mean_shifts <- function(data, sp, h, n){
  i = 1
  while( i< n+1){
    i = i + 1
    sp = next_step(data, sp, h)
  }
  return(sp)
}

## helper function to "cluster" starting points that end at same
## mode together into one mode. 
breakup_helper <- function(set, distance){
  i=2
  final_set <- c(set[1])
  while(i<=length(set)){
    # loop through each point
    # if it is far enough away from points already identified
    # as modes, then add it to the list, otherwise do nothing
    if(min(abs(final_set-set[i]))> distance){
      final_set <- append(final_set,set[i])
    }
    i= i + 1
  }
  return(final_set)
}

## find modes by iterating on k evenly spaced points in the range of
## the data. First make a vector of starting points (sp_set)
## then apply the mean_shifts to each n times.
mode_finding <- function(data, n, k, dist){
  sp_set <- seq(min(data), max(data),(max(data)-  min(data)) /k)
  dens = density(data)
  bw = dens$bw
  i = 1
  while(i<=k){
    sp_set[i] <- mean_shifts(data, sp_set[i], bw, n)
    i = i + 1
  }
  # use breakup_helper 
  modes <- breakup_helper(sp_set,dist )
  # plot density of the data from kernel regression alongside 
  # the mode estimations
  plot(dens)
  i = 1
  while(i<=length(modes)){
    abline(v = modes[i])
    i = i + 1
  }
  print(modes)
  return(modes)
}

# test data with multiple modes
# mixture of 3 normal distributions
i = 1
test_data <- numeric(500)
while(i<=500){
  x = runif(1)
  if(x<.25){
    test_data[i] = rnorm(1,-1,0.1)
  }
  else if(x<.5){
    test_data[i] = rnorm(1,1,0.1)
  }
  else if(x<.75){
    test_data[i] = rnorm(1,4,0.1)
  }
  else{
    test_data[i] = rnorm(1,8,0.1)
  }
  i=i+1
}

mode_finding(test_data,50,25,0.5)
