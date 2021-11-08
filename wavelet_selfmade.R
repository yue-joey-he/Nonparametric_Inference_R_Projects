## Implementation of Wavelet Regression

# will only be for datasets of 128 length for a 
# datapoint for each Haar function up to j=7


# Haar Mother wavelet
psi_f <- function(x){
  if(x<1/2 && x >= 0){
    return(-1)
  }
  else if(x>=1/2 && x<1){
    return(1)
  }
  else{
    return(0)
  }
}

## Haar functions that will form a basis on L_2(0,1)
## indexed by j,k with 0 <= k < 2^j
psi_jk <- function(j,k,x){
  return(2^(j/2) * psi_f(2^j * x - k))
}


# Project our observed values onto a 
# a Haar function psi_jk
# these will be estimates of coefficients for each function
Z_jk <- function(j,k, x_data, y_data){
  i = 1
  s = numeric(128)
  while(i<=128){
    s[i] <- psi_jk(j,k,x_data[i])
    i = i + 1
  }
  return(1/128 * sum(y_data * s))
}



# denoise our data somewhat by soft thresholding the 
# coefficients. Return the soft thresholded values
soft_thresh <- function(a, lambda){
  i = 1
  while(i<=length(a)){
    if(a[i] < -lambda){
      a[i] <- a[i] + lambda
    }
    else if(a[i] > lambda){
      a[i] <- a[i] - lambda
    }
    else{
      a[i] <- 0
    }
    i= i + 1
  }
  return(a)
}
# These thresholded values will make up the actual coefficients
# of our fit
beta_jk <- function(j,k,x_data, y_data, lambda){
  return(soft_thresh(Z_jk(j,k, x_data, y_data), lambda))
}

# fit the model to some data
wave_predict <- function(x_data, y_data, x){
  a_hat = mean(y_data)
  # choose lambda for smoothing based on MAD variance 
  # estimate on y data
  med_res <- y_data - median(y_data)
  MAD <- median(abs(med_res))
  sig_hat <-  MAD/.6745
  lambda <- sig_hat * sqrt(2*log(128)/128)
  
  # for x value where we want to make a prediction,we find the
  # coefficient beta_jk(x) for each j,k pair we want in the fit
  # then sum over beta_jk(x)*psi_jk(x) to get our fitted y value
  # for that x
  j = 0
  tot = 0
  while(j <= 6){
    k = 0
    while(k <= 2^j - 1){
      tot = tot + beta_jk(j,k, x_data, y_data, lambda)* psi_jk(j,k,x)
      k = k + 1
    }
    j = j + 1
  }
  # final values will be centered around 0,so we shift them by
  # the average of the observed y values
  a_hat = .75 * mean(y_data)
  return(a_hat + tot)
}
# try on some data
flies_data <- read.table('flies.dat', header =TRUE)
days <- flies_data$day[1:128]
days_scaled <- days/length(flies_data$day)
mortality_rate <- flies_data$mort.rate[1:128] + .01
log_mrate <- log(mortality_rate)

i = 1
predicted_vals = numeric(128)
while(i<=128){
  predicted_vals[i] <- wave_predict(days_scaled,log_mrate
                               ,days_scaled[i])
  i = i + 1
}

plot(days_scaled, log_mrate)
lines(days_scaled, predicted_vals, col = 'green')


