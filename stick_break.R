## Implementation of stick breaking to sample from a 
## Dirichlet Process

## Stick breaking with tuning alpha and N breaks
stick_break<- function(N, alpha){
  i = 1
  sticks <- c()
  stick_remain <- 1
  while(i<N){
    brk <- rbeta(1,1,alpha)
    sticks <- append(sticks,stick_remain * brk)
    stick_remain = stick_remain* (1-brk)
    i = i +  1
  }
  sticks <- append(sticks, stick_remain)
  return(sticks)
}


#stick breaking with standard normal as base distribution
norm_sample <- function(N, alpha){
  # sample N points from the normal distribution
  # as the locations of the sticks
  x_locs <- rnorm(N)
  # next perform the stick breaking process
  # then assign each x location its corresponding stick size
  y_locs <- stick_break(N, alpha)
  together <- data.frame(x_locs, y_locs)
  # order things by x location
  together <- together[order(together$x_locs),]
  return(together)
  
}
# construct a cdf from a given discrete pmf
# sampled from the Dirichlet Process DP(Norm(0,1), alpha)
cdf_builder <- function(N, alpha){
  mypdf <- norm_sample(N, alpha)
  x <- mypdf[,1]
  y <- mypdf[,2]
  total_sum <- 0
  cdf <- numeric(N)
  i = 1
  while(i<=N){
    total_sum = total_sum + y[i]
    cdf[i] <- total_sum
    i = i + 1
  }
  pt <- data.frame(x, cdf)
  return(pt)
  
}

# take k samples from DP(Norm(0,1),alpha), plot the
# sampled cdf, then compare against the actual CDF of
# Norm(0,1)
plotting_compare_norm <- function(N,k, alpha){
  i = 1
  tseq <- seq(-4,4,.001)
  plot(cdf_builder(N,alpha), type = 's')
  while(i<=k){
    lines(cdf_builder(N,alpha), type = 's')
    i = i + 1
    
  }
  lines(tseq, pnorm(tseq),col= 'red')
}


## calculating and plotting empirical CDFs for some data
# alongside confidence bands from DKW inequality

# construct an empirical CDF for some set of data
sample_cdf <- function(data){
  
  i = 1
  # count how many points fall at or below a point x
  # then normalize by total number of points
  count = 0
  xdat <- data[order(data)]
  cdf <- numeric(length(xdat))
  while(i<=length(xdat)){
    count = count + 1/length(xdat)
    cdf[i] <- count
    i = i + 1
    
  }
  
  return(data.frame(xdat,cdf))
}

## upper and lower bounds based on DKW inequality
## applicable to any empirically determined distribution 
## function. Can be used to create confidence bands for
## our estimates
u_bound <- function(data, n, alpha){
  u_b <- data + sqrt(1/(2*n) * log(2/alpha))
  u_b[u_b>1] <- 1
  return(u_b)
}

l_bound <- function(data, n, alpha){
  l_b <- data - sqrt(1/(2*n) * log(2/alpha))
  l_b[l_b<0] <- 0
  return(l_b)
}

q5_e <- function(n){
  a <- sample_cdf(rnorm(n, 5, 3))
  a_cdf <- a[,2]
  u <- u_bound(a_cdf, n, .95)
  l <- l_bound(a_cdf, n, .95)
  u_b <- data.frame(a[,1], u)
  l_b <- data.frame(a[,1], l)
  plot(a, type = 's')
  lines(u_b, type = 's', col = 'red')
  lines(l_b, type = 's', col = 'red')
  print(l_b)
  
}

# now construct estimates of a distribution based on 
# Bayesian updating using the stick breaking process

# we'd like to construct a posterior distribution based on 
# the data, which should take the form 
# DP(alpha* H/(alpha + n) + P_n/(alpha+n),alpha + n)
# with P_n a random sample from the data points.
posterior_sample <- function(n, alpha, data, sticks){
  q = rdunif(1, 1, n)
  b = rnorm(1)
  rand <- numeric(sticks)
  i = 1
  while(i<=sticks){
    q = rdunif(1, 1, n) # randomly select a data point
    b = rnorm(1) # prior is standard normal
    rand[i] <- (n*data[q] + alpha * b)/(n+alpha)
    i = i + 1
  }
  return(rand)
}
# now perform the stick breaking procedure again, using the
# randomly drawn point from the posterior base distribution
# as our x locations
posterior_stick_breaking <- function(sticks, alpha, n, data){
  x_locs <- posterior_sample(n, alpha, data, sticks)
  y_locs <- stick_break(sticks, alpha + n)
  together <- data.frame(x_locs, y_locs)
  together <- together[order(together$x_locs),]
  return(together)
  
}
# make our pdf from stick breaking into a cdf
posterior_cdf <- function(sticks, alpha, n, data){
  mypdf <- posterior_stick_breaking(sticks, alpha, n, data)
  x <- mypdf[,1]
  y <- mypdf[,2]
  total_sum <- 0
  cdf <- numeric(sticks)
  i = 1
  while(i<=sticks){
    total_sum = total_sum + y[i]
    cdf[i] <- total_sum
    i = i + 1
  }
  pt <- data.frame(x, cdf)
  return(pt)
}
## perform on some real data and compare with actual CDF

# 1) drawn from normal dist with shifted mean and different
#    variance
print('Norm mean 5, var 3 compare')
test_data <- rnorm(50,5,3)
post_cdf <- posterior_cdf(500, 5, 50,test_data)
curve(pnorm(x,5,3), from = -2, to = 12)
lines(post_cdf, type = 's', col = 'red')

# 2) drawn from uniform dist
print('Uniform dist compare')
test_data <- runif(50,-2,12)
post_cdf <- posterior_cdf(500, 5, 50,test_data)
curve(punif(x,-2,12), from = -2, to = 12)
lines(post_cdf, type = 's', col = 'red')

# 3) drawn from exponential dist
print('Exponential dist compare')
test_data <- rexp(100, 1)
post_cdf <- posterior_cdf(500, 1, 100,test_data)
curve(pexp(x,1), from = 0, to = 5)
lines(post_cdf, type = 's', col = 'red')
