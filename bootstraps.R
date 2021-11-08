## Implementation of Bootstrapping to estimate model parameters

library(purrr)

## We want to find x that maximizes 
## g(x) = b_0 + b_1*x + b_2 * x^2
## based on some data (X,Y) ***Model Y = g(X)+e

# we have true parameters
# b_0 = -1
# b_1 = 2
# b_2 = -1
g_x <- function(x){
  return(-1 + 2*x - x^2)
}
# so true max x should be at x=1

# for bootstrapping, we take a single sample of 100 points
# we then sample 100 points from that sample with replacement
# and calculate an estimate from each of those
# we can then use these bootstrapped estimates to find a
# better confidence interval for our maximizing x
xsampboot <- runif(100, 0,2)  
ysampboot <- g_x(xsampboot) + rnorm(100, 0, .04) # original sample
xsampbootsq <- xsampboot^2
theta_fit <- lm(ysampboot ~ xsampboot + xsampbootsq)
# compute an OLS estimate from this sample
# using estimates of b_1 and b_2
theta_hat <- -theta_fit$coefficients[2]/(2 * theta_fit$coefficients[3])

# now we bootstrap on the sample taken above
bootstrap_gx <- function(n){
  i = 1
  theta_hats <- numeric(1000)
  while(i<=n){
    # randomly select 100 points from the above sample 
    # with replacement! key to note here rdunif 
    # can and will repeat outputs sometimes
    boopler <- rdunif(100,1,100)
    hat_sample_x <- xsampboot[boopler]
    hat_sample_x_sq <- hat_sample_x^2
    hat_sample_y <- ysampboot[boopler]
    # repeat process of calculating OLS theta value
    # on our bootstrap sample
    hat_linfit <- lm(hat_sample_y ~ hat_sample_x + hat_sample_x_sq)
    theta_hats[i] <- -hat_linfit$coefficients[2]/(2 * hat_linfit$coefficients[3])
    i = i + 1
  }
  return(theta_hats)
}
# use the bootstrap estimates to calculate a variance
bootstrap_thetas <- bootstrap_gx(1000)
se_sq_bootstrap <- sum((bootstrap_thetas - theta_hat)^2)/1000
# now make a 95% confidence interval for the estimate
CI_bootstrap <- c(theta_hat - 1.96 * sqrt(se_sq_bootstrap)
                  , theta_hat + 1.96 * sqrt(se_sq_bootstrap))

plot(density(sample_thetas))
abline(v = CI_bootstrap[1])
abline(v = CI_bootstrap[2])
abline(v = theta_hat)
CI_bootstrap

# now repeat the process multiple times and check how
# often our 95% CI includes the true value 1

bootstrap_repetition <- function(n){
  i=1
  count = 0
  while(i<=n){
    # take a sample and calculate OLS fit estimate
    xsamp <- runif(100, 0,2)  
    ysamp <- g_x(xsamp) + rnorm(100, 0, .04)
    xsampsq <- xsamp^2
    t_fit <- lm(ysamp ~ xsamp + xsampsq)
    t_hat <- -t_fit$coefficients[2]/(2 * t_fit$coefficients[3])
    print(t_hat)
    # bootstrap for SE estimate
    bootstraps <- bootstrap_gx(1000)
    se_sq_boot <- sum((bootstraps - t_hat)^2)/1000
    # now make a 95% confidence interval for the estimate
    CI_bootstrap <- c(t_hat - 1.96 * sqrt(se_sq_boot)
                      , t_hat + 1.96 * sqrt(se_sq_boot))
    # if the true value 1 falls outside of the confidence 
    # interval count it as a miss
    print(CI_bootstrap)
    if(CI_bootstrap[1] > 1 | CI_bootstrap[2] < 1){
      
      count = count + 1
    }
    i=i+1
    print(i)
  }
  # return proportion of trials where true value
  # falls outside of the confidence interval
  return(count/n)
  
}

bootstrap_repetition(500)
