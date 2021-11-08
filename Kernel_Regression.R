## Implementation of Kernel Density Estimation and Kernel Regression
## with bandwidth choosing using leave-one-out cross validation

# density estimate at x for a gaussian kernel and tuning h
KDE_x <- function(x,data,h){
  return(sum(dnorm((x-data)/h))/(length(data)*h))
}
# finding optimal bandwidth h using Silverman's rule
Silverman <- function(data){
  return(.9 *length(data)^(-1/5)
         * min(sqrt(var(data)), IQR(data)/1.34))
}

data = rnorm(100)
# plot our estimates for some data
KDE_plot <- function(data){
  i=1
  x = seq(from=min(data),to=max(data), length.out = 1000)
  y = numeric(1000)
  h = Silverman(data)
  while(i<=1000){
    y[i] = KDE_x(x[i],data,h)
    i=i+1
  }
  plot(x,y, type ='s')

}
# test data
data <- rnorm(100,10,3)
KDE_plot(data)
# Kernel Regression of Y on X with tuning h
# Gaussian kernel
ker_reg <- function(x_data,y_data,h,x_input){
  return(sum(dnorm((x_input-x_data)/h)*y_data)
         /sum(dnorm((x_input-x_data)/h)))
}
# compute predicted y's on a range of x values
# for plotting purposes
ker_reg_x_range<- function(x_data,y_data,h, x_min,x_max,n){
  i=1
  x <- seq(x_min,x_max,length.out = n)
  pred_y <- numeric(n)
  while(i<=n){
    pred_y[i] = ker_reg(x_data,y_data,h,x[i])
    i=i+1
  }
  return(data.frame(x,pred_y))
}

# leave one out cross validation on the data for n
# values of h in the range h_min to h_max
LOOCV <- function(x_data,y_data,n, h_min, h_max){
  len = length(x_data)
  h_vals <- seq(h_min,h_max,length.out = n)
  score <- numeric(n)
  i=1
  while(i<=n){
    tot=0
    j= 1
    while(j<=len){
      tot = tot + (y_data[j] - ker_reg(x_data[-j],y_data[-j]
                                       ,h_vals[i],x_data[j]))^2
      # regress on the data without point (x_j,y_j)
      # plug x_j into that model, then check square difference
      # with y_j for the one-left-out model
      # sum for a total score
      j=j+1
      
    }
    score[i] = tot
    i=i+1
  }
  h_scores <- data.frame(h_vals,score)
  return(h_scores)
}

# generate some random data to test on
# using the nonlinear doppler function with some normal noise
doppler <- function(x){
  return(sqrt(x*(1-x))*sin((2.1*3.14)/(x+.05)))
}
x <- seq(0,1,length.out =1000)
y <- doppler(x) + rnorm(1000,0,.1)

## get the best h value from cross validation
loocv_result <- LOOCV(x,y,100,.001,.003)
plot(loocv_result)
loocv_result$h_vals[which.min(loocv_result$score)]
# use that to inform the model
model_predicted_values <-ker_reg_x_range(x,y
                                         ,loocv_result$h_vals[which.min(loocv_result$score)]
                                         ,0,1,500)
# plot true underlying function, generated data, and fitted curve
curve(doppler,from=0,to=1,col='blue')
points(x,y)
lines(model_predicted_values,col = 'red', type ='s')
