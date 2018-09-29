# Many local maxima
## Log-likelihood function 
x <- c(3.91, 4.85, 2.28, 4.06, 3.70, 4.04, 5.46, 3.53, 2.28, 1.96,
       2.53, 3.88, 2.22, 3.47, 4.82, 2.46, 2.99, 2.54, 0.52)
log_llh <- function (theta, sample_X) {
  log_llh <- 0
  for (i in 1:length(sample_X)) {
  log_llh <- log_llh + log(1-cos(sample_X[i] - theta)) - log(2*pi)
  }
  log_llh
}
## plot
library(ggplot2)
ggplot() + stat_function(aes(-pi:pi), fun = log_llh, args = list(sample_X = x)) +
  labs(title = expression(paste("Log-likelihood function of ", theta)),
       x = expression(theta), y = "Value of log-likelihood function") +
  theme(plot.title = element_text(hjust = 0.5))

## MOM estimation of theta
theta_mom <- asin(mean(x) - pi)
theta_mom

## MLE for theta with initial value being estimation of theta by MOM
### First derivative of loglikelihood
first_derv <- function(theta, sample_X) {
  first_derv <- 0
  for (i in 1:length(sample_X)){
    first_derv <- first_derv - 
      (sin(sample_X[i] - theta)/(1 - cos(sample_X[i] - theta)))
  }
  first_derv
}

### Second derivative of loglikelihood
second_derv <- function(theta, sample_X) {
  second_derv <- 0
  for (i in 1:length(sample_X)){
    second_derv <- second_derv + (1/(cos(sample_X[i] - theta) - 1))
  }
  second_derv
}

### Newton–Raphson method
newton <- function(init, pre=.Machine$double.neg.eps, maxrun=200) {
  n <- 1
  xt <- init
  while (n < maxrun){
    fx <- first_derv(xt, x)
    fx_d <- second_derv(xt, x)
    if (fx == 0) {break}
    ht <- -fx/fx_d
    xt1 <- xt + ht
    if (abs(xt1-xt) < pre) {break}
    xt <- xt1
    n <- n+1
  }
  return(c(initial = init, root = xt, iter = n))
}

newton(theta_mom)
result2 <- as.data.frame(matrix(0,2,3))
init2 <- c(-2.7, 2.7)
for (i in 1:length(init2)) {
  result2[i,] <- newton(init2[i])
}
colnames(result2) <- c("Initial value", "Root", "Iteration #")
library(pander)
pander(result2, caption = "The result of Newton-Raphson method optimization")

## 200 starting values
init3 <- seq(-pi, pi, length.out = 200)
result3 <- as.data.frame(matrix(0,200,4))
for (i in 1:length(init3)) {
  result3[i,2:4] <- newton(init3[i])
  result3[i,1] <- i
}
colnames(result3) <- c("#", "Initial value", "Root", "Iteration #")
ggplot(result3, aes(result3[,2],result3[,3])) + geom_point(aes(col = "r")) + 
  labs(title = "Root of estimation with different initial values ", 
       x = "Initial values", y = expression(paste("Estimation of ", theta))) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none")
### Group
group_root <- result3
group_root[,3] <- round(group_root[,3],digits = 5)
library(gsubfn)
library(proto)
library(RSQLite)
library(sqldf)
group_root <- sqldf(
  'SELECT min([#]), [Initial value], Root, [Iteration #] 
  FROM group_root 
  GROUP BY Root'
  )
for (i in 1: dim(group_root)[1]) {
  if (i == dim(group_root)[1]) {
    group_root[i,1] <- paste(group_root[i,1], " - 200")
  } else {
    group_root[i,1] <- paste(group_root[i,1], " - ", as.numeric(group_root[i+1,1])-1)
    }
}
group_root <- group_root[,c(1,3)]
colnames(group_root) <- c("# Init", "Root")
pander(group_root, caption = "Groups with unique outcome of optimization")


# Modeling beetle data
## Gauss_Newton
beetles <- data.frame(
  days    = c(0,  8,  28,  41,  63,  69,   97, 117,  135,  154),
  beetles = c(2, 47, 192, 256, 768, 896, 1120, 896, 1184, 1024))
N0 <- beetles[which(beetles[,1]==0),2]
fomula1 <- beetles ~ K*N0/(N0+(K-N0)*exp(-r*days))
gauss_newton <- nls(fomula1, data = beetles, start = list(K=1184, r=0.5), trace = T)
gauss_newton
## Contour
### Squared Error
sqr_error <- function(r, K, sample = beetles) {
  sqr_error <- 0
  for (i in 1:dim(beetles)[1]) {
    sqr_error <- sqr_error + (K*N0/(N0+(K-N0)*exp(-r*sample[i,1]))-sample[i,2])^2
  }
  sqr_error
}
## Contour plot
K <- seq(10, 2000, 10)
r <- seq(0.01, 1, 0.01)
plot_data <- as.data.frame(matrix(0,(length(K)*length(r)),3))
colnames(plot_data) <- c("K", "r", "Squared Error")
for (i in 1:length(K)) {
  for (j in 1:length(r)) {
    plot_data[j+(i-1)*length(r),1] <- K[i]
    plot_data[j+(i-1)*length(r),2] <- r[j]
    plot_data[j+(i-1)*length(r),3] <- sqr_error(r[j], K[i])
  }
}
ggplot(plot_data, aes(x = K, y = r)) + 
  geom_contour(aes(z=plot_data[,3], col = ..level..), bins = 20) + 
  labs(title = "Contour Plot") + theme(plot.title = element_text(hjust = 0.5))

## BFGS
log_llh2 <- function (theta,   sample = beetles) {
  r <- theta[1]
  K <- theta[2]
  sigma_sqr <- theta[3]
  log_llh2 <- 0 
  for (i in 1: dim(sample)[1]) {
    log_llh2 <- - (log_llh2 + log(1/(sample[i,2]*(sigma_sqr*2*pi)^0.5)) -
      (log(sample[i,2])-log(K*N0/(N0+(K-N0)*exp(-r*sample[i,1]))))^2/(2*sigma_sqr))
  }
  log_llh2
}
result4 <- optim(c(0.2, 1200, 100), fn = log_llh2, method = "BFGS",hessian = T)
result4
var<- diag(solve(-result4$hessian))
var
