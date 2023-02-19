library(BayesTree)
library(tgp)
library(GPFDA)
library(TonyColors)
library(GPBart)
library(mvtnorm)

source("Funcs.R")

set.seed(1237)
N <- 100
N_test <- 1000
P <- 10
X <- FriedSamp(n = N, dim_x = P)
y_train <- X$Y
X_train <- X$X
X <- FriedSamp(n = N_test, dim_x = P)
y_test <- X$Y
X_test <- X$X

library(codetools)
fit_gp <- FitGP(X_train = X_train, y_train = y_train, X_test = X_test, rescale = "None", nburn = 500, nsave = 500)
f_hat <- fit_gp$f_test_mean

library(tgp)
fit_tgp <- bgp(X = as.data.frame(X_train), Z = y_train, XX = as.data.frame(X_test), corr = "exp")
f_hat_3 <- fit_tgp$ZZ.mean

library(bartMachine)
fit <- bartMachine(X = as.data.frame(X_train), y = y_train, num_burn_in = 1000, num_iterations_after_burn_in = 4000)
f_hat_4 <- predict(fit, as.data.frame(X_test))

RMSE(f_hat, Fried(X_test))
RMSE(f_hat_3, Fried(X_test))
RMSE(f_hat_4, Fried(X_test))


foo2 <- bgpllm(X = as.data.frame(X_train), Z = y_train, XX=as.data.frame(X_test))
foo3 <- btgp(X = as.data.frame(X_train), Z = y_train, XX=as.data.frame(X_test))
foo4 <- btgpllm(X = as.data.frame(X_train), Z = y_train, XX=as.data.frame(X_test))

plot(y_test, foo$ZZ.mean)
abline(a=0,b=1)

e <- (foo4$ZZ.mean + 0.5) * (U-L) + L - Fried(X_test)
sqrt(mean(e * e))

ham <- predict(foo, as.data.frame(X_test))


for(i in 1:100) {
  
  N_train <- 100
  N_test <- 1000
  N <- N_train + N_test
  idx_train <- 1:N_train
  idx_test <- (N_train+1):N
  P <- 10
  X <- matrix(runif(N * P), nrow = N)
  X_train <- X[idx_train, ]
  X_test <- X[idx_test, ]

  
  y <- Fried(X) + rnorm(N)
  y_train <- y[idx_train]
  y_test <- y[idx_test]
  
  R <- matrix(0, nrow = N, ncol=P)
  

  
  bo_2 <- bart(x.train = X_train, y.train = y_train, x.test = X_test, ndpost = 1000, nskip = 200)
  hist_sigma <- hist(bo_2$sigma)
  sigma_factor <- max(hist_sigma$density)
  sigma_min <- min(bo_2$sigma)
  sigma_max <- max(bo_2$sigma)
  
  

  
  y <- Standardize(y)
  U <- y$U
  L <- y$L
  y <- y$y
  sigma_true <- 1 / (U - L)
  y_train <- y[idx_train]
  y_test  <- y[idx_test]
  
  for(p in 1:P) {
    R[,p] <- rank(X[,p])
  }
  
  tau <- 1/4
  Ed <- 2
  p <- 1 / Ed
  # M <- GetRankKern(R)
  M <- exp(-as.matrix(dist(X)^2))
  Sigma <- tau * tau * (M + 0.00001 * diag(N))
  # Sigma <- tau * tau * p * M / (1 - M * (1-p))
  
  Sigma_train <- Sigma[idx_train, idx_train]
  Sigma_test <- Sigma[idx_test, idx_train]
  
  nu <- 3
  q <- 0.9
  lambda <- GetLambda(y_train, X_train, nu, q)$out
  
  loglik <- function(sigma) {  
    K <- sigma^2 * diag(N_train) + Sigma_train
    
    
    
    out <- dmvnorm(x=y_train, mean=numeric(N_train), sigma=K, log=TRUE)
    
    out <- out + dgamma(1.0 / sigma^2, shape=nu/2, rate=nu*lambda/2, log=TRUE) - 3 * log(sigma)
    # out <- out + dgamma(sigma, shape = 10, scale = 1 / (U - L) / 10, log=TRUE)
    
    return(out)
  }
  
  posterior <- function(sigma) {
    f <- Vectorize(loglik)
    out <- f(sigma)
    return(exp(out - max(out)))
  }
  
  hist(bo_2$sigma, xlim=c(sigma_min/3, 3*sigma_max), freq=FALSE)
  plot(function(x) posterior(x / (U-L)) * sigma_factor, from=sigma_min/3, to=3*sigma_max, add=TRUE, col='red')
  plot(function(x) posterior(x / (U-L)) * 7, from=.1, to=.6, add=TRUE, col='red')
  
  sigma_est <- (U - L) * optimize(f = loglik, interval = c(.001, 1), maximum = TRUE)$max
  
  
  K <- (sigma_est)^2 * diag(N_train) + Sigma_train
  K_inv <- solve((sigma_est)^2 * diag(N_train) + Sigma_train)
  
  y_pred <- Sigma_test %*% K_inv %*% y_train
  
  RMSE <- function(y, X) {
    y_norm <- y + 0.5
    y_norm <- y_norm * (U - L) + L
    e <- y_norm - Fried(X)
    sqrt(mean(e * e))
  }
  
  foo <- c(foo, RMSE(y_pred, X_test))
  sigest <- c(sigest, sigma_est)
  
  bo <- bart(x.train = X_train, y.train = y_train, x.test = X_test,nskip = 5000,ndpost = 5000)
  
}

CI <- function(x) {
  mean(x) + c(-1, 1) * 2 * sd(x) / sqrt(length(x))
}

CI(foo)
mean(foo)

plot(y_test, y_pred)




posterior <- function(sigma) {
  out <- Vectorize(loglik)(sigma)
  out <- out - max(out)
  return(exp(out))
}



plot(loglik, from=.01, to=.1, n=50)
abline(v=sigma_true)

alpha <- Sigma_inv %*% y
loglik <- sum(alpha * alpha) - sum(diag(solve(K_y)))

rank(X[,1])
