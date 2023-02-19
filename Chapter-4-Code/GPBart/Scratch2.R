x <- as.factor(c("a","b","a","c"))
x

X_1 <- matrix(runif(100 * 10), nrow=100)
X_2 <- matrix(sample(1:10, size = 100*10, replace = TRUE), nrow=100)

out <- GetRankKernCat(X_1, X_2)

as.numeric(x)


data_dir <- c("../BartGP/BartGP-R/data")
dat_files <- paste0(data_dir, "/", list.files(data_dir))

baseball <- read.csv(dat_files[32])
head(baseball)

cat_idx <- c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,18,19)
qual_idx <- setdiff(1:dim(baseball)[2], cat_idx)

N <- dim(baseball)[1]
N_train <- floor(2 * N / 3)
N_test <- N - N_train
P_cols <- 1:(dim(baseball)[2]-1)
P_cols <- 2:(dim(baseball)[2]-1)
y_col <- dim(baseball)[2]

set.seed(1235)
idx_train <- sample(1:N, N_train)
idx_test <- setdiff(1:N, idx_train)

CrossVal <- function() {
    idx_train <- sample(1:N, N_train)
    idx_test <- setdiff(1:N, idx_train)
    
    X_train <- as.matrix(baseball[idx_train, P_cols])
    X_test <- as.matrix(baseball[idx_test, P_cols])
    y_train <- baseball[idx_train, y_col]
    y_test <- baseball[idx_test, y_col]
    
    foo <- bart(x.train = X_train, y.train = as.double(y_train), x.test = X_test)
  #   foo2 <- bgp(X = X_train, Z = as.double(y_train), XX = X_test)
    foo_bart <- BartGP(X_train = X_train, y_train = y_train, X_test = X_test)
    
    rmse_bart <- sqrt(mean((foo$yhat.test.mean - y_test)^2))
  #   rmse_bgp <- sqrt(mean((foo2$ZZ.mean - y_test)^2))
    rmse_bartgp <- sqrt(mean((foo_bart - y_test)^2))
    
    return(c(rmse_bartgp=rmse_bartgp, rmse_bart=rmse_bart))
  #   return(c(rmse_bart=rmse_bart, rmse_bgp=rmse_bgp, rmse_bartgp=rmse_bartgp))
    
  }

CrossVal2 <- function() {
  idx_train <- sample(1:N, N_train)
  idx_test <- setdiff(1:N, idx_train)
  
  X_train <- as.matrix(baseball[idx_train, P_cols])
  X_test <- as.matrix(baseball[idx_test, P_cols])
  y_train <- baseball[idx_train, y_col]
  y_test <- baseball[idx_test, y_col]
  
#   foo <- bart(x.train = X_train, y.train = as.double(y_train), x.test = X_test)
  #   foo2 <- bgp(X = X_train, Z = as.double(y_train), XX = X_test)
  foo_bart <- BartGP(X_train = X_train, y_train = y_train, X_test = X_test)
  
#   rmse_bart <- sqrt(mean((foo$yhat.test.mean - y_test)^2))
  #   rmse_bgp <- sqrt(mean((foo2$ZZ.mean - y_test)^2))
  rmse_bartgp <- sqrt(mean((foo_bart - y_test)^2))
  
  return(c(rmse_bartgp=rmse_bartgp))
  #   return(c(rmse_bart=rmse_bart, rmse_bgp=rmse_bgp, rmse_bartgp=rmse_bartgp))
  
}

X_train <- as.matrix(baseball[idx_train, P_cols])
X_test <- as.matrix(baseball[idx_test, P_cols])
y_train <- baseball[idx_train, y_col]
y_test <- baseball[idx_test, y_col]

foo <- bart(x.train = X_train, y.train = as.double(y_train), x.test = X_test)
foo2 <- bgpllm(X = X_train, Z = as.double(y_train), XX = X_test)
foo_bart <- BartGP(X_train = X_train, y_train = y_train, X_test = X_test)

sqrt(mean((foo$yhat.test.mean - y_test)^2))
sqrt(mean((foo2$ZZ.mean - y_test)^2))
sqrt(mean((foo_bart - y_test)^2))

BartGP <- function(X_train, y_train, X_test) {
  
  y_stand <- Standardize(y = y_train)
  
  U <- y_stand$U
  L <- y_stand$L
  y_stand <- y_stand$y
  
  X <- rbind(X_train, X_test)
  R <- matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
  
  idx_train <- 1:dim(X_train)[1]
  idx_test <- (dim(X_train)[1]+1):(dim(X_train)[1]+dim(X_test)[1])
  
  P <- dim(X)[2]
  for(p in 1:P) {
    R[,p] <- rank(X[,p])
  }

  tau <- 1/(2 * 2)
  Ed <- 3
#   p <- 1 / Ed
  Omega <- GetRankKern(R)
  N <- dim(X)[1]
  N_train <- dim(X_train)[1]
  #   M <- exp(-as.matrix(dist(X)^2))
#   M <- p * M / (1 - (1 - p) * M) + diag(N) * .00001
  M <- exp(Ed * (Omega - 1))
  Sigma <- tau * tau * M
  
  
  Sigma_train <- Sigma[idx_train, idx_train]
  Sigma_test <- Sigma[idx_test, idx_train]
  
  nu <- 3
  q <- 0.9
  lambda <- GetLambda(y_stand, X_train, nu, q)$out
  
  loglik <- function(sigma) {  
    K <- sigma^2 * diag(N_train) + Sigma_train
    
    
    
    out <- dmvnorm(x=y_stand, mean=numeric(N_train), sigma=K, log=TRUE)
    
    out <- out + dgamma(1.0 / sigma^2, shape=nu/2, rate=nu*lambda/2, log=TRUE) - 3 * log(sigma)
    # out <- out + dgamma(sigma, shape = 10, scale = 1 / (U - L) / 10, log=TRUE)
    
    return(out)
  }
  
  loglik_2 <- function(sigma, tau, lambda) {
    K <- tau * tau * exp(lambda * (Omega[idx_train, idx_train] - 1)) + sigma * sigma * diag(N_train)
    out <- dmvnorm(x=y_stand, mean=numeric(N_train), sigma=K, log=TRUE)
    
    # out <- out + dgamma(1.0 / sigma^2, shape=nu/2, rate=nu*lambda/2, log=TRUE)
    
    return(out)
  }
  
  sigma_est <- optimize(f = loglik, interval = c(.001, 1), maximum = TRUE)$max * 2
  sigma_est_2 <- optim(par = c(8,tau=8,lambda=8), fn = function(x) -loglik_2(x[1],x[2],x[3]), method = "L-BFGS-B", lower=c(0.01,0.01,0.01))
  
  sigma_est <- sigma_est_2$par[1]
  tau_est <- sigma_est_2$par[2]
  lambda_est <- sigma_est_2$par[3]
  Sigma_est <- tau_est* tau_est * exp(lambda_est * Omega - 1)
  Sigma_train <- Sigma_est[idx_train, idx_train]
  Sigma_test <- Sigma_est[idx_test, idx_train]
  
  K <- (sigma_est)^2 * diag(N_train) + Sigma_train
  K_inv <- solve((sigma_est)^2 * diag(N_train) + Sigma_train)
  
  y_pred <- Sigma_test %*% K_inv %*% y_stand
  
  y_pred <- y_pred + 0.5
  y_pred <- y_pred * (U - L)
  y_pred <- y_pred + L
  
  return(list(y_pred=y_pred, sigma_est=sigma_est_2))
  
}



Setup <- function() {
  N_train <<- 100
  N_test <<- 1000
  N <<- N_train + N_test
  idx_train <<- 1:N_train
  idx_test <<- (N_train+1):N
  P <<- 10
  X <<- matrix(runif(N * P), nrow = N)
  X_train <<- X[idx_train, ]
  X_test <<- X[idx_test, ]
  y <<- Fried(X) + rnorm(N)
  y_train <<- y[idx_train]
  y_test <<- y[idx_test]
}

Setup()
y_pred_1  <- BartGP(X_train, y_train, X_test)
sigma_est <- y_pred_1$sigma_est
y_pred_1  <- y_pred_1$y_pred
bo_1      <- bart(x.train = X_train, y.train = as.double(y_train), x.test = X_test)

rmse_1 <- sqrt(mean((y_pred_1 - y_test)^2))
rmse_2 <- sqrt(mean((bo_1$yhat.test.mean - y_test)^2))
print(rmse_1 / rmse_2)
