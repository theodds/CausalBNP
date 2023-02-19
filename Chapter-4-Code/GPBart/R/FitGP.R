FitGP <- function(X_train, y_train, X_test,
                  nburn              = 50,
                  nsave              = 100,
                  nthin              = 1,
                  cov                = "OU",
                  rescale            = "Rank",
                  nu                 = 3,
                  q                  = 0.9,
                  update_lambda      = 1,
                  update_sigma_sq_mu = 1)
{

  ## Constants
  N_train <- length(y_train)
  N_test  <- dim(X_test)[1]
  N       <- N_train + N_test
  P       <- dim(X_test)[2]

  ## Preprocess Y
  L <- min(c(y_train))
  U <- max(c(y_train))
  yy_train <- (y_train - L) / (U - L) - .5
  ## yy_test  <- (y_test - L) / (U - L) - .5
  

  ## Initial values for things
  lambda <- 1.2
  sigma_sq_mu <- .25^2

  ## Empirical Prior for sigma
  tmp <- GetLambda(yy_train, X_train, nu, q)
  sigma_sq <- mean(tmp$sigma^2)
  shape_tau <- tmp$shape
  rate_tau  <- tmp$rate

  X <- rbind(X_train, X_test)
  
  ## Rescale X
  XX <- X
  stopifnot(rescale %in% c("Rank", "Uniform", "None"))
  if(rescale == "Rank") {
    for(p in 1:P) {
      XX[,p] <- rank(X[,p]) / N
    }
  }
  if(rescale == "Uniform") {
    for(p in 1:P) {
      u <- max(XX[,p])
      ell <- min(XX[,p])
      XX[,p] <- (XX[,p] - ell) / (u - ell)
    }
  }
  XX_train <- XX[1:N_train, ]
  XX_test <- XX[-c(1:N_train), ]

  ## Calculate K_train and K_test
  if(cov == "OU") {
    K <- as.matrix(dist(XX, method = "manhattan")) / P
  } else {
    K <- as.matrix(dist(XX, method = "euclidean")) / P
  }
  K_train <- K[1:N_train, 1:N_train]
  K_test  <- K[-c(1:N_train), 1:N_train]

  fit <- GPFit(yy_train, XX_train, XX_test, K_train, K_test, sigma_sq, sigma_sq_mu,
               lambda, rate_tau, shape_tau, nburn, nthin, nsave,
               update_sigma_sq_mu, update_lambda)

  fit$sigma_sq <- fit$sigma_sq * (U - L)^2
  fit$sigma_sq_mu <- fit$sigma_sq_mu * (U - L)^2
  fit$f_test <- (fit$f_test + 0.5) * (U - L) + L
  fit$f_test_mean <- rowMeans(fit$f_test)
  fit$v <- fit$v * (U - L)
  fit$midrange <- (U + L) / 2
  fit$sigma <- sqrt(fit$sigma_sq)
  fit$sigma_mu <- sqrt(fit$sigma_sq_mu)

  class(fit) <- "GPFit"

  return(fit)
}

GetLambda <- function(y, X, nu=3, q=.9) {
  fit      <- lm(y ~ X)
  sigma_sq <- summary(fit)$sigma^2
  sigma    <- sqrt(sigma_sq)
  
  out <- qgamma(1 - q, shape=nu/2, rate=nu/2)
  out <- out * sigma_sq
  
  test_it <- sqrt(1.0 / rgamma(10000, nu/2, rate=out*nu/2))
  conf <- sum(test_it < sigma)
  
  out_list <- list(out=out, test_it=test_it, sigma=sigma, shape = nu / 2, 
                   rate = out * nu / 2, conf=conf)
  
  return(out_list)
}
