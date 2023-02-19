FitERT <- function(X_train,
                   y_train,
                   X_test,
                   nburn              = 50,
                   nsave              = 100,
                   nthin              = 1,
                   nert               = 500,
                   alpha              = 0.95,
                   beta               = 2,
                   nu                 = 3,
                   q                  = 0.9,
                   update_lambda      = 0,
                   update_sigma_sq_mu = 1,
                   rescale            = "Rank",
                   alpha_lasso        = .1,
                   k = 2
                   )
{
  N_train <- dim(X_train)[1]
  P <- dim(X_train)[2]
  X <- rbind(X_train, X_test)
  
  L <- min(y_train)
  U <- max(y_train)
  yy_train <- (y_train - L) / (U - L) - 0.5
  
  ## Rescale X
  XX <- X
  N  <- dim(XX)[1]
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
  
  XXX <- list()
  for(i in 1:nert) {
    XXX[[i]] <- TreeCompress(alpha, beta, rep(1/P, P), XX)
    if(i %% 100 == 0) cat("Finishing random tree", i, "\n")
  }
  XXX <- do.call(cbind, XXX)
  XXX_train <- XXX[1:N_train, ]
  XXX_test <- XXX[-c(1:N_train), ]
  gc()

  fit <- cv.glmnet(XXX_train, y_train, alpha=alpha_lasso)
  fit$f_hat <- predict(fit, XXX_test)

  return(fit)
  
  ## K_train <- -log(XXX_train %*% t(XXX_train) / nert)
  ## K_test  <- -log(XXX_test %*% t(XXX_train) / nert)
  
  ## lambda <- 1
  ## sigma_sq_mu = .25^2
  
  ## tmp <- GetLambda(yy_train, X_train, nu, q)
  ## sigma_sq <- mean(tmp$sigma^2)
  ## shape_tau <- tmp$shape
  ## rate_tau  <- tmp$rate

  ## fit <- GPFit(
  ##          yy_train,
  ##          XX_train,
  ##          XX_test,
  ##          K_train,
  ##          K_test,
  ##          sigma_sq,
  ##          sigma_sq_mu,
  ##          lambda,
  ##          rate_tau,
  ##          shape_tau,
  ##          nburn,
  ##          nthin,
  ##          nsave,
  ##          update_sigma_sq_mu,
  ##          update_lambda
  ##        )
  
  ## fit$sigma_sq <- fit$sigma_sq * (U - L)^2
  ## ## fit$sigma_sq_mu <- fit$sigma_sq_mu * (U - L)^2
  ## fit$f_test <- (fit$f_test + 0.5) * (U - L) + L
  ## fit$f_test_mean <- rowMeans(fit$f_test)
  ## fit$v <- fit$v * (U - L)
  ## fit$midrange <- (U + L) / 2
  ## fit$sigma <- sqrt(fit$sigma_sq)
  ## fit$sigma_mu <- sqrt(fit$sigma_sq_mu)

  ## class(fit) <- "ERT"
  ## return(fit)
  
  ## library(MASS)
  ## Sigma_yy <- XX_train %*% t(XX_train)

  ## sigma_sq_mu <- (1 / (2 * k * sqrt(nert)))^2
  ## ## objective <- function(x) {
  ## ##   Sigma_hat <- x[1] / nert * Sigma_yy + x[2] * diag(length(yy_train))
  ## ##   out       <- -0.5 * N_train * log(2 * pi)
  ## ##   v         <- solve(Sigma_hat, yy_train)
  ## ##   out       <- out - 0.5 * determinant(Sigma_hat)$modulus
  ## ##   out       <- out - 0.5 * yy_train %*% v

  ## ##   if(x[2] < 1E-16 | x[1] < 1E-16) out <- -Inf

  ## ##   return(-out)
  ## ## }

  ## objective <- function(x) {
  ##   Sigma_hat <- sigma_sq_mu * Sigma_yy + x[2] * diag(length(yy_train))
  ##   out       <- -0.5 * N_train * log(2 * pi)
  ##   v         <- solve(Sigma_hat, yy_train)
  ##   out       <- out - 0.5 * determinant(Sigma_hat)$modulus
  ##   out       <- out - 0.5 * yy_train %*% v

  ##   if(x[2] < 1E-16 | x[1] < 1E-16) out <- -Inf

  ##   return(-out + abs(x[1]))
  ## }
  
  ## empirical_bayes <- optim(c(1, 1), objective)
  ## ## sigma_sq_mu     <- empirical_bayes$par[1]
  ## sigma_sq        <- empirical_bayes$par[2]
  ## Sigma_hat       <- sigma_sq_mu * Sigma_yy + sigma_sq * diag(length(yy_train))
  ## alpha_hat       <- sigma_sq_mu * t(XX_train) %*% solve(Sigma_hat, yy_train)
  ## ## Sigma_hat       <- sigma_sq_mu / nert * Sigma_yy + sigma_sq * diag(length(yy_train))
  ## ## alpha_hat       <- sigma_sq_mu / nert * t(XX_train) %*% solve(Sigma_hat, yy_train)
  ## f_hat           <- XX_test %*% alpha_hat
  
  ## fit <- list()
  ## fit$empirical_bayes <- empirical_bayes
  ## fit$alpha_hat <- alpha_hat
  ## fit$f_hat     <- f_hat
  ## fit$sigma_sq  <- sigma_sq * (U - L)^2

  ## fit$f_hat <- (fit$f_hat + 0.5) * (U - L) + L

  return(fit)
}
                   
