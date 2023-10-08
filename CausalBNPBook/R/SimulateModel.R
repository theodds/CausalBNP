MakeFakeData <- function(K = 10, J = 6, N = 200, alpha = 3) {

  ## Simulate weights
  omega <- rdirichlet(1, rep(alpha / K, K))

  ## Simulate PARAFAC missingness parameters
  theta <- matrix(rbeta(J * K, 1.4, .6), K, J)

  ## Draw hyperparameters a and rho
  a <- array(data = 1 / rgamma(n = 2 * 2 * J, shape = 1, rate = 1),
             dim = c(2, 2, J))
  rho <- array(data = runif(n = 2 * 2 * J), dim = c(2,2,J))

  ## Simulate mixture component parameters
  beta <- array(0, dim = c(2,2,J,K))
  for(r in 1:2) {
    for(y in 1:2) {
      for(j in 1:J) {
        beta[r,y,j,] <- rbeta(n = K,
                              shape1 = a[r,y,j] * rho[r,y,j],
                              shape2 = a[r,y,j] * (1-rho[r,y,j]))
      }
    }
  }

  ## Simulate latent classes
  latent_class <- sample(x = 1:K, size = N, replace = TRUE, prob = omega)

  ## Simulate data

  prob_matrix <- matrix(runif(N * J), N, J)
  R <- ifelse(prob_matrix < theta[latent_class, ], 1, 0)

  Y <- matrix(data = 0, nrow = N, ncol = J)
  for(i in 1:N) {
    k <- latent_class[i]
    r <- 1
    y <- 1
    Y[i,1] <- rbinom(n = 1, size = 1, prob = beta[r, y, 1, k])
    for(j in 2:J) {
      r <- R[i,j-1] + 1
      y <- Y[i,j-1] + 1
      Y[i,j] <- rbinom(n = 1, size = 1, prob = beta[r,y,j,k])
    }
  }



  Yobs <- ifelse(R == 1, Y, NA)

  return(list(Y = Y, Yobs = Yobs, R = R, latent_class = latent_class,
              beta = beta, a = a, rho = rho, theta = theta, omega = omega))

}
