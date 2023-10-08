#' Initialize the Parafac model
#'
#' Initializes the nonignorable Parafac model by sampling from the prior.
#'
#' @param N the number of observations
#' @param J the number of time points
#' @param K the number of latent classes to truncate to
#' @param alpha Concentration parameter for the Dirichlet process approximation
#'
#' @return An object of type parafacstate, holding the state of the Markov chain.

InitParafac <- function(N, J, K, alpha = 1) {

  beta <- matrix(0, nrow = K, ncol = J)
  gamma <- matrix(0, nrow = K, ncol = J)

  omega <- rdirichlet(n = 1, alpha = rep(alpha / K, K))
  rho_beta <- runif(J)
  rho_gamma <- runif(J)
  Z_beta <- runif(J)
  Z_gamma <- runif(J)
  a_beta <- K * Z_beta / (1 - Z_beta)
  a_gamma <- K * Z_gamma / (1 - Z_gamma)

  for(j in 1:J) {
    beta[,j] <- rbeta(K, rho_beta[j] * a_beta[j], (1 - rho_beta[j]) * a_beta[j])
    gamma[,j] <- rbeta(K, rho_gamma[j] * a_gamma[j], (1 - rho_gamma[j]) * a_gamma[j])
  }

  C <- sample(x = 1:K, size = N, replace = TRUE, prob = omega)

  Y <- matrix(0, nrow = N, ncol = J)
  R <- matrix(0, nrow = N, ncol = J)
  for(n in 1:N) {
    Y[n,] <- rbinom(n = J, size = 1, prob = beta[C[n], ])
    R[n,] <- rbinom(n = J, size = 1, prob = gamma[C[n], ])
  }

  Y <- ifelse(R == 1, Y, NA)

  success_counts_Y <- matrix(0, nrow = K, ncol = J)
  failure_counts_Y <- matrix(0, nrow = K, ncol = J)
  success_counts_R <- matrix(0, nrow = K, ncol = J)
  failure_counts_R <- matrix(0, nrow = K, ncol = J)

  for(k in 1:K) {
    Y_counts <- Y[C == k, , drop = FALSE]
    R_counts <- R[C == k, , drop = FALSE]
    success_counts_Y[k,] <- colSums(Y_counts, na.rm = TRUE)
    failure_counts_Y[k,] <- colSums(1 - Y_counts, na.rm = TRUE)
    success_counts_R[k,] <- colSums(R_counts, na.rm = TRUE)
    failure_counts_R[k,] <- colSums(1 - R_counts, na.rm = TRUE)
  }


  state <- list(Y = Y,
                R = R,
                beta = beta,
                gamma = gamma,
                omega = omega,
                alpha = alpha,
                C = C,
                K = K,
                J = J,
                N = N,
                rho_beta = rho_beta,
                rho_gamma = rho_gamma,
                a_beta = a_beta,
                a_gamma = a_gamma,
                sigma_a = K,
                success_counts_Y = success_counts_Y,
                failure_counts_Y = failure_counts_Y,
                success_counts_R = success_counts_R,
                failure_counts_R = failure_counts_R,
                class_counts = tabulate(C, K)
  )

  class(state) <- "parafacstate"

  return(state)

}

UpdateBeta <- function(state) {
  state %>% with({
    UpdateBetaCpp(success_counts = success_counts_Y,
                  failure_counts = failure_counts_Y,
                  col_shape_1 = rho_beta * a_beta,
                  col_shape_2 = (1 - rho_beta) * a_beta
    )
  }) -> beta_up
  state$beta <- beta_up
  return(state)
}

UpdateGamma <- function(state) {
  state %>% with({
    UpdateBetaCpp(success_counts = success_counts_R,
                  failure_counts = failure_counts_R,
                  col_shape_1 = rho_gamma * a_gamma,
                  col_shape_2 = (1 - rho_gamma) * a_gamma)
  }) -> gamma_up
  state$gamma <- gamma_up
  return(state)
}

UpdateOmega <- function(state) {
  state %>% with({
    rdirichlet(n = 1, alpha = alpha / K + class_counts)
  }) -> omega_up
  state$omega <- omega_up
  return(state)
}

GetSufficientStatistics <- function(state) {
  state %>% with({
    UpdateSufficient(Y = Y, R = R, C = C, K = K)
  }) -> out
  state$success_counts_Y <- out$success_counts_Y
  state$success_counts_R <- out$success_counts_R
  state$failure_counts_Y <- out$failure_counts_Y
  state$failure_counts_R <- out$failure_counts_R
  state$class_counts <- out$class_counts

  return(state)
}

UpdateRhoGamma <- function(state, w = .04) {
  rho_gamma_up <- numeric(state$J)
  for(j in 1:state$J) {

    x <- state$gamma[state$class_counts > 0, j]

    log_posterior <- function(y) {
      sum(dbeta(x = x,
                shape1 = y * state$a_gamma[j],
                shape2 = (1 - y) * state$a_gamma[j], log = TRUE))
    }
    rho_gamma_up[j] <- uni.slice(x0 = state$rho_gamma[j],
                                 g = log_posterior,
                                 w = w, lower = 0, upper = 1)
  }
  state$rho_gamma <- rho_gamma_up
  return(state)
}

UpdateRhoBeta <- function(state, w = .04) {
  rho_beta_up <- numeric(state$J)
  for(j in 1:state$J) {

    x <- state$beta[state$class_counts > 0, j]

    log_posterior <- function(y) {

      sum(dbeta(x = x,
                shape1 = y * state$a_beta[j],
                shape2 = (1 - y) * state$a_beta[j], log = TRUE))
    }
    rho_beta_up[j] <- uni.slice(x0 = state$rho_beta[j],
                                g = log_posterior,
                                w = w, lower = 0, upper = 1)
  }
  state$rho_beta <- rho_beta_up
  return(state)
}

UpdateABeta <- function(state, w = .1) {
  a_beta_up <- numeric(state$J)
  for(j in 1:state$J) {

    x <- state$beta[state$class_counts > 0, j]

    log_posterior <- function(y) {
      a <- state$sigma_a * y / (1 - y)
      shape1 <- a * state$rho_beta[j]
      shape2 <- a * (1 - state$rho_beta[j])
      out <- sum(dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE))
      return(out)

    }
    z0 <- state$a_beta[j] / (state$a_beta[j] + state$sigma_a)
    z <- uni.slice(x0 = z0, g = log_posterior, w = w, lower = 0, upper = 1)
    a_beta_up[j] <- state$sigma_a * z / (1 - z)
  }

  state$a_beta <- a_beta_up
  return(state)
}

UpdateAGamma <- function(state) {
  a_gamma_up <- numeric(state$J)
  for(j in 1:state$J) {

    x <- state$gamma[state$class_counts > 0, j]

    log_posterior <- function(y) {
      a <- state$sigma_a * y / (1 - y)
      shape1 <- state$rho_gamma[j] * a
      shape2 <- (1 - state$rho_gamma[j]) * a
      out <- sum(dbeta(x, shape1 = shape1, shape2 = shape2, log = TRUE))

      return(out + log(state$sigma_a) - 2 * log(state$sigma_a + y))
    }

    z0 <- state$a_gamma[j] / (state$a_gamma[j] + state$sigma_a)
    z <- uni.slice(x0 = z0, g = log_posterior, w = state$sigma_a, lower = 0, upper = 1)
    a_gamma_up[j] <- state$sigma_a * z / (1 - z)
  }

  state$a_gamma <- a_gamma_up
  return(state)
}

UpdateClass <- function(state) {
  state %>% with({
    UpdateClassCpp(Y, R, log(beta), log(1 - beta), log(gamma), log(1 - gamma), log(omega))
  }) -> out

  state$C <- out$C + 1
  state$loglik_data <- out$loglik_data

  return(state)
}

#' Fit nonignorable PARAFAC model to the observed data
#'
#' Samples from the posterior of the nonignorable PARAFAC model by Gibbs sampling
#'
#' @param Y a binary N x J matrix with missing values denoted by NA. The observed data.
#' @param R a binary N x J matrix of missingness indicators.
#' @param sigma_a The shrinkage parameter for uniform shrinkage prior;
#'        very roughly, can be set to the expected number of nonempty classes,
#'        which a priori is roughly alpha * (digamma(alpha + N) - digamma(alpha)) + 1
#' @param alpha The dirichlet process concentration parameter, = 1 by default
#' @param K the number of latent classes in the truncated Dirichlet approximation. Defaults to sqrt(nrow(Y))
#' @param state An existing object of type parafacstate. If provided, (Y,R,K) are ignored and samples are based on this state.
#' @param nburn The number of warmup iterations for the algorithm
#' @param nsave The number of samples to be collected
#' @param nthin The thinning interval
#'
#' @return Returns an object of type parafac, consisting of the parameters samples from the posterior.

ParafacMNAR <- function(Y, R, sigma_a, alpha = 1, K = NULL, state = NULL, nburn = 1000, nsave = 1000, nthin = 1) {

  J <- N <- NULL

  if(is.null(K)) {
    K <- floor(sqrt(nrow(Y)))
  }

  if(is.null(state)) {
    state <- InitParafac(nrow(Y), ncol(Y), K)
    state$Y <- Y
    state$R <- R
    state$sigma_a <- sigma_a
    state$alpha <- 1
  }

  output <- with(state, {
    list(beta = array(NA, c(nsave, K, J)),
         gamma = array(NA, c(nsave, K, J)),
         omega = matrix(NA, nrow = nsave, ncol = K),
         alpha = numeric(nsave),
         C = matrix(NA, nrow = nsave, ncol = N),
         class_counts = matrix(NA, nrow = nsave, ncol = K),
         loglik_data = matrix(NA, nrow = nsave, ncol = N),
         loglik = numeric(nsave),
         K = K,
         J = J,
         N = N,
         rho_beta = matrix(NA, nrow = nsave, ncol = J),
         rho_gamma = matrix(NA, nrow = nsave, ncol = J),
         a_beta = matrix(NA, nrow = nsave, ncol = J),
         a_gamma = matrix(NA, nrow = nsave, ncol = J))
  })

  GibbsIterate <- function(y) {
    y %>% UpdateClass() %>% GetSufficientStatistics() %>% UpdateBeta() %>%
      UpdateGamma() %>% UpdateOmega() %>% UpdateRhoGamma() %>% UpdateRhoBeta %>%
      UpdateAGamma() %>% UpdateABeta() -> out
    return(out)
  }

  for(i in 1:nburn) {
    state <- GibbsIterate(state)

    if(i %% 100 == 0) cat("Finishing burnin iteration", i, "\n")
  }

  for(i in 1:nsave) {
    for(b in 1:nthin) {
      state <- GibbsIterate(state)
    }

    output$beta[i,,] <- state$beta
    output$gamma[i,,] <- state$gamma
    output$omega[i,] <- state$omega
    output$alpha[i] <- state$alpha
    output$C[i,] <- state$C
    output$class_counts[i,] <- state$class_counts
    output$loglik_data[i,] <- state$loglik_data
    output$loglik[i] <- sum(state$loglik_data)
    output$rho_beta[i,] <- state$rho_beta
    output$rho_gamma[i,] <- state$rho_gamma
    output$a_beta[i,] <- state$a_beta
    output$a_gamma[i,] <- state$a_gamma

    if(i %% 100 == 0) cat("Finishing save iteration", i, "\n")

  }

  class(output) <- "parafac"

  return(output)

}
