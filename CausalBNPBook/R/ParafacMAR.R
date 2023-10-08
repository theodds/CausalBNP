
#' Fit ignorable PARAFAC model to the observed data
#' 
#' Samples from the posterior of the ignorable PARAFAC model by Gibbs sampling
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
ParafacMAR <- function(Y, R, sigma_a, alpha = 1, K = NULL, state = NULL, nburn = 1000, nsave = 1000, nthin = 1) {

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
      a_beta = matrix(NA, nrow = nsave, ncol = J)
    )
  })

  GibbsIterate <- function(y) {
    y %>% UpdateClassMAR() %>% GetSufficientStatisticsMAR() %>% UpdateBeta() %>%
      UpdateOmega() %>% UpdateRhoBeta() %>% UpdateABeta() -> out
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
    output$omega[i,] <- state$omega
    output$alpha[i] <- state$alpha
    output$C[i,] <- state$C
    output$class_counts[i,] <- state$class_counts
    output$loglik_data[i,] <- state$loglik_data
    output$loglik[i] <- sum(state$loglik_data)
    output$rho_beta[i,] <- state$rho_beta
    output$a_beta[i,] <- state$a_beta

    if(i %% 100 == 0) cat("Finishing save iteration", i, "\n")

  }

  class(output) <- "parafacmar"
  
  return(output)

}

UpdateClassMAR <- function(state) {
  state %>% with({
    UpdateClassMARCpp(Y, R, log(beta), log(1 - beta), log(omega))
  }) -> out

  state$C <- out$C + 1
  state$loglik_data <- out$loglik_data

  return(state)
}

GetSufficientStatisticsMAR <- function(state) {
  state %>% with({
    UpdateSufficient(Y = Y, R = R, C = C, K = K)
  }) -> out
  state$success_counts_Y <- out$success_counts_Y
  state$failure_counts_Y <- out$failure_counts_Y
  state$class_counts <- out$class_counts

  return(state)
}
