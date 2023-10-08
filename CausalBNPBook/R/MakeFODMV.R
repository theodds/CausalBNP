#' Complete Case Missing Value G-computatin
#' 
#' Conducts G-computation under the observed data missing value assumption
#' 
#' @param chain Object of type parafac, which contains in particular the Markov chain with samples of the posterior
#' @param n_sim The number of samples to take during the G-computation; defaults to 100 times the MCMC length
#' @param j_0 The timepoint at which to compute the expected value, by default this is the last observation time
#' @param xi_f A function which samples from the prior distribution of the tilting parameters.
#' @param iterations the iterations of the chain to look at, by default this is 1:I where I is the number of samples collected
#' 
#' @return Returns a data frame contaning
#' \itemize{
#' \item iteration: the iteration the sample is taken
#' \item mean: the marginal mean of Y_{j_0}
#' \item std_error: the standard error of the estimate of Y_{j_0}
#' }

GcompODMV <- function(chain, n_sim = NULL, j_0 = NULL, xi_f = function() 0, iterations = NULL) {
  if(is.null(iterations)) iterations <- 1:nrow(chain$beta)
  J <- chain$J
  if(is.null(j_0)) j_0 <- J
  
  f <- MakeFODMV(chain, n_sim, j_0, J, xi_f)
  results <- lapply(iterations, f)
  out <- tibble(iteration = iterations)
  out$mean <- sapply(results, function(x) x$mu_g)
  out$std_error <- sapply(results, function(x) x$sd_g)
  
  return(out)
}

MakeFODMV <- function(chain, n_sim = NULL, j_0 = 8, J = 8, xi_f = function() 0) {

  if(is.null(n_sim)) {
    n_sim <- length(chain$loglik) * 100
  }

  f <- function(i) {
    g_comp <- with(chain, ParafacGcomp( N_sim = n_sim,
                                     J = J,
                                     j_0 = j_0,
                                     omega = omega[i,],
                                     log_omega = log(omega[i,]),
                                     gamma = gamma[i,,],
                                     log_gamma =log(gamma[i,,]) ,
                                     log_1_m_gamma = log(1 - gamma[i,,]),
                                     beta = beta[i,,],
                                     log_beta = log(beta[i,,]),
                                     log_1_m_beta = log(1 - beta[i,,]),
                                     xi = xi_f()))

    mu_g <- mean(g_comp)
    sd_g <- sd(g_comp) / sqrt(length(g_comp))

    if(i %% 100 == 0) cat("Finishing gcomp iter", i, "\n")

    return(list(mu_g = mu_g, sd_g = sd_g))
  }
  f
}
