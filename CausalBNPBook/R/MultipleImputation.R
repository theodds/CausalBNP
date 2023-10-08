#' Perform multiple imputations using a Dirichlet process mixture of multinomial models
#'
#' This function performs multiple imputations on a dataset using a specified
#' method and a parafac object obtained from the ParafacMNAR (or ParafacMAR)
#' model fitting functions. The imputations are performed for a specified number
#' of times, and the results are returned in a list of data frames.
#'
#' @param Y A binary N x J matrix with missing values denoted by NA. The observed data.
#' @param R A binary N x J matrix of missingness indicators.
#' @param chain A list containing the parameters sampled from the posterior in a DPM of multinomial models (an object of type "parafac" or "parafacmar").
#' @param num_impute The number of imputations to perform. Defaults to 20.
#' @param method The method to use for imputation. One of "CCMV", "TLO", "NIP", or "MAR". Defaults to "CCMV". See monograph for details.
#' @param j_0 The time to impute data, required for "NIP" and "TLO" methods. If NULL, an error is thrown when using "NIP" or "TLO" methods.
#' @param xi A numeric vector of sensitivity parameters of length num_impute, used in the "TLO" (tilted last occasion) method. Defaults to a vector of zeros.
#'
#' @return A list of length num_impute, where each element is the result of an imputation.
#'
#' @seealso
#' \code{\link{ParafacMNAR}} and \code{\link{ParafacMAR}} for fitting the PARAFAC model to obtain the `chain` parameter.
#'
#' @details
#' The function performs multiple imputations based on the specified method and
#' the parameters obtained from a PARAFAC model fitting procedure. The `chain`
#' parameter should be an object of type "parafac" obtained from a function like
#' `ParafacMNAR`, or for imputing under MAR should be of type "parafacmar".
#'
#' The function supports four methods for imputation: "CCMV", "TLO", "NIP", and
#' "MAR". Depending on the method chosen, different parameters from the `chain`
#' object are used in the imputation procedure.
ParafacMI <- function(Y, R, chain, num_impute = 20, method = "CCMV", j_0 = NULL, xi = NULL) {
  out <- list()
  iters <- floor(seq(from = 1, to = length(chain$loglik), length = num_impute))

  if(is.null(j_0) & method == "NIP") stop("Must include time to impute data for NIP")
  if(is.null(j_0) & method == "TLO") stop("Must include time to impute data for TLO")
  if(!is.null(chain$gamma) & method == "MAR") stop("This looks like a chain for MNAR, but using MAR imputation; fit the model under MAR to use this")
  if(is.null(xi)) xi <- rep(0, num_impute)
  stopifnot(method %in% c("CCMV", "TLO", "NIP", "MAR"))
  stopifnot(is.numeric(xi)); stopifnot(length(xi) == num_impute)

  for(n in 1:num_impute) {
    i <- iters[n]
    omega <- chain$omega[i, ]
    log_omega <- log(omega)
    gamma <- chain$gamma[i,,]
    if(method != "MAR") {
      log_gamma <- log(gamma)
      log_1_m_gamma <- log(1 - gamma)
    }
    beta <- chain$beta[i,,]
    log_beta <- log(beta)
    log_1_m_beta <- log(1 - beta)

    if(method == "CCMV") {
      out[[n]] <- CCMVMI(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta)
    } else if(method == "NIP") {
      out[[n]] <- NIPMI(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, j_0)
    } else if(method == "TLO") {
      out[[n]] <- TLOMI(Y, R, omega, log_omega, gamma, log_gamma, log_1_m_gamma, beta, log_beta, log_1_m_beta, xi[n], j_0)
    } else {
      out[[n]] <- MARMI(Y = Y, R = R, omega = omega, log_omega = log_omega, beta = beta, log_beta = log_beta, log_1_m_beta = log_1_m_beta)
    }
  }
  return(out)
}
