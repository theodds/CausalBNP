#' Extracts posterior samples of the marginal means
#' 
#' Computes the marginal mean for each time point, under the MAR assumption. 
#' 
#' @param chain An object of type parafacmar, obtained from fitting the PARAFAC 
#' model under an ignorable MAR assumption
#' 
#' @return Returns an N x J matrix where N is the number of posterior samples 
#' collected and J is the number of time points. The entry (i,j) corresponds to 
#' the sample of the marginal mean at iteration i and time j.

ExtractMARMeans <- function(chain) {
  
  stopifnot(class(chain) == "parafacmar")

  N_iter <- length(chain$omega[,1])

  f <- function(i) {
    t(chain$beta[i,,]) %*% chain$omega[i,]
  }

  out <- t(sapply(1:N_iter, FUN = f))

  return(out)

}
