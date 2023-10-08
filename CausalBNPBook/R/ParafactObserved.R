#' Compute observed data statistics for nonignorable parafac model
#'
#' Computes posterior mean and credible intervals for the marginal means of the
#' observed data model, and compres to nonparametric intervals. This gives a
#' sanity check on the model; if the mode disagrees substantially with the
#' nonparametric fit, this suggests potential lack of fit.
#'
#' @param chain an object of type parafac, obtained from fitting the nonignorable parafac model
#' @param Y binary NxJ matrix of observed responses.
#' @param R binary NxJ matrix of missingness indicators
#' @param do_plot boolean; if true, creates a diagnostic plot checking that the prior reproduces nonparametric inference on the identified marginal parameters. The dark error bars give 95% credible intervals for the parameters, while orange bars give 95 percent nonparametric confidence intervals.
#'
#' @return Returns a list consisting of
#' \itemize{
#' \item R_means: IxJ matrix, where (i,j)th entry is the mean of R_j at iteration i.
#' \item Y_means: IxJ matrix, where (i,j)th entry is the mean of the observed Y_j at iteration i.
#' \item Y_posterior_mean: J component vector, with jth entry given by the posterior mean of the observed data mean at time j.
#' \item Y_lcl: Lower confidence limit (95\%) for the observed data posterior mean at time j.
#' \item Y_ucl: Upper confidence limit (95\%) for the observed data posterior mean at time j.
#' \item R_posterior_mean: J component vector, with jth entry given by the posterior mean of the probability of missingness at time j.
#' \item R_ucl: See Y_ucl, but for R.
#' \item R_lcl: See R_lcl, but for R.
#' \item Y_obs: The empirical observed data means.
#' \item Y_obs_lcl: Nonparametric lower confidence limit (95\%) for the observed data mean at time j.
#' \item Y_obs_ucl: Nonparametric upper confidence limit (95\%) for the observed data mean at time j.
#' \item N_obs: J vector giving the number of observed responses at time j.
#' }

ParafacObserved <- function(chain, Y, R, do_plot = TRUE) {

  stopifnot(class(chain) == "parafac")

  ## Mean for each R_j:
  ## sum_k omega_k gamma_{kj}

  N_iter <- length(chain$omega[,1])
  R_means <- t(sapply(1:N_iter, function(i) t(chain$gamma[i,,]) %*% chain$omega[i,]))

  R_posterior_mean <- colMeans(R_means)
  R_lcl <- apply(R_means, 2, function(x) quantile(x, .025))
  R_ucl <- apply(R_means, 2, function(x) quantile(x, 0.975))

  ## Mean for each Y_j
  ## sum_k omega_k gamma_{kj} beta_{kj} / sum_k omega_k gamma_{kj}

  Y_means <- t(sapply(1:N_iter, function(i) t(chain$gamma[i,,] * chain$beta[i,,]) %*% chain$omega[i,]))
  Y_means <- Y_means / R_means
  Y_posterior_mean <- colMeans(Y_means)
  Y_lcl <- apply(Y_means, 2, function(x) quantile(x, 0.025))
  Y_ucl <- apply(Y_means, 2, function(x) quantile(x, 0.975))



  Y_obs <- colMeans(Y, na.rm = TRUE)
  N_obs <- apply(Y, MARGIN = 2, FUN = function(x) sum(!is.na(x)))
  sd_obs <- apply(Y, MARGIN = 2, FUN = function(x) sd(x, na.rm = TRUE))
  sd_obs <- sd_obs / sqrt(N_obs)

  Y_obs_lcl <- Y_obs - qnorm(.975) * sd_obs
  Y_obs_ucl <- Y_obs + qnorm(.975) * sd_obs
  R_obs <- colMeans(R)
  R_obs_lcl <- R_obs - qnorm(.975) * sqrt(R_obs * (1 - R_obs) / N_obs)
  R_obs_ucl <- R_obs + qnorm(.975) * sqrt(R_obs * (1 - R_obs) / N_obs)

  if(do_plot) {

    time_points <- 1:ncol(Y)

    p <- ggplot()
    p <- p + geom_line(aes(x = time_points, y = Y_obs), color = 'darkorange1')
    p <- p + geom_point(aes(x = time_points, y = Y_posterior_mean))
    p <- p + geom_errorbar(aes(x = time_points, ymax = Y_ucl, ymin = Y_lcl), width = .3)
    p <- p + geom_errorbar(aes(x = time_points + .3, ymax = Y_obs_ucl, ymin = Y_obs_lcl), width = .15, color = 'darkorange1')
    p <- p + theme_bw()

    p_2 <- ggplot()
    p_2 <- p_2 + geom_line(aes(x = time_points, y = R_obs), color = 'darkorange1')
    p_2 <- p_2 + geom_point(aes(x = time_points, y = R_posterior_mean))
    p_2 <- p_2 + geom_errorbar(aes(x = time_points, ymax = R_ucl, ymin = R_lcl), width = .3)
    p_2 <- p_2 + geom_errorbar(aes(x = time_points + .3, ymax = R_obs_ucl, ymin = R_obs_lcl), width = .15, color = 'darkorange1')
    p_2 <- p_2 + theme_bw()

    grid.arrange(p, p_2, nrow = 1)

    # plot(Y_posterior_mean, pch = 1)
    # points(colMeans(Y, na.rm = TRUE), pch = 10, col = 'blue')

  }

  return(list(R_means = R_means,
              Y_means = Y_means,
              Y_posterior_mean = Y_posterior_mean,
              Y_lcl = Y_lcl,
              Y_ucl = Y_ucl,
              R_posterior_mean = R_posterior_mean,
              R_lcl = R_lcl,
              R_ucl = R_ucl,
              Y_obs = colMeans(Y, na.rm = TRUE),
              N_obs = apply(Y, MARGIN = 2, FUN = function(x) sum(!is.na(x))),
              Y_obs_lcl = Y_obs_lcl,
              Y_obs_ucl = Y_obs_ucl
              ))

}
