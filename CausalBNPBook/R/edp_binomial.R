#' Enriched Dirichlet Process for Binary Outcome
#'
#' This function implements an enriched Dirichlet process for observational studies
#' with a binary outcome, possibly missing covariates.
#'
#' @param Y Numeric vector representing the binary outcome variable.
#' @param A Numeric vector indicating the treatment assignment.
#' @param L_binary Data frame or matrix of binary covariates.
#' @param L_continuous Data frame or matrix of continuous covariates.
#' @param R_binary Data frame or matrix with missingness indicators for binary covariates.
#' @param R_continuous Data frame or matrix with missingness indicators for continuous covariates.
#' @param num_burn Integer. Number of burn-in iterations. Default is 1000.
#' @param num_thin Integer. Thinning rate for MCMC. Default is 1.
#' @param num_save Integer. Number of iterations to save after burn-in. Default is 1000.
#' @param num_mc Integer. Number of Monte Carlo simulations. Default is 2000.
#' @param thin_mc_int Integer. Interval for thinning in Monte Carlo simulations. Default is 100.
#' @param M Integer. Number of clusters for the Dirichlet Process. Default is 1000.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{summary}}{ A data frame with summary statistics for the variables:
#'     \code{alpha_theta}, \code{alpha_omega}, and \code{psi_rr}.}
#'   \item{\code{alpha_theta}}{ The draws of \code{alpha_theta} post-burnin.}
#'     \item{\code{alpha_omega}}{ The draws of \code{alpha_omega} post-burnin.}
#'     \item{\code{psi_rr}}{The draws of \code{psi_rr}, the causal relative risk,
#'     from the chain post-burnin.}
#' }
#' The variable \code{psi_rr} denotes the causal relative risk, which is the object of primary interest, whereas \code{alpha_omega} and \code{\alpha_theta} are the concentration parameters from the enriched Dirichlet proces.
#'
edp_binomial <- function(Y, A, L_binary, L_continuous, R_binary, R_continuous,
                         num_burn = 1000,
                         num_thin = 1,
                         num_save = 1000,
                         num_mc = 2000,
                         thin_mc_int = 100,
                         M = 1000) {

  stopifnot(is.matrix(A))
  stopifnot(is.matrix(L_binary))
  stopifnot(is.matrix(L_continuous))

  gibbs_total <- num_burn + num_thin * num_save
  MCnum <- num_mc
  n <- length(Y)
  burn_in <- num_burn

  # Define number of binary covariates
  p1 <- ncol(L_binary)

  # Define number of continuous covariates
  p2 <- ncol(L_continuous)

  # Define number of treatment variables
  ptx <- ncol(A)

  # Calculate number of covariates
  p <- 1 + ptx + p1 + p2

  ## Define L matrix: binary columns come first
  L     <- cbind(L_binary, L_continuous)
  L_mar <- cbind(R_binary, R_continuous)

  ## Define Functions
  rinvchisq <- function(n, df, scale = 1/df) {
    df <- rep(df, len = n)
    scale <- rep(scale, len = n)
    if (any(df <= 0))
      stop("The df parameter must be positive.")
    if (any(scale <= 0))
      stop("The scale parameter must be positive.")
    x <- (df * scale)/rchisq(n, df = df)
    return(x)
  }

  # Define function to calculate inverse of a matrix by calling chol2inv to
  # invert the matrix from its Choleski decomposition
  # inv <- function(x) {
  #   chol2inv(chol(x))
  # }

  # Define function to update beta coefficients in regression. Input current
  # value of beta coefficient, a vector of predictors(including intercept)
  # and outcome.  Currently uses Normal(mean=0,sd=2) prior for beta
  # coefficients. Change this within the function if needed.
  ## FLAG
  update_beta_coeff <- function(beta_coeff_current, x_temp, y_temp) {

    proposed <- as.numeric(rmvn(1, beta_coeff_current, diag(c(rep(0.01, p)))))

    loglike_prop <- sum(dbinom(y_temp, 1, prob = plogis(x_temp %*% proposed),
                               log = TRUE)) + dmvn(proposed, beta_coeff_init, (diagbetacov0), log = TRUE) +
      dmvn(beta_coeff_current, proposed, diag(c(rep(0.01, p))), log = TRUE)

    loglike_old <- sum(dbinom(y_temp, 1, prob = plogis(x_temp %*% beta_coeff_current),
                              log = TRUE)) + dmvn(beta_coeff_current, beta_coeff_init, (diagbetacov0),
                                                  log = TRUE) + dmvn(proposed, beta_coeff_current, diag(c(rep(0.01, p))),
                                                                     log = TRUE)

    # prevent program from quitting
    loglike_prop <- replace(loglike_prop, which(loglike_prop == -Inf), -999999)
    loglike_old <- replace(loglike_old, which(loglike_old == -Inf), -999999)


    ans <- min(exp(loglike_prop - loglike_old), 1)

    beta_coeff_draw <- rep(0, p)

    if (runif(1, 0, 1) < ans)
      beta_coeff_draw <- proposed else beta_coeff_draw <- beta_coeff_current

    return(beta_coeff_draw)
  }


  # Define function to update variance parameter for covariate, N-Inv-chi^2
  # model
  update_variance <- function(x_temp) {

    df_new <- nu0 + length(x_temp)

    if (length(x_temp) == 1) {
      varx <- 0
    }

    if (length(x_temp) > 1) {
      varx <- var(x_temp)
    }

    numer <- nu0 * tau0 + (length(x_temp) - 1) * varx + (c0 * length(x_temp)/(c0 +
                                                                                length(x_temp))) * (mean(x_temp) - mu0)^2

    newval <- rinvchisq(1, df_new, numer/df_new)

    return(newval)
  }

  # Define function to update mean parameter for covariate
  update_mean <- function(x_temp, tau_current) {
    newvar <- 1/(c0/tau_current + length(x_temp)/tau_current)
    newmean <- (mu0 * c0/tau_current + mean(x_temp) * length(x_temp)/tau_current) * newvar
    newval <- rnorm(1, newmean, sqrt(newvar))
    return(newval)
  }

  # Define function to update alpha_theta prior Gamma(alpha_a0,alpha_b0)
  update_alpha_theta <- function(number_of_clusters, alpha_current) {
    eta <- rbeta(1, alpha_current + 1, n)
    pieta <- (number_of_clusters/(n * (1 - log(eta))))/(1 + number_of_clusters/(n *
                                                                                  (1 - log(eta))))
    whichmix <- rbinom(1, 1, pieta)
    alpha_new <- whichmix * rgamma(1, (alpha_a0 + number_of_clusters), (alpha_b0 -
                                                                          log(eta))) + (1 - whichmix) * rgamma(1, (alpha_a0 + number_of_clusters - 1),
                                                                                                               (alpha_b0 - log(eta)))
    return(alpha_new)
  }

  # Define function to update alpha_omega using Metropolis Hastings
  update_alpha_omega <- function(alpha_current) {
    sort_unique_y <- sort(unique(s[, 1]))
    ss <- numeric(length(sort_unique_y))
    for (j in 1:length(sort_unique_y)) {
      ss[j] <- sum(s[, 1] == sort_unique_y[j])
    }
    like_current <- dgamma(alpha_current, alpha_a0, alpha_b0) * (alpha_current^(nrow(unique(s)) - k)) *
      prod((alpha_current + ss) * beta(alpha_current + 1, ss))
    proposed <- rgamma(1, 2, 1)
    like_proposed <- dgamma(proposed, alpha_a0, alpha_b0) * (proposed^(nrow(unique(s)) - k)) *
      prod((proposed + ss) * beta(proposed + 1, ss))
    ratio <- (like_proposed)/(like_current)
    alpha_new <- ifelse((runif(1, 0, 1) < ratio), proposed, alpha_current)

    return(alpha_new)
  }


  # Standardize continuous covariates (important when choosing priors)
  L[, (p1 + 1):(p1 + p2)] <- scale(L[, (p1 + 1):(p1 + p2)])

  # Define design matrix
  X <- cbind(A, L)
  X_matrix <- cbind(1, X)


  # Set initial values
  # ------------------------------------------------------------

  # Parameters in outcome regression
  ## FLAG
  beta_coeff_init <- glm(Y ~ X, family = "binomial")$coef
  # tau_init <- 1

  # Hyperparameters

  # for modeling binary covariates
  a0 <- 1
  b0 <- 1

  # for modeling continuous covariates
  nu0 <- 2
  tau0 <- 1
  c0 <- 0.5
  mu0 <- 0

  # for outcome regression coefficients
  ## FLAG
  # bayesreg <- glm(Y ~ X, family = "binomial")
  # betacov0 <- vcov(bayesreg)
  diagbetacov0 <- diag(c(rep(4, p)))

  # for concentration parameters
  alpha_a0 <- 1
  alpha_b0 <- 1

  # Set initial values for cluster membership.
  # Use k-means to find 3 whole data clusters
  s <- matrix(nrow = n, ncol = 2, 1)
  s[, 1] <- kmeans(cbind(Y, A, L), 3)$cluster

  # Make 2 x-clusters within each y-cluster
  for (j in 1:3) {
    s[s[, 1] == j, 2] <- kmeans(cbind(A[s[, 1] == j], L[s[, 1] == j, ]), 2)$cluster
  }

  # an alternate clustering method
  distyx <- dist(cbind(Y, A, L))
  hi <- hclust(distyx, method = "ward.D")
  s[, 1] <- cutree(hi, k = 3)

  for (j in 1:3) {
    distx <- dist(cbind(A[s[, 1] == j], L[s[, 1] == j, ]))
    hi <- hclust(distx, method = "ward.D")
    groups <- cutree(hi, k = 2)
    s[s[, 1] == j, 2] <- groups
  }


  # Calculate number of y-clusters
  k <- length(unique(s[, 1]))

  # Update values of beta_coeffs for outcome regressions
  beta_coeff <- matrix(nrow = k, ncol = ncol(X_matrix))
  for (i in 1:k) {
    beta_coeff_current <- as.numeric(rmvn(1, beta_coeff_init, (diagbetacov0)))
    beta_coeff[i, ] <- update_beta_coeff(beta_coeff_current, X_matrix[s[, 1] ==
                                                                        unique(s[, 1])[i], ], Y[s[, 1] == unique(s[, 1])[i]])

  }

  # Update covariate parameters
  # Within each cluster parameters are stored in a long vector nk x 1
  # Binary covariates have 1 parameter (x_prob_param)
  # Continuous covariates have 2 parameters (x_mean_param, x_var_param)

  # Initialize vectors
  # nk is the total number of clusters (x and y clusters)
  nk <- nrow(unique(s))
  x_prob_param <- matrix(nrow = nk, ncol = (p1 + ptx))

  if (p2 > 0) {
    x_mean_param <- matrix(nrow = nk, ncol = p2)
    x_var_param <- matrix(nrow = nk, ncol = p2)
  }

  # Initialize matrices for use in Monte Carlo integration
  max.num.xclusters <- max(s[, 2])
  x_prob_param.mat <- array(0, c((k + 1), (max.num.xclusters + 1), ptx + p1))

  if (p2 > 0) {
    x_mean_param.mat <- array(0, c((k + 1), (max.num.xclusters + 1), p2))
    x_var_param.mat <- array(0, c((k + 1), (max.num.xclusters + 1), p2))
  }

  # Update parameter for binary covariates
  uniquey <- unique(s[, 1])
  sorteduniquey <- uniquey[order(uniquey)]
  for (i in 1:(ptx + p1)) {
    # beta(1,1) prior
    count <- 1
    for (j in 1:k) {
      for (l in 1:length(unique(s[s[, 1] == sorteduniquey[j], 2]))) {
        x_temp <- X[(s[, 1] == sorteduniquey[j]), i]
        s_temp <- s[(s[, 1] == sorteduniquey[j]), 2]
        x_temp <- x_temp[s_temp == sort(unique(s_temp))[l]]
        # posterior is beta
        x_prob_param.mat[j, l, i] <- rbeta(1, (sum(x_temp) + a0), (length(x_temp) -
                                                                     sum(x_temp) + b0))
        x_prob_param[count, i] <- x_prob_param.mat[j, l, i]
        count <- count + 1
      }
    }
  }

  # update parameters for continuous covariates
  if (p2 > 0) {
    for (i in 1:p2) {
      # beta(1,1) prior
      count <- 1
      for (j in 1:k) {
        for (l in 1:length(unique(s[s[, 1] == sorteduniquey[j], 2]))) {
          x_temp <- X[(s[, 1] == sorteduniquey[j]), (p1 + ptx + i)]
          s_temp <- s[(s[, 1] == sorteduniquey[j]), 2]
          x_temp <- x_temp[s_temp == sort(unique(s_temp))[l]]

          x_var_param.mat[j, l, i] <- update_variance(x_temp)
          x_var_param[count, i] <- x_var_param.mat[j, l, i]

          x_mean_param.mat[j, l, i] <- update_mean(x_temp, x_var_param.mat[j, l, i])
          x_mean_param[count, i] <- x_mean_param.mat[j, l, i]
          count <- count + 1
        }
      }
    }
  }

  # Calculate h0i and h0y for use in cluster function in Gibbs Sampler
  # Average covariate distribution over prior for each x_i
  h0x <- matrix(nrow = n, ncol = ncol(X), 0)
  denom <- (sqrt(2 * pi)/sqrt(c0)) * gamma(nu0/2) * (2/(nu0 * tau0))^(nu0/2)
  cn <- c0 + 1
  nun <- nu0 + 1
  # Continuous covariates
  for (i in 1:n) {
    for (j in 1:p2) {
      x <- X[i, (p1 + ptx + j)]
      taun <- (1/nun) * (nu0 * tau0 + (c0/(c0 + 1)) * (mu0 - x)^2)
      num <- (sqrt(2 * pi)/sqrt(cn)) * gamma(nun/2) * (2/(nun * taun))^(nun/2)
      h0x[i, (p1 + ptx + j)] <- num/(denom * sqrt(2 * pi))  #sqrt(2*pi) is from the data part of likelihood
    }
  }
  # Binary covariates (including treatment) Beta-Binomial
  for (i in 1:n) {
    for (j in 1:(p1 + ptx)) {
      h0x[i, j] <- beta(a0 + X[i, j], b0 - X[i, j] + 1)
    }
  }

  # Initialize
  h0xa1 <- h0x
  h0xa0 <- h0x

  # Integrate Kernel for treatment variable over distribution of probability parameter
  # Approximate integral with average
  # For treatment variable, kernel is binomial(pi). Prior for pi is beta(a0,b0)
  # h0xa1 is when treatment is set to 1
  # h0xa0 is when treatment is set to 0
  for (i in 1:n) {
    for (j in 1:(ptx)) {
      piprior <- rbeta(100000, a0, b0)
      h0xa1[i, j] <- mean(dbinom(1, 1, piprior))
      h0xa0[i, j] <- mean(dbinom(0, 1, piprior))
    }
  }

  # Take product (covariates are assumed to be locally independent).
  # Result is vector of size n
  h0i <- apply(h0x, 1, prod)
  # h0ia1 <- apply(h0xa1, 1, prod)
  # h0ia0 <- apply(h0xa0, 1, prod)

  # Find E0(Y|A,L,beta_coeff) integrating over prior for beta coefficient
  # using Monte Carlo integration
  ## FLAG
  beta_coeff_prior <- beta_coeff_init
  h0ya1 <- numeric(n)
  h0ya0 <- numeric(n)
  h0y <- numeric(n)
  for (i in 1:100000) {
    beta_coeff_prior <- as.numeric(rmvn(1, beta_coeff_init, (diagbetacov0)))

    h0ya1 <- h0ya1 + plogis(cbind(1, 1, X[, ((ptx + 1):(ptx + p1 + p2))]) %*%
                              beta_coeff_prior)
    h0ya0 <- h0ya0 + plogis(cbind(1, 0, X[, ((ptx + 1):(ptx + p1 + p2))]) %*%
                              beta_coeff_prior)
    h0y <- h0y + dbinom(Y, 1, plogis(X_matrix %*% beta_coeff_prior))
  }
  h0ya1 <- h0ya1/100000
  h0ya0 <- h0ya0/100000
  h0y <- h0y/100000

  # Will need the following for calculating causal effects using MC
  # integration.  Any calculations involving X are done in the Gibbs Sampler
  h0xa1.mc <- matrix(nrow = M, ncol = (ptx + p1 + p2), 0)
  h0xa0.mc <- matrix(nrow = M, ncol = (ptx + p1 + p2), 0)

  for (i in 1:M) {
    for (j in 1:ptx) {
      piprior <- rbeta(MCnum, a0, b0)
      h0xa1.mc[i, j] <- mean(dbinom(1, 1, piprior))
      h0xa0.mc[i, j] <- mean(dbinom(0, 1, piprior))
    }
  }

  h0ya1.mc <- numeric(M)
  h0ya0.mc <- numeric(M)

  beta_coeff_prior.mc <- matrix(nrow = MCnum, ncol = p, 0)
  for (i in 1:MCnum) {
    beta_coeff_prior.mc[i, ] <- as.numeric(rmvn(1, beta_coeff_init, (diagbetacov0)))
  }


  # Initialize alpha parameters
  alpha_theta <- 2
  alpha_omega <- 2

  # Make vectors to store draws from Gibbs Sampler
  psi_rr <- numeric(gibbs_total)
  psi_rd <- numeric(gibbs_total)

  alpha_theta_draws <- numeric(gibbs_total)
  alpha_omega_draws <- numeric(gibbs_total)

  # Start Gibbs Sampler
  # -------------------------------------------------------
  for(gibbsreps in seq_len(gibbs_total)) {
    uniques <- unique(s)
    # Sort clusters
    sortedunique <- uniques[order(uniques[, 1], uniques[, 2]), , drop = FALSE]

    # Update cluster membership --------------------------------------------
    cluster_res <- cluster(Y, X, X_matrix, s[, 1], s[, 2], beta_coeff, x_prob_param,
                           x_mean_param, x_var_param, alpha_omega, alpha_theta, h0y, h0i, sortedunique,
                           c0, mu0, nu0, tau0, a0, b0, beta_coeff_init, diagbetacov0, p1, ptx, p2)

    # Store cluster membership output from cluster function
    s <- cbind(cluster_res$Sy, cluster_res$Sx)
    beta_coeff <- cluster_res$beta
    x_prob_param <- cluster_res$pipars
    x_mean_param <- cluster_res$mupars
    x_var_param <- cluster_res$sigpars

    # Make vector of y-clusters
    uniquey <- unique(s[, 1])
    sorteduniquey <- uniquey[order(uniquey)]

    # Calculate number of y-clusters
    k <- length(uniquey)

    # Make matrix of clusters
    uniques <- unique(s)

    # Find largest number of x clusters
    max.num.xclusters <- max(s[, 2])

    # Store number of x-clusters in each y-cluster in a vector
    maxslj <- rep(0, k)
    for (j in 1:k) {
      maxslj[j] <- sum(uniques[, 1] == sorteduniquey[j])
    }

    # Calculate number of subjects in each y-cluster and store in vector nj.
    # Calculate number of subjects in each x-cluster and store in matrix nlj.
    # Use j to index y-clusters and l to index x-clusters.
    nj <- rep(0, k)
    nlj <- matrix(0, nrow = k, ncol = (max.num.xclusters))
    for (j in 1:k) {
      nj[j] <- sum(s[, 1] == sorteduniquey[j])
      s_temp2 <- s[s[, 1] == sorteduniquey[j], 2]
      s_temp2 <- sort(unique(s_temp2))
      for (l in 1:maxslj[j]) {
        nlj[j, l] <- sum(s[, 1] == sorteduniquey[j] & s[, 2] == s_temp2[l])
      }
    }

    # End update of cluster membership
    # --------------------------------------------

    ## FLAG
    # Update values of beta coefficients in outcome regressions-------------
    for (j in 1:k) {
      beta_coeff_current <- beta_coeff[j, ]
      beta_coeff[j, ] <- update_beta_coeff(beta_coeff_current,
                                           X_matrix[s[, 1] == j, ],
                                           Y[s[, 1] == j])
    }
    # End update of beta coefficients in outcome regressions-------------------

    # Update cluster specific parameters for covariate distributions.----------
    # Binary covariates have 1 parameter (x_prob_param).
    # Continuous covariates have 2 parameters (x_mean_param,x_var_param).
    # Also store parameters in a matrix for use in MC integration later.
    # Vectors of parameters(x_prob_param, x_var_param, x_mean_param) is for cluster function
    x_prob_param.mat <- array(0, c((k + 1), (max.num.xclusters + 1), ptx + p1))

    if (p2 > 0) {
      x_mean_param.mat <- array(0, c((k + 1), (max.num.xclusters + 1), p2))
      x_var_param.mat <- array(0, c((k + 1), (max.num.xclusters + 1), p2))
    }

    # Update parameters for binary covariates
    # prior is beta(a0,b0), posterior is beta
    for (i in 1:(ptx + p1)) {
      count <- 1
      for (j in 1:k) {
        for (l in 1:maxslj[j]) {
          x_temp <- X[(s[, 1] == sorteduniquey[j]), i]
          s_temp <- s[(s[, 1] == sorteduniquey[j]), 2]
          x_temp <- x_temp[s_temp == sort(unique(s_temp))[l]]

          x_prob_param.mat[j, l, i] <- rbeta(1, sum(x_temp) + a0, length(x_temp) -
                                               sum(x_temp) + b0)
          x_prob_param[count, i] <- x_prob_param.mat[j, l, i]
          count <- count + 1
        }
      }
    }

    # Update parameters for continuous predictors using functions
    if (p2 > 0) {
      for (i in 1:p2) {
        count <- 1
        for (j in 1:k) {
          for (l in 1:maxslj[j]) {
            x_temp <- X[(s[, 1] == sorteduniquey[j]), (p1 + ptx + i)]
            s_temp <- s[(s[, 1] == sorteduniquey[j]), 2]
            x_temp <- x_temp[s_temp == sort(unique(s_temp))[l]]

            x_var_param.mat[j, l, i] <- update_variance(x_temp)
            x_var_param[count, i] <- x_var_param.mat[j, l, i]

            # posterior for mu. prior for mu|sigma^2 mean 0 prior sample size 2
            x_mean_param.mat[j, l, i] <- update_mean(x_temp, x_var_param.mat[j,
                                                                             l, i])
            x_mean_param[count, i] <- x_mean_param.mat[j, l, i]
            count <- count + 1
          }
        }
      }
    }
    # End update of cluster-specific parameters for covariates----------

    # Update concentration parameters----------------------------------
    alpha_theta <- update_alpha_theta(length(unique(s[, 1])), alpha_theta)
    alpha_theta_draws[gibbsreps] <- alpha_theta

    alpha_omega <- update_alpha_omega(alpha_omega)
    alpha_omega_draws[gibbsreps] <- alpha_omega

    # End update of concentration parameters--------------------------

    sorteduniquey <- sort(unique(s[, 1]))

    uniquey <- unique(s[, 1])
    uniques <- unique(s)
    # nunique <- nrow(uniques)
    # sortuniquey <- sort(uniquey)
    sortuniqueyx <- uniques[order(uniques[, 1], uniques[, 2]), , drop = FALSE]
    # nuniquey <- length(sortuniquey)


    # End update of all parameters in BNP model--------------------------------

    # Update covariates with missing data-------------------------------------

    # Find the row of the sorted cluster matrix that each subject is a member of.
    rowx <- numeric(n)
    for (i in 1:n) {
      rowx[i] <- which(s[i, 1] == sortuniqueyx[, 1] & s[i, 2] == sortuniqueyx[, 2])
    }

    # Binary covariates: assume binomial distribution
    if(p1>0){
      for (r in 1:p1) {
        # Find missing observations
        tempx1 <- X_matrix[L_mar[, r] == 1, ]
        nmiss <- dim(tempx1)[1]

        # Calculate probability that missing binary covariate is 1.
        # Set missing covariate to 0 and 1.
        ## FLAG
        tempx0 <- tempx1
        tempx1[, (r + 1 + ptx)] <- 1
        tempx0[, (r + 1 + ptx)] <- 0
        tempprobs <- cbind((1 - x_prob_param[rowx[L_mar[, r] == 1], (r + ptx)]) *
                             plogis(rowSums(tempx0 * (beta_coeff[s[L_mar[, r] == 1, 1], ]))),
                           x_prob_param[rowx[L_mar[, r] == 1], (r + ptx)] * plogis(rowSums(tempx1 *
                                                                                             (beta_coeff[s[L_mar[, r] == 1, 1], ]))))
        # Draw missing values from binomial distribution
        L[L_mar[, r] == 1, r] <- rbinom(nmiss, 1, tempprobs[, 2]/rowSums(tempprobs))
      }
    }

    # Continuous covariates are updated using Metropolis Hastings
    if(p2>0){
      for (r in (p1 + 1):(p1 + p2)) {
        proposed <- rnorm(sum(L_mar[, r] == 1), L[L_mar[, r] == 1, r], 0.1)
        propX_matrix <- X_matrix[L_mar[, r] == 1, ]
        propX_matrix[, (r + 1 + ptx)] <- proposed

        loglikeprop <- dnorm(proposed, x_mean_param[rowx[L_mar[, r] == 1],
                                                    (r - p1)], sqrt(x_var_param[rowx[L_mar[, r] == 1], (r - p1)]),
                             log = TRUE) + dbinom(Y[L_mar[, r] == 1], 1, plogis(rowSums(propX_matrix *
                                                                                          (beta_coeff[s[L_mar[, r] == 1, 1], ]))), log = TRUE) + dnorm(L[L_mar[,
                                                                                                                                                               r] == 1, r], proposed, 0.1, log = TRUE)

        loglikeold <- dnorm(L[L_mar[, r] == 1, r], x_mean_param[rowx[L_mar[,
                                                                           r] == 1], (r - p1)], sqrt(x_var_param[rowx[L_mar[, r] == 1], (r -
                                                                                                                                           p1)]), log = TRUE) + dbinom(Y[L_mar[, r] == 1], 1, plogis(rowSums(X_matrix[L_mar[,
                                                                                                                                                                                                                            r] == 1, ] * (beta_coeff[s[L_mar[, r] == 1, 1], ]))), log = TRUE) +
          dnorm(proposed, L[L_mar[, r] == 1, r], 0.1, log = TRUE)

        # Prevent R from quitting
        loglikeprop <- replace(loglikeprop, which(loglikeprop == -Inf), -999999)
        loglikeold <- replace(loglikeold, which(loglikeold == -Inf), -999999)

        indswitch <- runif(sum(L_mar[, r] == 1), 0, 1) < exp(loglikeprop -
                                                               loglikeold)

        L[L_mar[, r] == 1, r] <- L[L_mar[, r] == 1, r] * (1 - indswitch) +
          indswitch * proposed
      }
    }

    # End update of missing covariates------------------------------


    # Now update variables that depend on updated covariates-----------------------

    # Update X
    # -------------------------------------------------------------------------

    X <- cbind(A, L)
    X_matrix <- cbind(1, X)

    # Calculate h0i and h0y for use in cluster function in Gibbs Sampler
    # Average covariate distribution over prior for each x_i

    # Initialize
    h0x <- matrix(nrow = n, ncol = ncol(X), 0)
    denom <- (sqrt(2 * pi)/sqrt(c0)) * gamma(nu0/2) * (2/(nu0 * tau0))^(nu0/2)
    cn <- c0 + 1
    nun <- nu0 + 1
    # Continuous covariates
    for (i in 1:n) {
      for (r in 1:p2) {
        x <- X[i, (p1 + ptx + r)]
        taun <- (1/nun) * (nu0 * tau0 + (c0/(c0 + 1)) * (mu0 - x)^2)
        num <- (sqrt(2 * pi)/sqrt(cn)) * gamma(nun/2) * (2/(nun * taun))^(nun/2)
        h0x[i, (p1 + ptx + r)] <- num/(denom * sqrt(2 * pi))  #sqrt(2*pi) is from the data part of likelihood
      }
    }

    # Binary covariates (including treatment) Beta-Binomial
    for (i in 1:n) {
      for (r in 1:(p1 + ptx)) {
        h0x[i, r] <- beta(a0 + X[i, r], b0 - X[i, r] + 1)
      }
    }

    # h_0(x_i) is n by 1
    h0i <- apply(h0x, 1, prod)


    # Find E0(Y|A,L,beta_coeff) integrating over prior for beta coefficient
    # using Monte Carlo integration
    beta_coeff_prior <- beta_coeff_init
    h0ya1 <- numeric(n)
    h0ya0 <- numeric(n)
    h0y <- numeric(n)
    for (i in 1:MCnum) {
      beta_coeff_prior <- as.numeric(rmvn(1, beta_coeff_init, (diagbetacov0)))
      h0ya1 <- h0ya1 + plogis(cbind(1, 1, X[, (2:(ptx + p1 + p2))]) %*% beta_coeff_prior)
      h0ya0 <- h0ya0 + plogis(cbind(1, 0, X[, (2:(ptx + p1 + p2))]) %*% beta_coeff_prior)
      h0y <- h0y + dbinom(Y, 1, plogis(X_matrix %*% beta_coeff_prior))
    }
    h0ya1 <- h0ya1/MCnum
    h0ya0 <- h0ya0/MCnum
    h0y <- h0y/MCnum

    # Calculate causal effect using Monte Carlo Integration to integrate out
    # confounders.  Don't need to do every iteration so compute every 100th
    # draw after burn_in)
    # Start----------------------------------------------------------------------

    ##FLAG: Need to control how Y is simulated, so this seems important
    if(gibbsreps >= num_burn & gibbsreps %% thin_mc_int) {

      # Initialize matrices and vectors
      beta_coeff.mc <- matrix(nrow = (k + 1), ncol = ncol(X_matrix))
      beta_coeff.mc[1:k, ] <- beta_coeff
      s.mc <- matrix(1, nrow = M, ncol = 2)
      L.mc <- matrix(NA, nrow = M, ncol = (p1 + p2))
      prob_ycluster <- rep(NA, (k + 1))

      # Assign M observations one of the current y-clusters or to a new y-cluster
      # using draws from a multinomial distribution with probabilities
      # 'prob_ycluster'

      pdenom <- n + alpha_theta
      for (j in 1:k) {
        prob_ycluster[j] <- nj[j]/pdenom
      }
      prob_ycluster[(k + 1)] <- alpha_theta/pdenom

      s.mc[, 1] <- sample(1:(k + 1), size = M, replace = T, prob = prob_ycluster)

      # Store indices of non-empty y-clusters (some y-clusters could be empty and
      # there could be a new cluster)
      yclusters.mc <- sort(unique(s.mc[, 1]))

      # Store number of y clusters including empty ones and possibly a new one
      k.mc <- max(yclusters.mc)

      # Store number of observations in each y-cluster (should add to a total of
      # M)
      Nj.mc <- rep(0, k.mc)

      # If a new y-cluster was opened, draw parameters for this cluster from
      # priors.

      if (k.mc > k) {
        for (rep in 1:(ptx + p1)) {
          x_prob_param.mat[k.mc, 1, rep] <- rbeta(n = 1, shape1 = a0, shape2 = b0)
        }

        if (p2 > 0) {
          for (rep in 1:p2) {
            x_var_param.mat[k.mc, 1, rep] <- rinvchisq(n = 1, df = nu0, scale = tau0)
            x_mean_param.mat[k.mc, 1, rep] <- rnorm(n = 1, mean = mu0,
                                                    sd = sqrt(1/c0 * x_var_param.mat[k.mc, 1, rep]))
          }
        }
        beta_coeff.mc[k.mc, ] <- as.numeric(rmvn(1, beta_coeff_init,
                                                 (diagbetacov0)))
      }


      # Draw x-clusters within each y-cluster
      for (j in yclusters.mc) {
        # Number of observations (out of M total) in each y-cluster
        Nj.mc[j] <- sum(s.mc[, 1] == j)
        if (j != (k + 1)) {
          pxcluster <- rep(1, (maxslj[j] + 1))
          pxdenom <- nj[j] + alpha_omega
          for (l in 1:maxslj[j]) {
            pxcluster[l] <- nlj[j, l]/pxdenom
          }
          pxcluster[(maxslj[j] + 1)] <- alpha_omega/pxdenom
          s.mc[s.mc[, 1] == j, 2] <- sample(1:(maxslj[j] + 1), size = Nj.mc[j],
                                            replace = T, pxcluster)

          # If a new x-cluster was opened, draw parameter from prior
          if (max(s.mc[s.mc[, 1] == j, 2]) == (maxslj[j] + 1)) {
            for (rep in 1:(ptx + p1)) {
              x_prob_param.mat[j, (maxslj[j] + 1), rep] <- rbeta(n = 1,
                                                                 shape1 = a0, shape2 = b0)
            }

            if (p2 > 0) {
              for (rep in 1:p2) {
                x_var_param.mat[j, (maxslj[j] + 1), rep] <- rinvchisq(n = 1,
                                                                      df = nu0, scale = tau0)
                x_mean_param.mat[j, (maxslj[j] + 1), rep] <- rnorm(n = 1,
                                                                   mean = mu0, sd = sqrt(1/c0 * x_var_param.mat[j, (maxslj[j] +
                                                                                                                      1), rep]))
              }
            }
          }
        }
      }

      # Draw covariates for each of the M observations-------------------
      # Calculate largest number of x-clusters
      max.num.xclusters.mc <- max(s.mc[, 2])
      # Initialize matrix
      Nlj.mc <- matrix(0, nrow = k.mc, ncol = max.num.xclusters.mc)

      for (j in yclusters.mc) {
        xcluster.mc <- sort(unique(s.mc[s.mc[, 1] == j, 2]))
        for (l in xcluster.mc) {
          mask <- s.mc[, 1] == j & s.mc[, 2] == l
          Nlj.mc[j, l] <- sum(mask)
          if (p1 > 0) {
            for (rep in 1:p1) {
              L.mc[mask, rep] <- rbinom(Nlj.mc[j, l], 1, x_prob_param.mat[j, l, (ptx + rep)])
            }
          }
          if (p2 > 0) {
            for (rep in (p1 + 1):(p1 + p2)) {
              L.mc[mask, rep] <- rnorm(n = Nlj.mc[j, l],
                                       mean = x_mean_param.mat[j, l, (rep - p1)],
                                       sd = sqrt(x_var_param.mat[j, l, (rep - p1)]))
            }
          }
        }
      }  #End draw of covariates-----------------------------------

      # Compute E(Y|A=0, L=l,theta,omega,s) and E(Y|A=0, L=l,theta,omega,s) using
      # M new draws for covariates.-----------------------------------------

      # Set treatment to 0 and 1, and add column of 1's.
      newXa1.mc <- cbind(1, 1, L.mc)  #set treatment to 1
      newXa0.mc <- cbind(1, 0, L.mc)  #set treatment to 0

      X.mc <- cbind(1, L.mc)

      # Calculate h0ia1/0 and h0ya1/0 for use in calculating causal effect.
      # Average covariate distribution over prior for each x_i.

      # Binary covariates (not including treatment)
      for (i in 1:M) {
        if (p1 > 0) {
          for (j in 1:p1) {
            beta_draw <- beta(a0 + X.mc[i, (ptx + j)], b0 - X.mc[i, (ptx + j)] + 1)
            h0xa1.mc[i, (ptx + j)] <- beta_draw
            h0xa0.mc[i, (ptx + j)] <- beta_draw
          }
        }
        # Continuous covariates
        if (p2 > 0) {
          for (j in 1:p2) {
            x.mc <- X.mc[i, (p1 + ptx + j)]
            taun.mc <- (1/nun) * (nu0 * tau0 + (c0/(c0 + 1)) * (mu0 - x.mc)^2)
            num.mc <- (sqrt(2 * pi)/sqrt(cn)) * gamma(nun/2) * (2/(nun *
                                                                     taun.mc))^(nun/2)
            frac.mc <- num.mc/(denom * sqrt(2 * pi))  #sqrt(2*pi) is from the data part of likelihood
            h0xa1.mc[i, (p1 + ptx + j)] <- frac.mc
            h0xa0.mc[i, (p1 + ptx + j)] <- frac.mc
          }
        }
      }

      # Take product (covariates are assumed to be locally independent).
      # Result is vector of size M
      h0ia1.mc <- apply(h0xa1.mc, 1, prod)
      h0ia0.mc <- apply(h0xa0.mc, 1, prod)

      # Find E0(Y|A,L,beta_coeff) integrating over P_0beta.
      # beta_coeff_prior.mc was calculated before Gibbs sampling loop
      for (i in 1:MCnum) {
        h0ya1.mc <- h0ya1.mc + plogis(cbind(1, 1, X.mc[, ((ptx + 1):(ptx +
                                                                       p1 + p2))]) %*% beta_coeff_prior.mc[i, ])
        h0ya0.mc <- h0ya0.mc + plogis(cbind(1, 0, X.mc[, ((ptx + 1):(ptx +
                                                                       p1 + p2))]) %*% beta_coeff_prior.mc[i, ])
      }
      h0ya1.mc <- h0ya1.mc/MCnum
      h0ya0.mc <- h0ya0.mc/MCnum

      # Calculate causal effects --------------------------------------
      # Calculate E[Y^a] using MC integration to integrate E(Y|A=a,L=l,thetastar,omegastar,s)
      # over covariate distribution

      # Initialize
      parta1.mc <- 0
      parta0.mc <- 0
      denompart1a1.mc <- 0
      denompart1a0.mc <- 0

      for (j in yclusters.mc) {
        # (Use yclusters.mc rather than 1:k.mc in case some y-clusters are empty)
        sum2a1.mc <- 0
        sum2a0.mc <- 0

        xcluster.mc <- sort(unique(s.mc[s.mc[, 1] == j, 2]))
        # Calculate K(a,l;omega*_l|j) for binary covariates
        for (l in xcluster.mc) {
          prodxa1.mc <- 1
          prodxa0.mc <- 1
          for (rep in 1:(ptx + p1)) {
            prodxa1.mc <- prodxa1.mc * dbinom(newXa1.mc[, (rep + 1)],
                                              1, x_prob_param.mat[j, l, rep])
            prodxa0.mc <- prodxa0.mc * dbinom(newXa0.mc[, (rep + 1)],
                                              1, x_prob_param.mat[j, l, rep])
          }

          # Calculate K(a,l;omega*_l|j) for continuous covariates
          prodx2.mc <- 1
          if (p2 > 0) {
            for (rep in 1:p2) {
              prodx2.mc <- prodx2.mc * dnorm(L.mc[, (p1 + rep)], x_mean_param.mat[j, l, rep],
                                             sqrt(x_var_param.mat[j, l, rep]))
            }
          }
          # Calculate sum^{k_j}_{l=1}\frac{n_{l|j}}{alpha_omega+n_j}*K(a,l;omega*_l|j)
          #(This is part of the denominator of E(Y|A=a,L=l,thetastar,omegastar, s))
          sum2a1.mc <- sum2a1.mc + (Nlj.mc[j, l]/(alpha_omega + Nj.mc[j])) *
            prodxa1.mc * prodx2.mc
          sum2a0.mc <- sum2a0.mc + (Nlj.mc[j, l]/(alpha_omega + Nj.mc[j])) *
            prodxa0.mc * prodx2.mc

        }
        # Calculate sum^k_(j=1) w_j(a,l)E(y|a,l,theta_j*)
        # (This is part of the numerator of E(Y|A=a,L=l,thetastar,omegastar, s))
        parta1.mc <- parta1.mc + (Nj.mc[j]/(M + alpha_theta)) * ((alpha_omega/(alpha_omega +
                                                                                 Nj.mc[j])) * h0ia1.mc + sum2a1.mc) * plogis(newXa1.mc %*%  beta_coeff.mc[j, ])
        parta0.mc <- parta0.mc + (Nj.mc[j]/(M + alpha_theta)) * ((alpha_omega/(alpha_omega +
                                                                                 Nj.mc[j])) * h0ia0.mc + sum2a0.mc) * plogis(newXa0.mc %*% beta_coeff.mc[j, ])

        # Calculate sum^k_(j=1) w_j(a,l) (This is part of the denominator of
        # E(Y|A=a,L=l,thetastar,omegastar, s))
        denompart1a1.mc <- denompart1a1.mc + (Nj.mc[j]/(M + alpha_theta)) *
          ((alpha_omega/(alpha_omega + Nj.mc[j])) * h0ia1.mc + sum2a1.mc)
        denompart1a0.mc <- denompart1a0.mc + (Nj.mc[j]/(M + alpha_theta)) *
          ((alpha_omega/(alpha_omega + Nj.mc[j])) * h0ia0.mc + sum2a0.mc)
      }

      # Numerator of E(Y|A=a,L=l,thetastar,omegastar, s)
      numa1.mc <- parta1.mc + (alpha_theta/(alpha_theta + M)) * h0ia1.mc * h0ya1.mc
      numa0.mc <- parta0.mc + (alpha_theta/(alpha_theta + M)) * h0ia0.mc *  h0ya0.mc

      # Denominator of E(Y|A=a,L=l,thetastar,omegastar, s)
      denoma1.mc <- denompart1a1.mc + (alpha_theta/(alpha_theta + M)) *  h0ia1.mc
      denoma0.mc <- denompart1a0.mc + (alpha_theta/(alpha_theta + M)) *  h0ia0.mc

      # Calculate average causal effects Average causal relative risk: E[Y^1]/E[Y^0]
      psi_rr[gibbsreps] <- mean(numa1.mc/denoma1.mc)/mean(numa0.mc/denoma0.mc)

      # Average causal risk difference: E[Y^1] - E[Y^0]
      psi_rd[gibbsreps] <- mean(numa1.mc/denoma1.mc) - mean(numa0.mc/denoma0.mc)

      # Error handling
      if (is.finite(psi_rr[gibbsreps]) != TRUE) {
        stop()
      }

    } ## MC Int loop

    print(gibbsreps)

  } ## Gibbs loop

  # Calculate summary statistics from posterior distributions

  # Calculate median of posterior for concentration parameters
  at <- alpha_theta_draws[burn_in:gibbs_total]
  ao <- alpha_omega_draws[burn_in:gibbs_total]

  median_alpha_theta <- median(alpha_theta_draws[burn_in:gibbs_total])
  median_alpha_omega <- median(alpha_omega_draws[burn_in:gibbs_total])

  sd_alpha_theta <- sd(at)
  sd_alpha_omega <- sd(ao)
  lower_alpha_theta <- quantile(at, 0.025)
  upper_alpha_theta <- quantile(at, 0.975)
  lower_alpha_omega <- quantile(ao, 0.025)
  upper_alpha_omega <- quantile(ao, 0.975)

  # Calculate median of posterior for average relative risk
  median_psi_rr <- median(psi_rr[psi_rr != 0])
  sd_psi_rr <- sd(psi_rr[psi_rr != 0])

  # Calculate width of credible interval
  width_psi_rr <- abs(min(quantile(psi_rr[psi_rr != 0], c(0.025, 0.975))) - max(quantile(psi_rr[psi_rr != 0],
                                                                                         c(0.025, 0.975))))

  lower_psi_rr <- quantile(psi_rr[psi_rr != 0], 0.025)
  upper_psi_rr <- quantile(psi_rr[psi_rr != 0], 0.975)

  # Calculate median of posterior for average risk difference
  median_psi_rd <- median(psi_rd[psi_rd != 0])
  sd_psi_rd <- sd(psi_rd[psi_rd != 0])

  # Calculate width of credible interval
  width_psi_rd <- abs(min(quantile(psi_rd[psi_rd != 0], c(0.025, 0.975))) - max(quantile(psi_rd[psi_rd != 0], c(0.025, 0.975))))

  # Save results to a vector
  # results <- cbind(median_alpha_theta, median_alpha_omega, median_psi_rr, sd_psi_rr,
  #                  cov_psi_rr, width_psi_rr, median_psi_rd, sd_psi_rd, cov_psi_rd, width_psi_rd)
  results <- data.frame(
    Name = c("alpha_theta", "alpha_omega", "psi_rr"),
    Median = c(median_alpha_theta, median_alpha_omega, median_psi_rr),
    SE = c(sd_alpha_theta, sd_alpha_omega, sd_psi_rr),
    Lower = c(lower_alpha_theta, lower_alpha_omega, lower_psi_rr),
    Upper = c(upper_alpha_theta, upper_alpha_omega, upper_psi_rr)
  )

  ## Return result

  return(list(summary = results,
              alpha_theta_draws = at,
              alpha_omega_draws = ao,
              psi_rr_draws = psi_rr[psi_rr != 0]))

}
