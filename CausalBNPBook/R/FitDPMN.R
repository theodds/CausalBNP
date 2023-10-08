#' Fit Dirichlet Process Mixture of Multivariate Normals Model
#'
#' This function fits a Dirichlet process mixture of multivariate normals outcome
#' model to longitudinal outcome data subject to informative missingness. The model
#' is designed to account for continuous pretreatment predictors by prepending them
#' to the vector of outcome variables. Discrete variables can be handled through
#' stratification.
#'
#' @details
#' The observed data model being fitted has the following form:
#'
#' \deqn{
#'  \begin{array}{rl}
#'   Z_i &\sim \mbox{Categorical}(\omega), \\
#'   R_{i1} &\equiv 1, \\
#'   [Y_{i1} \mid Z_i] &\sim N(\mu^{(Z_i)}_1, \sigma^{2(Z_i)}_1), \\
#'   \mbox{logit}\Pr(R_{i(j+1)} =  1 \mid R_{ij} = 1, \overline{Y}_{i(j-1)}, Z_i)
#'     &= \zeta^{(Z_i)}_j + \gamma^{(Z_i)}_{j1} \,  Y_{ij} + \,
#'         + \gamma_{j2}^{(Z_i)} Y_{i(j-1)}, \\
#'   [Y_{ij} \mid R_{ij} = 1, \overline{Y}_{i(j-1)}, Z_i] &\sim
#'     N\left\{\mu^{(Z_i)}_j +
#'         \sum_{\ell = 1}^{j-1} \phi^{(Z_i)}_{j\ell}(Y_{i\ell} - \mu_\ell^{(Z_i)}),
#'         \sigma^2\right\}
#'  \end{array}
#' }
#'
#' This model does not account for any predictors. However, it's possible to
#' account for continuous pretreatment predictors by prepending them to the
#' vector \eqn{Y_i = (Y_{i1}, \ldots, Y_{iJ})} of outcome variables. Discrete
#' variables can be allowed through stratification as well.

#' @param Y A matrix of outcome measurements. Each row corresponds to a set of measurements
#'          for a given subject, and each column corresponds to a specific measurement time.
#' @param K An integer indicating the number of mixture components in the truncated Dirichlet
#'          process mixture. Defaults to 20.
#' @param n.adapt An integer specifying the number of adaptation iterations in JAGS. Defaults to 1000.
#' @param n.save An integer specifying the number of iterations to save. Defaults to 5000.
#' @param thin An integer specifying the thinning interval for the MCMC samples. Defaults to 1.
#' @param n.chains An integer specifying the number of MCMC chains. Defaults to 1.
#'
#' @return A list object of class runjags from the package runjags. This object contains
#'         MCMC samples and other relevant information related to the fit. Users can refer to the
#'         [runjags package documentation](https://cran.r-project.org/web/packages/runjags/vignettes/quickjags.html)
#'         for further details on working with the returned object.
#'
#' @examples
#' \dontrun{
#'   library(CausalBNPBook)
#'   library(tidyverse)
#'   data(GH)
#'   Y1 <- GH %>% filter(trt == 1) %>% select(Y1:Y3) %>% as.matrix()
#'   fit_dpmn_1 <- FitDPMN(Y1, K = 20)
#' }
#'
#' @seealso [runjags package documentation](https://cran.r-project.org/web/packages/runjags/vignettes/quickjags.html)
#' @author Antonio R. Linero
#' @export
FitDPMN <- function(Y, K = NULL,
                    n.adapt = 1000, n.save = 5000, thin = 1, n.chains = 1) {


  R <- !is.na(Y) + 0

  ## Fits the nonparametric model that does not share across
  ## treatment groups. Works by stratifying; we do not pass in
  ## treatment indicators.

  get.R.1m <- function(R) {
    stopifnot(is.matrix(R))
    J <- dim(R)[2]
    R.1m <- 1 - R
    for(j in J:2) R.1m[,j] <- ifelse(R[,j-1] == 1, R.1m[,j], NA)
    return(R.1m)
  }

  load.module("glm")
  load.module("dic")
  J <- dim(Y)[2]
  N <- dim(Y)[1]
  n.iter <- n.save * thin
  if(is.null(K)) {
    K <- 20
    warning("\nUsing Default K = 20\n")
  }
  Y[R == 0] <- NA
  R.1m <- get.R.1m(R)
  datum <- list(
    K = K,
    Y = Y,
    J = J, N = N,
    ## R_1M = R.1m,
    R  = R,
    gtt = get.c.vars(Y)
  )


  collect <- c("Y.new", "R.new", "mass", "a", "phi",
               "tau", "xi", "zeta", "gamma", "C")

  out <- run.jags(model = noshare_str, monitor = collect, data = datum,
                  adapt = n.adapt, burnin = n.adapt, sample = n.save, n.chains = n.chains)

  return(out)
}

get.c.vars <- function(Y) {

  J <- dim(Y)[2]
  ## Get the parameters used in the GARP prior

  s <- prelim.norm(Y)
  e <- em.norm(s)
  params <- getparam.norm(s,e)

  G <- params$sigma
  gtt <- diag(G)
  gtt.tm1 <- gtt
  for(j in 2:J) {
    gtt.tm1[j] <- gtt[j] - G[j,1:(j-1),drop=FALSE] %*%
      solve(G[1:(j-1),1:(j-1)]) %*% G[1:(j-1),j,drop=FALSE]
  }
  return(gtt.tm1)
}


noshare_str <-
  "

model {

## Hyperparameters
tau.mass <- 1 / 4
tau.mu.j <- 1 / 5 / 5
tau.mu.phi <- 1 / 2.5 / 2.5


mass <- 1
for(k in 1:K) {
mass.vec[k] <- mass / K
}
xi ~ ddirich(mass.vec)


# Prior for a
for(j in 1:J) {
  for(k in 1:K) {
    a[k,j] ~ dnorm(mu[j], tau.a[j])
  }
}


## Prior on phi[k,j,l]
for(k in 1:K) {
  for(j in 1:J) {
    for(l in 1:J) {
      mm[k,j,l] <- ifelse(l==(j-1), mu.phi, 0)
      phi[k,j,l] ~ dnorm(mm[k,j,l], tau.phi[j])
    }
  }
}




lambda ~ dgamma(4,4)
shape ~ dgamma(5,1)
for(j in 1:J) {
  muu[j] <- lambda / gtt[j] * 2
  rate[j] <- shape / muu[j]
}
for(k in 1:K) {
  tau[k,1] ~ dgamma(shape, rate[1])
  for(j in 2:J) {
    tau[k,j] ~ dgamma(shape + 1, rate[j])
  }
}

## Prior for mu[j]
for(j in 1:J) {
  mu[j] ~ dnorm(0, .1)
}
## Prior for sigma.a[j] and tau.a[j]
shape.a ~ dgamma(J,1)
lambda.a ~ dgamma(4, 4)
for(j in 1:J) {
  sigma.a.raw[j] ~ dgamma(shape.a,1)
  sigma.a[j] <- sigma.a.raw[j] * sqrt(gtt[j] / 2) / shape.a * lambda.a
  tau.a[j] <- 1 / sigma.a[j] / sigma.a[j]
}
## Prior for mu.phi
mu.phi ~ dt(0, tau.mu.phi, 1)
## Prior for tau.phi
scale.phi ~ dgamma(1,2 / gtt[1])
shape.phi ~ dgamma(J, 1)
for(j in 1:J) {
sigma.phi.raw[j] ~ dgamma(shape.phi, 1)
sigma.phi[j] <- sigma.phi.raw[j] * scale.phi / shape.phi
tau.phi[j] <- 1/pow(sigma.phi[j],2)
}
## Prior for psi
psi ~ dt(0, tau.mu.j, 1)
## Prior for omega.inv
omega.inv ~ dgamma(0.1, 0.1)

## Prior on zeta and gamma
for(k in 1:K) {
for(j in 1:J) {
gamma[k,j] ~ dnorm(mu.gamma, tau.gamma)
zeta[k,j] ~ dnorm(mu.zeta, tau.zeta)
}
}
## And associated hypers
mu.zeta ~ dt(0, tau.mu.j, 1)
tau.zeta ~ dgamma(0.1, 0.1)
mu.gamma ~ dt(0, tau.mu.phi, 1)
tau.gamma ~ dgamma(0.1, 0.1)

## Finally, model the data
for(n in 1:N) {
C[n] ~ dcat(xi[])
Y[n,1] ~ dnorm(a[C[n],1], tau[C[n],1])
for(j in 2:J) {
m[n,j] <- a[C[n], j] + inprod(phi[C[n],j,1:(j-1)], Y[n,1:(j-1)])
Y[n,j] ~ dnorm(m[n,j], tau[C[n],j])
}
## LEAVE OFF WITH DROPOUT SITUATION
for(j in 1:(J-1)) {
prob[n,j] <- (1 - ilogit(zeta[C[n], j]
+ gamma[C[n],j] * Y[n,j])) *
R[n,j]
R[n,j+1] ~ dbern(prob[n,j])
}

}

## Predict New
C.new ~ dcat(xi[])
tmp <- 1
R.new[1] ~ dbern(tmp)
Y.new[1] ~ dnorm(a[C.new,1], tau[C.new, 1])
for(j in 2:J) {
m.new[j] <- a[C.new, j] + inprod(phi[C.new,j,1:(j-1)], Y.new[1:(j-1)])
Y.new[j] ~ dnorm(m.new[j], tau[C.new, j])
}
for(j in 1:(J-1)) {
prob.new[j] <- (1 - ilogit(zeta[C.new, j] +
gamma[C.new,j] * Y.new[j])) * R.new[j]
R.new[j+1] ~ dbern(prob.new[j])
}
}


"

#' G-Computation for Dirichlet Process Mixture of Multivariate Normals Model
#'
#' This function implements the \eqn{g}-computation algorithm to generate samples
#' of the causal parameters based on a Dirichlet process mixture of multivariate
#' normals model fitted using the `FitDPMN` function. The most important component
#' is the matrix change_from_baseline, which contains the change in the outcome
#' between time 1 and time \eqn{J}.
#'
#' @param fit A list object returned from the `FitDPMN` function, typically of class
#'            `runjags` containing MCMC samples and other relevant information related to the fit.
#' @param Y A matrix of outcome measurements, similar to the one provided to the `FitDPMN` function.
#' @param K An integer indicating the number of mixture components in the truncated Dirichlet process mixture.
#' @param sens.param A matrix of sensitivity parameters of size \eqn{M \times J} where the
#'                   \eqn{ (m,j)^{\text{th}} } entry is a value of \eqn{ \xi_j } to use at iteration \eqn{ m }.
#' @param nsim An integer specifying the number of simulated observations for the \eqn{g}-computation. Defaults to 100 times the sample size.
#'
#' @return
#' A list object containing the following elements:
#' \itemize{
#'   \item \code{means}: The mean parameters \eqn{E(Y_{ij})} for each treatment and time.
#'   \item \code{var_of_mean}: The Monte Carlo variance from the \eqn{ g }-computation.
#'   \item \code{obs.means}: The observed data mean at each time, useful for model checking.
#'   \item \code{drop.probs}: The probability of dropout at each time, useful for model checking.
#'   \item \code{haz}: The hazard of dropout at each time, useful for model checking.
#'   \item \code{drop.probs}: The cumulative probability of dropout at each time, i.e., \eqn{\Pr(R_{ij} = 1, R_{i(j+1)} = 0)}, useful for model checking.
#'   \item \code{change_from_baseline}: The change from baseline \eqn{E(Y_{iJ} - Y_{i1})}.
#'   \item \code{change_from_baseline_var}: The Monte Carlo variance of the change from baseline.
#' }
#'
#' @details
#' This function uses the package's compiled code to calculate the MNAR means under the assumption of no sharing
#' of the baseline distribution. The function extracts relevant parameters from the MCMC samples and performs the
#' \eqn{ g }-computation to generate the desired samples of causal parameters.
#'
#' @seealso \code{\link{FitDPMN}}
#' @author Antonio R. Linero
#' @export
#'
#' @examples
#' \dontrun{
#'   library(CausalBNPBook)
#'   library(tidyverse)
#'   data(GH)
#'   Y1 <- GH %>% filter(trt == 1) %>% select(Y1:Y3) %>% as.matrix()
#'   fit_dpmn_1 <- FitDPMN(Y1, K = 20)
#'   gcomp_dpmn_1 <- GCompDPMN(fit_dpmn_1, Y = Y1, K = 20, sens.param = array(0, c(dim(fit_dpmn_1$mcmc[[1]])[1], ncol(Y1))))
#'
#'   ## Plot the change from baseline for the first treatment group:
#'   qplot(gcomp_dpmn_1$change_from_baseline[,1])
#' }
GCompDPMN <- function(fit, Y, K, sens.param, nsim = NULL) {

  mcmc <- fit$mcmc
  if(is.null(nsim)) nsim <- 100 * nrow(Y)

  noshare.impute <- function(a, sigma, phi,
                             zeta, gamma, xi, sens.param,
                             num.sim = nsim) {

    ## Wrapper for the compiled code to calculate the MNAR means when
    ## NOT sharing the baseline distribution

    J <- dim(a)[2]
    K <- dim(a)[1]
    num.iter <- dim(a)[3]

    if(missing(sens.param))
      sens.param <- array(0, c(num.iter, J))

    .Call("noshare_impute", a, sigma, phi, zeta, gamma, xi, num.sim,
          num.iter, J, K, sens.param)
  }

  get.means <- function(cs, num.save, sens.param, K, J) {

    ## USES THE PACKAGE to calculate the means. Does this by
    ## extracting the relevant paramters calling the compiled code.

    ## Returns matrix of mnar means

    a <- jags.extract.2(cs[[1]], "a", c(K,J), num.save)
    num.iter <- dim(a)[3]
    sigma <- 1 / sqrt(jags.extract.2(cs[[1]], "tau", c(K,J), num.save))
    phi <- jags.extract.2(cs[[1]], "phi", c(K,J,J), num.save)
    zeta <- jags.extract.2(cs[[1]], "zeta", c(K,J), num.save)
    gamma <- jags.extract.2(cs[[1]], "gamma", c(K,J), num.save)
    xi <- jags.extract.2(cs[[1]], "xi", K, num.save)
    if(missing(sens.param))
      sens.param <- array(0, c(num.iter, J))

    noshare.impute(a,sigma,phi,zeta,gamma,xi,sens.param)
  }

  jags.extract.2 <- function(js, str, d, n.iter) {

    ## Same as jags.extract but BETTER because it puts the iteration
    ## in the last component instead of the first, i.e. a[i,j,k]
    ## returns a[i,j] for iteration k

    n.dim <- length(d)
    if(n.dim==0) {
      return(js[,str])
    }
    out <- array(0, c(d, n.iter))
    if(n.dim==1) {
      for(j in 1:d) {
        str2 <- paste(str, "[", j, "]", sep = '')
        out[j,] <- js[,str2]
      }
    }
    if(n.dim==2) {
      for(j1 in 1:d[1]) {
        for(j2 in 1:d[2]) {
          str2 <- paste(str, "[", j1, ",", j2, "]", sep = '')
          out[j1,j2,] <- js[,str2]
        }
      }
    }
    if(n.dim==3) {
      for(j1 in 1:d[1]) {
        for(j2 in 1:d[2]) {
          for(j3 in 1:d[3]) {
            str2 <- paste(str, "[", j1, ",", j2, ",", j3, "]", sep = '')
            out[j1,j2,j3,] <- js[,str2]
          }
        }
      }
    }
    return(out)
  }

  cs <- mcmc
  num.save <- nrow(cs[[1]])
  J <- ncol(Y)

  return(get.means(cs, num.save, sens.param, K, J))


}

