sse <- function(x) sum(x^2)
RMSE <- function(x, y) sqrt(sse(x - y) / length(x))

GetLambda <- function(y, X, nu=3, q=.9) {
  fit      <- lm(y ~ X)
  sigma_sq <- summary(fit)$sigma^2
  sigma    <- sqrt(sigma_sq)
  
  out <- qgamma(1 - q, shape=nu/2, rate=nu/2)
  out <- out * sigma_sq
  
  test_it <- sqrt(1.0 / rgamma(10000, nu/2, rate=out*nu/2))
  conf <- sum(test_it < sigma)
  
  out_list <- list(out=out, test_it=test_it, sigma=sigma, shape = nu / 2, 
                   rate = out * nu / 2, conf=conf)
  
  return(out_list)
}

Fried <- function(X) {
  10 * sin(pi * X[, 1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
}

FriedSamp <- function(n, dim_x) {
  X <- matrix(runif(n * dim_x), n, dim_x)
  mu <- Fried(X)
  Y <- rnorm(n, mu)
  
  return(list(X=X, Y=Y))
}

Standardize <- function(y) {
  L <- min(y)
  U <- max(y)
  
  return(list(y=(y - L) / (U - L) - 0.5, L=L,U=U))
}

Unstand <- function(y, L, U) {
  out <- (y + .5) * (U - L) + L
  return(out)
}