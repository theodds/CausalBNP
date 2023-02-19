## Load ----

library(GPBart)

set.seed(473920)

P <- 200
p <- rep(1/P,P)
alpha <- 0.95
beta <- 2
cov_select <- p
num_tree <- 50
sigma_mu <- 3/(2*sqrt(num_tree))
N <- 1000


## Simulate X ----

X <- matrix(runif(N * P), nrow = N)

## Sample Mu ----

sample_Phi <- function() {
  out <- replicate(num_tree, TreeCompress(alpha,beta,p,X))
  my_features <- do.call(cbind,out)
  beta <- rnorm(n = ncol(my_features), sd = sigma_mu)
  mu <- my_features %*% beta
}

mu_1 <- sample_Phi()
mu_0 <- sample_Phi()
e    <- pnorm(sample_Phi())
A    <- ifelse(runif(length(e)) < e, 1, 0)

cor_val <- cor(mu_1, e)


do_sim <- function() {
  X <- matrix(runif(N * P), nrow = N)
  mu_1 <- sample_Phi()
  mu_0 <- sample_Phi()
  e    <- pnorm(sample_Phi())
  A    <- ifelse(runif(length(e)) < e, 1, 0)
  cor_val <- cor(mu_1, e)
}

out <- replicate(1000, do_sim())

hist(out)
