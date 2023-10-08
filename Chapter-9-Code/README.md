---
title: "Enriched Dirichlet Processes"
author: "Antonio R. Linero"
date: "2023-09-13"
output: html_document
editor_options: 
  markdown: 
    wrap: 80
  chunk_output_type: console
---

We illustrate the use of the Enriched Dirichlet process (EDP) for causal
inference using data simulated under settings of Roy et al. The original Github
repository is available [at this
link.](https://github.com/jasonroy0/EDP_causal/tree/master).

First, we load the required packages:


```r
library(CausalBNPBook)
library(tidyverse)
library(miscTools)
```

The following function generates some simulated data for us to use; it
corresponds to Scenario 1 under MAR in the manuscript by [Roy et
al.](https://pubmed.ncbi.nlm.nih.gov/29579341/):


```r
generate_data <- function(n) {
  # number of binary covariates
  p1 <- 2
  # number of continuous covariates
  p2 <- 2
  
  L <- matrix(0, nrow = n, ncol = (p1 + p2))
  A <- rep(0, n)
  Y <- rep(0, n)
  
  # Generate covariates
  L[, 1] <- rbinom(n, 1, 0.2)
  L[, 2] <- rbinom(n, 1, 0.3 + 0.2 * L[, 1])
  
  L[, 3] <- rnorm(n, L[, 1] - L[, 2], 1)
  L[, 4] <- rnorm(n, 1 + 0.5 * L[, 1] + 0.2 * L[, 2] - 0.3 * L[, 3], 2)
  
  
  # Generate treatment A given L
  prob_a <- plogis(-0.4 + L[, 1] + L[, 2] + L[, 3] - 0.4 * L[, 4])
  A <- rbinom(n, 1, prob_a)
  
  # Generate outcome Y
  prob_y <- plogis(-0.5 + 0.78 * A - 0.5 * L[, 1] - 0.3 * L[, 2] + 0.5 * L[, 3] - 0.5 * L[, 4])
  Y <- rbinom(n, 1, prob_y)
  
  ## Generate MAR Missingness
  # Generate missing at random data in covariates
  L_mar <- matrix(0, nrow = n, ncol = p1 + p2)
  L_mar[, 1] <- rbinom(n, 1, plogis(-1.5 + L[, 2] + Y))
  L_mar[, 2] <- rbinom(n, 1, plogis(-1.5 + L[, 3] + A))
  L_mar[, 3] <- rbinom(n, 1, plogis(-1 - A + Y))
  L_mar[, 4] <- rbinom(n, 1, plogis(-0.4 - L[, 1] - L[, 2]))
  
  return(data.frame(Y = Y, A = A, L = L, L_mar = L_mar))
}
```

We use this to generate some simulated data as:

```r
my_data <- generate_data(250)

Y <- my_data$Y
A <- my_data$A %>% as.matrix()
L <- my_data %>% select(starts_with("L.")) %>% as.matrix()
R <- my_data %>% select(starts_with("L_")) %>% as.matrix()

L_binary <- L[,1:2]
L_continuous <- L[,3:4]
R_binary <- R[,1:2]
R_continuous <- R[,3:4]

num_burn <- 1000 ## Make larger in practice
num_thin <- 1
num_save <- 1000 ## Make larger in practice
num_mc <- 1000 ## ~10000 used by Roy et al.
thin_mc_int <- 2 ## 10 used by Roy et al.
M <- 1000 ## 1000 used by Roy et al.
```

The `edp_binomial` fuction fits an Enriched Dirichlet process with a discrete
(binary) outcome, as we have above. It requires specifying both the (potentially
missing) covariates `L_continuous` and `L_binary`, which are assumed to have
been *initialized* prior to feeding into `edp_binomial`, as well as `R_binary`
and `R_continuous`, which give the missing indicators. This is, of course, in
addition to the binary outcome `Y` and the binary treatment `A`.

In addition to these quantities, the user supplies the usual parameters of the
Markov chain (`num_burn`, `num_thin`, and `num_save`) and some parameters for
performing the $g$-computation. The `thin_mc_int` variable gives the thinning
interval on performing the $g$-computation. The `num_mc` variable gives the
number of Monte Carlo samples to take when performing the $g$-computation; this
is used in the process of integrating out the base measure. The quantity `M`
gives the number of randomly-selected covariates to integrate over when
performing $g$-computation.

The following code fits the model:


```r
edp_fit <-
  edp_binomial(
    Y = Y,
    A = A,
    L_binary = L_binary,
    L_continuous = L_continuous,
    R_binary = R_binary,
    R_continuous = R_continuous,
    num_burn = num_burn,
    num_thin = num_thin,
    num_save = num_save,
    num_mc = num_mc,
    thin_mc_int = thin_mc_int,
    M = M
  )
```

The following output gives a summary of the overall fit:

```r
print(edp_fit$summary)
```

```
##          Name    Median        SE      Lower     Upper
## 1 alpha_theta 0.5762908 0.3895720 0.13633709 1.6463300
## 2 alpha_omega 0.1044235 0.1364042 0.03784082 0.4612751
## 3      psi_rr 1.7002281 0.2580635 1.20438618 2.2063953
```

The true relative risk for this simulation setting is 1.5, which is comfortably
within the credible interval in this case.

