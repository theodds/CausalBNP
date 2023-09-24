# CausalBNP

This repository contains code for replicating some of the methods described by
[Daniels, Linero, and
Roy](https://www.routledge.com/Bayesian-Nonparametrics-for-Causal-Inference-and-Missing-Data/Daniels-Linero-Roy/p/book/9780367341008).
It is structured as follows:

-   Each chapter for which there is code has a directory. Inside this directory
    is a markdown document describing how to reproduce certain methods used in
    each chapter.
    
-   This repository also contains an R package __CausalBNPBook__ that contains
    routines for fitting the models. It can be installed either by (i)
    downloading the repository and installing the package from source (e.g., by
    opening the .Rproj object in RStudio and using the "Build" pane to install
    it) or by (ii) installing directly from the pepository using the
    `install_github` function in the __devtools__ package (
    `devtools::install_github("theodds/CausalBNP", subdir = "CausalBNPBook")`).
    
-   R Markdown files for each of the sets of code is available in the
    __CausalBNPBook__ package as well, in the `inst/scripts/` folder. Knitting
    this code should "just work" as a way to reproduce the analyses, provided
    that the necessary packages are installed.

-   Prior to installing __CausalBNPBook__ it is recommended to install the
    dependencies. Unfortunately, the package __DPpackage__ is required to
    reproduce the examples from Chapter 13 (which uses the __BNPMediation__
    package, which is also only available on Github), and this package is no
    longer available on CRAN. An old version can be installed from source by
    running the following code:
  
    ```{r, eval = FALSE}
    str_dppack <- stringr::str_c("https://cran.r-project.org/src/contrib/",
                                 "Archive/DPpackage/DPpackage_1.1-6.tar.gz")
    devtools::install_url(str_dppack)
    devtools::install_github("lit777/BNPMediation")
    ```

# Status of the Page

The software given here uses a wide variety of packages for model fitting, as it
is comprised of code from different papers over many years. For example, some
Dirichlet process methods are fit using __rjags__, some are fit using
__DPpackage__, and some (e.g., the EDP model) are fit using custom C++ code. A
more unified package is being developed, and will be made available here when it
is complete. The main function of the __CausalBNPBook__ is to streamline the use
of this preexisting code, adding documentation and examples of usage for
estimating causal effects.

This page will be updated over time to include more (and better) software and
vignettes as they are created. At the moment, software implementing the
following methods is available:

-   The DP-BART approach to causal quantile estimation (Chapter 8).
-   An enriched Dirichlet process mixture model for causal inference with missing
    covariates, with binary outcomes (Chapter 9).
-   Dirichlet process observed data models for continuous outcomes subject to
    missing data (Chapter 11).
-   Dirichlet process observed data models for binary outcomes subject to
    non-monotone missingness (Chapter 12).
-   Dirichlet process mixture models for mediation analysis with continuous
    outcomes (Chapter 13).
-   Bayesian additive regression tree models for mediation analysis, assuming both
    a continuous mediator and a continuous outcome (Chapter 14).
-   A DDP-GP model for semicompeting risks (Chapter 15).

__Currently missing models used in the book are:__

- Dependent Dirichlet process mixtures of Gaussian process (DDP-GP) models for
  causal inference (Chapter 10).
- The enriched Dirichlet process mixture model for causal inference with missing
  covariates, with continuous outcomes (Chapter 11).

# Future Plans

- A more thorough implementation of enriched Dirichlet process mixture models as
  a package (intended to be available on CRAN) is currently in development, and
  will be added to this repository when it is ready.
- A more thorough/cohesive package combining all of the models here is in
  development. 
- Examples showing how to fit DDP-GP models 
