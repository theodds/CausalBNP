# CausalBNP

This repository contains code for replicating some of the methods described by
[Daniels, Linero, and
Roy](https://www.routledge.com/Bayesian-Nonparametrics-for-Causal-Inference-and-Missing-Data/Daniels-Linero-Roy/p/book/9780367341008).
It is structured as follows:

- Each chapter for which there is code has a directory. Inside this directory is
  a markdown document describing how to reproduce certain methods used in each
  chapter.
  
- This repository also contains an R package __CausalBNPBook__ that contains
  routines for fitting the models. It can be installed either by (i) downloading
  the repository and installing the package from source (e.g., by opening the
  .Rproj object in RStudio and using the "Build" pane to install it) or by (ii)
  installing directly from the pepository using the `install_github` function in
  the __devtools__ package (
  `devtools::install_github("theodds/CausalBNP", subdir = "CausalBNPBook")`).
  
- 

# Status of the Page
