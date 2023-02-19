## Load ----

library(tikzDevice)
library(tidyverse)
library(splines)

# meps <- read.csv(file = "data/meps.csv") %>% as_tibble %>% mutate(y = log(y))
old_meps <- tibble(age = seq(from = 20, to = 100, length = 500))
new_meps <- tibble(age = seq(from = 20, to = 100, length = 500))

## Construct Basis ----

spline_basis <- ns(x = old_meps$age, df = 10)
new_meps_basis <- ns(x = new_meps$age, 
                     knots = attr(spline_basis, "knots"),
                     Boundary.knots = attr(spline_basis, "Boundary.knots"))

my_plot <- function() {
  xlim <- range(new_meps$age)
  ylim <- range(new_meps_basis)
  plot(xlim, ylim, type = 'n', xlab = "$x$", ylab = "$\\psi_b(x)$")
  for(j in 1:ncol(new_meps_basis)) {
    lines(new_meps$age, new_meps_basis[,j])
  }
}

tikz("Spline-Basis.tex", width = 5, height = 3, standAlone = TRUE)
par(mar = c(5,4,1,1))
my_plot()
dev.off()
tools::texi2pdf("Spline-Basis.tex", clean = TRUE)
