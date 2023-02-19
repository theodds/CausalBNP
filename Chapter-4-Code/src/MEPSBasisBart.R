## Load ----

library(tikzDevice)
library(tidyverse)
library(splines)

meps <- read.csv(file = "data/meps.csv") %>% as_tibble %>% mutate(y = log(y))

## Fit splines ----

fitted_splines <- lm(y ~ ns(age, df = 5), data = meps)
fitted_linear <- lm(y ~ age, data = meps)

pchisq(deviance(fitted_linear) - deviance(fitted_splines), 
       df = df.residual(fitted_linear) - df.residual(fitted_splines),
       lower.tail=FALSE)

## Get prediction bands ----

new_meps <- tibble(age = seq(from = 20, to = 100, length = 500))
predicted_splines <- predict(fitted_splines, newdata = new_meps, se.fit = TRUE)
predicted_lm <- predict(fitted_linear, newdata = new_meps, se.fit = TRUE)

plot_splines <- function() {
  pal <- viridis::viridis(10)
  plot(new_meps$age, predicted_splines$fit, col = 3, type = 'l', lwd = 1.5, 
       xlab = "Age", ylab = "Log-Income", ylim = range(predicted_lm$fit))
  lines(new_meps$age, with(predicted_splines, fit - 2 * se.fit), lty = 2, 
        col = 8, lwd = 1.9)
  lines(new_meps$age, with(predicted_splines, fit + 2 * se.fit), lty = 2,
        col = 8, lwd = 1.9)
}

plot_splines()

## Fit LM ---- 

plot_line <- function() {
  pal <- viridis::viridis(10)
  plot(new_meps$age, predicted_lm$fit, col = 3, type = 'l', lwd = 1.5, 
       xlab = "Age", ylab = "Log-Income")
  lines(new_meps$age, with(predicted_lm, fit - 2 * se.fit), lty = 2, 
        col = 8, lwd = 1.9)
  lines(new_meps$age, with(predicted_lm, fit + 2 * se.fit), lty = 2,
        col = 8, lwd = 1.9)
}

plot_line()

## Do the same for BART ----
# 
# fitted_bart <- BART::wbart(x.train = as.matrix(as.numeric(meps$age)),
#                             y.train = meps$y, 
#                             x.test = as.matrix(new_meps$age), 
#                            ndpost = 5000, nskip = 5000)
# 
# ## Plot the BART ----
# 
# quant_bart <- apply(X = fitted_bart$yhat.test, MARGIN = 2, 
#                     FUN = function(x) quantile(x,c(.025,.975)))
# 
# preds_bart <- tibble(mean = fitted_bart$yhat.test.mean, 
#                      LCL = quant_bart[1,],
#                      UCL = quant_bart[2,],
#                      age = new_meps$age)
# 
# plot_bart <- function() {
#   pal <- viridis::viridis(10)
#   with(preds_bart, {
#     plot(age, mean, col = 3, type = 'l', lwd = 1.5, 
#          xlab = "Age", ylab = "Log-Income")
#     lines(age, LCL, lty = 2, col = 8, lwd = 1.9)
#     lines(age, UCL, lty = 2, col = 8, lwd = 1.9)
#   })
# }
# 
# plot_bart()

## Plot Both ----

tikz("SplineLine.tex", standAlone = TRUE, width = 5*1.2, height = 2.5*1.2)

par(mfrow = c(1,2), mar = c(5,4,1,1))

plot_splines()
plot_line()
# plot_bart()

dev.off()
tools::texi2pdf("SplineLine.tex", clean = TRUE)


