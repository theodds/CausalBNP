## Load required packages ----

library(CausalBNPBook)
library(tidyverse)

devtools_installed <- require(devtools, quietly = TRUE)
bnpqte_installed <- require(BNPqte, quietly = TRUE)

if(!devtools_installed) {
  install.packages("devtools")
}

if(!bnpqte_installed) {
  install.packages("plotly")
  devtools::install_github("https://github.com/chujiluo/BNPqte")
  require(BNPqte)
}

## Set up data ----

set.seed(20983)

meps <- CausalBNPBook::meps
idx <- sample(seq_len(nrow(meps)), 1000)

y <- log(meps$y + 1000)[idx]
x <- model.matrix(~age + sex + race + loginc + marriage + seat_belt + edu - 1,
                  data = meps[idx,])
a <- meps$smoke[idx]


## Getting the qte ----

meps_qte <- qte(y = y, x = x, treatment = a)

## Plot the result ----

results <- tibble(
  logexp = rep(meps_qte$grid, each = nrow(meps_qte$control.cdfs)),
  iteration = rep(seq_len(nrow(meps_qte$control.cdfs)), length(meps_qte$grid)),
  trt_cdf = as.numeric(meps_qte$treatment.cdfs),
  ctr_cdf = as.numeric(meps_qte$control.cdfs)
)

results <- results %>% mutate(delta = trt_cdf - ctr_cdf)

results_sum <- results %>%
  group_by(logexp) %>%
  summarise(delta_mean = mean(delta),
            delta_min = quantile(delta, 0.025),
            delta_max = quantile(delta, 0.975))

ggplot(results_sum, aes(x = logexp)) +
  geom_line(aes(y = delta_mean)) +
  geom_ribbon(aes(ymin = delta_min, ymax = delta_max), alpha = 0.2)

## Visualizing Posterior Differences in Quantiles ----

get_qt_df <- function(j) {
  tibble(
    variable = c("Treated", "Control", "Difference"),
    mean = c(meps_qte$treatment.quantiles.avg[j],
             meps_qte$control.quantiles.avg[j],
             meps_qte$qtes.avg[j]),
    lcl = c(meps_qte$treatment.quantiles.ci[j,1],
            meps_qte$control.quantiles.ci[j,1],
            meps_qte$qtes.ci[j,1]),
    ucl = c(meps_qte$treatment.quantiles.ci[j,2],
            meps_qte$control.quantiles.ci[j,2],
            meps_qte$qtes.ci[j,2]),
    quantile = names(meps_qte$qtes.avg)[j]
  )
}

qt_df <- do.call(rbind, lapply(1:5, get_qt_df))

qt_df %>% filter(variable == "Difference") %>%
  ggplot(aes(x = quantile)) +
  geom_point(aes(y = mean)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl)) +
  xlab("Quantile") +
  ylab("QTE") +
  theme_bw()



