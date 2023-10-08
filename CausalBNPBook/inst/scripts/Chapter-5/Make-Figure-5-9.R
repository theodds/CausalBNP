## Load Dataset ----

library(CausalBNPBook)
library(tidyverse)
library(splines)
library(BART)

## Preprocess meps ----

data(meps)

meps <- meps %>%
  as_tibble %>%
  mutate(y = log(y), sex = factor(sex), race = factor(race),
         marriage = factor(marriage), seat_belt = factor(seat_belt),
         edu = factor(edu))

## Create testing set ----

meps_test_0 <- meps %>% mutate(smoke = 0) %>% select(-y)
meps_test_1 <- meps %>% mutate(smoke = 1) %>% select(-y)
meps_test_01 <- rbind(meps_test_0, meps_test_1)

## Fit BART using smoking as a covariate ----

bart_with_covariates <- wbart(x.train = meps %>% select(-y) %>% as.data.frame(),
                              y.train = meps$y,
                              x.test = meps_test_01 %>% as.data.frame())


## Get Bayesian Bootstrap Weights ----

bb_weight <- MCMCpack::rdirichlet(n = nrow(bart_with_covariates$yhat.test),
                                  alpha = rep(1,nrow(meps)))

mean_0 <- rowSums(bb_weight * bart_with_covariates$yhat.test[,1:nrow(meps)])
mean_1 <- rowSums(bb_weight * bart_with_covariates$yhat.test[,-(1:nrow(meps))])

## Fit BART by stratifying ----

meps_train_0 <- meps %>% filter(smoke == 0) %>% select(-y) %>% as.data.frame()
meps_train_1 <- meps %>% filter(smoke == 1) %>% select(-y) %>% as.data.frame()
bart_fit_0 <- wbart(meps_train_0,
                    meps$y[meps$smoke == 0],
                    x.test = as.data.frame(meps_test_0))
bart_fit_1 <- wbart(meps_train_1,
                    meps$y[meps$smoke == 1],
                    x.test = as.data.frame(meps_test_1))

## Get estimates ----

mean_0_strat <- rowSums(bb_weight * bart_fit_0$yhat.test)
mean_1_strat <- rowSums(bb_weight * bart_fit_1$yhat.test)

## Plot effect estimates ----

effect_est <- mean_1 - mean_0
effect_est_strat <- mean_1_strat - mean_0_strat

df_to_plot <- tibble(
  effect = c(effect_est, effect_est_strat),
  method = rep(c("Covariate", "Stratified"), each = length(effect_est))
)

ggplot(df_to_plot, aes(x = effect)) +
  geom_histogram(color = 'white') +
  facet_wrap(~method) +
  xlab("$\\Delta$") +
  ylab("Frequency") +
  theme_bw()

print(mean(effect_est > 0))
print(mean(effect_est_strat > 0))



