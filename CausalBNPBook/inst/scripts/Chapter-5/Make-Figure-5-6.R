## Load packages ----

library(CausalBNPBook)
library(tidyverse)
library(splines)
library(BART)
library(rpart)

set.seed(843592734)

## Preprocess the meps dataset ----

data(meps)

meps <- meps %>% as_tibble %>% mutate(y = log(y))
meps_age <- meps %>% select(age, y)
meps_age_test <- tibble(age = seq(from = min(meps$age),
                                  to = max(meps$age),
                                  length = 200))

## Fit CART ----

fitted_cart <- rpart(y ~ age, data = meps_age, cp = .001)
meps_age_test$cart_pred <- predict(fitted_cart, meps_age_test)

## Fit BART ----

fitted_bart <- wbart(x.train = select(meps_age, -y) %>% as.data.frame(),
                     y.train = meps_age$y,
                     x.test = meps_age_test %>% as.data.frame())

## Plot results ----

meps_age_test$bart_pred <- fitted_bart$yhat.test.mean

p_1 <- ggplot(meps_age_test, aes(x = age, y = cart_pred)) +
  geom_line() +
  xlab("Age") +
  ylab("CART Prediction") +
  theme_bw()

p_2 <- ggplot(meps_age_test, aes(x = age, y = bart_pred)) +
  geom_line() +
  xlab("Age") +
  ylab("BART Prediction") +
  theme_bw()

gridExtra::grid.arrange(p_1, p_2, nrow = 1)

