## Load Dataset ----

library(tikzDevice)
library(tidyverse)
library(splines)
library(BART)
library(xgboost)
library(randomForest)
library(gbm)

set.seed(843592734)

meps <- read.csv(file = "data/meps.csv") %>% 
  as_tibble %>% 
  mutate(y = log(y)) %>%
  select(y, age)

test_df <- tibble(age = seq(from = 20, to = 100, length = 500))

## Fit RF ----

fitted_rf <- randomForest(x = meps %>% select(age), y = meps$y)
pred_rf <- predict(fitted_rf, test_df)

## Fit boosting ----

fitted_boost <- gbm(y ~ age, data = meps, distribution = "gaussian")
pred_boost <- predict(fitted_boost, test_df, n.trees = 200)

## Plot results ----

results_df <- tibble(age = rep(test_df$age,2), 
                     fage = c(pred_rf, pred_boost), 
                     method = rep(c("Random Forests", "Boosting"), 
                                  each = nrow(test_df)))

tikz("BagBoost.tex", width = 5.5, height = 3*.8, standAlone =  TRUE)

ggplot(results_df, aes(x = age, y = fage)) + 
  geom_line() + 
  facet_wrap(~method) + 
  xlab("Age") +
  ylab("Prediction") +
  theme_bw()

dev.off()
tools::texi2pdf("BagBoost.tex", clean = TRUE)