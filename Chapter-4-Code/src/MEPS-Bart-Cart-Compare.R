## Load stuff ----

library(tikzDevice)
library(tidyverse)
library(splines)
library(BART)
library(rpart)

set.seed(843592734)

meps <- read.csv(file = "data/meps.csv") %>% as_tibble %>% mutate(y = log(y))
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

tikz("BartCartCompare.tex", width = 5, height = 3, standAlone =  TRUE)
gridExtra::grid.arrange(p_1, p_2, nrow = 1)
dev.off()
tools::texi2pdf("BartCartCompare.tex", clean = TRUE)

## Short CV Experiment ----

# folds <- caret::createFolds(meps_age$y)
# 
# eval_folds <- function(c) {
#   n_folds <- length(folds)
#   out <- numeric(length(folds))
#   for(i in 1:n_folds) {
#     train_set <- meps_age[-folds[[i]],]
#     test_set <- meps_age[folds[[i]],]
#     fitted <- rpart(y ~ age, data = train_set, cp = c)
#     preds <- predict(fitted, test_set)
#     out[i] <- sum((test_set$y - preds)^2)
#   }
#   sqrt(sum(out) / nrow(meps_age))
# }
# 
# eval_folds(.01)
