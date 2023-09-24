## Load Dataset ----

library(CausalBNPBook)
library(tidyverse)
library(splines)
library(BART)
library(rpart)
library(rpart.plot)

set.seed(843592734)

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

## Fit BART by stratifying ----

meps_train_0 <- meps %>% filter(smoke == 0) %>% select(-y) %>% as.data.frame()
meps_train_1 <- meps %>% filter(smoke == 1) %>% select(-y) %>% as.data.frame()
bart_fit_0 <- wbart(meps_train_0,
                    meps$y[meps$smoke == 0],
                    x.test = as.data.frame(meps_test_0))
bart_fit_1 <- wbart(meps_train_1,
                    meps$y[meps$smoke == 1],
                    x.test = as.data.frame(meps_test_1))

## Exact Effects ----

pred_diff <- bart_fit_1$yhat.test.mean - bart_fit_0$yhat.test.mean
meps_test_ns <- meps_test_0 %>% select(-smoke) %>% mutate(delta = pred_diff)
rpart_diff <- rpart(delta ~ ., data = meps_test_ns)

rpart.plot(rpart_diff)
