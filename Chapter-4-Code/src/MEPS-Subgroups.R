## Load Dataset ----

library(tikzDevice)
library(tidyverse)
library(splines)
library(BART)

set.seed(843592734)

meps <- read.csv(file = "data/meps.csv") %>% as_tibble %>% mutate(y = log(y))

## Others ----

meps_subgroup <- meps %>% filter(
  age > 50, 
  sex == "male", 
  race == "white", 
  seat_belt == "always"
  )

x_smoke_1 <- meps_subgroup %>% filter(smoke == 1) %>% select(y) %>% unlist()
x_smoke_0 <- meps_subgroup %>% filter(smoke == 0) %>% select(y) %>% unlist()

length(x_smoke_1)
length(x_smoke_0)
t.test(x_smoke_1, x_smoke_0, paired = FALSE)

x_glob_1 <- meps %>% filter(smoke == 1) %>% select(y) %>% unlist()
x_glob_0 <- meps %>% filter(smoke == 0) %>% select(y) %>% unlist()



t.test(x_glob_1, x_glob_0, paired = FALSE)

## yougn ----

meps_subgroup_young <- meps %>% filter(
  age < 32, 
  sex == "male", 
  race == "white", 
  seat_belt == "always"
)

x_y_1 <- meps_subgroup_young %>% filter(smoke == 1) %>% select(y) %>% unlist()
x_y_0 <- meps_subgroup_young %>% filter(smoke == 0) %>% select(y) %>% unlist()

length(x_y_1)
length(x_y_0)
t.test(x_y_1, x_y_0,paired = FALSE)
