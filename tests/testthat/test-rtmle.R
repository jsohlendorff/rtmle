library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
library(targets)


set.seed(112)
ld <- simulate_long_data(n = 10000,number_epochs = 20,beta = list(A_on_Y = -.2,A0_on_Y = -0.3,A0_on_A = 6),register_format = TRUE)
## w_outcome <- widen_outcome(outcome_name = "Y",pop = pop,intervals = seq(0,2000,30.45*6),como = list("Y" = outcome_data),id = "id")
## w_treatment <- widen_treatment(treatment = treatment_data,treatment_name = "A",pop = pop,intervals = seq(0,2000,30.45*6),id = "id")
## w_como <- widen_treatment(treatment_name = "L",treatment = como_data,pop = pop,intervals = seq(0,2000,30.45*6),id = "id")

# rtmle
## tar_source("~/tmp/rtmle/R/")

x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
add_data(x) <- ld
prepare_data(x,intervals = seq(0,2000,30.45*6)) <- "A"
protocol(x,name = "always_A",variable = "A") <- 1
model(x) <- "additive"
x$estimator <- "tmle"
system.time(x <- run_rtmle(x))
x$estimates
all.equal(summary(tfit)$estimate,x$estimate)
baseline_data(x) <- pop[,.(id,age,sex)]



