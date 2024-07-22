library(testthat)
library(rtmle)
library(data.table)
library(prodlim)
library(targets)
setwd("~/research/SoftWare/rtmle/")
tar_source("R/")
tar_source("~/research/Methods/TMLE_for_breakfast/Ltmle/R/")
source("~/research/Epidemi/Reddie/LEADER/functions/run_ltmle.R")
source("~/research/Epidemi/Reddie/LEADER/functions/summary.runLtmle.R")

set.seed(118)
ld <- simulate_long_data(n = 10000,number_epochs = 20,beta_sum_A_on_Y = -.1,beta_A0_on_Y = 0)
ld[,table(terminal_event)]
ld[,table(table(id))]
plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld[!duplicated(id)]),cause = "Y",xlim = c(0,365.25*3))
plot(prodlim(Hist(terminal_time,terminal_event,cens.code = "C")~A_0,data = ld[!duplicated(id)]),cause = "D",xlim = c(0,365.25*3))
ld[,nth := 1:.N,by = "id"]
ld[,num := .N,by = "id"]
pop <- ld[nth == num,.(id,start = 0, death_date = terminal_time,censored_date = terminal_time,terminal_event,sex,age)]
pop[terminal_event == 0,death_date := NA]
pop[terminal_event == 1,censored_date := NA]
pop[,terminal_event := NULL]
outcome_data <- ld[!duplicated(id) & terminal_event == "Y",.(id,date = terminal_time,stop = NA)]
como_baseline <- ld[!duplicated(id) & L_0 == 1,.(id,date = 0)]
como_data <- ld[event == "L",.(id,date = time)]
como_data <- rbind(como_baseline,como_data)
# random baseline treatment
treatment_baseline <- ld[!duplicated(id) & A_0 == 1,.(id,date = 0)]
treatment_data <- ld[event == "A",.(id,date = time)]
treatment_data <- rbind(treatment_baseline,treatment_data)
w_outcome <- widen_outcome(outcome_name = "Y",pop = pop,intervals = seq(0,2000,30.45*6),como = list("Y" = outcome_data),id = "id")
w_treatment <- widen_treatment(treatment = treatment_data,treatment_name = "A",pop = pop,intervals = seq(0,2000,30.45*6),id = "id")
w_como <- widen_treatment(treatment_name = "L",treatment = como_data,pop = pop,intervals = seq(0,2000,30.45*6),id = "id")
# ltmle
system.time(tfit <- run_ltmle(name_outcome="Y",time_horizon=3,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=pop[,.(id,sex,age)],timevar_data=w_como,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = rep(1,3),name_id = "id",verbose=FALSE,gcomp = FALSE))
system.time(tfit1 <- run_ltmle(name_outcome="Y",time_horizon=3,reduce = FALSE,regimen_data=list("A" = w_treatment),outcome_data=list("Y" = w_outcome),baseline_data=pop[,.(id,sex,age)],timevar_data=w_como,SL.library="glm",censor_others = FALSE,gbounds=c(0,1),abar = list(rep(1,3),rep(0,3)),name_id = "id",verbose=FALSE,gcomp = TRUE))
## summary(tfit)
summary(tfit1)

# rtmle
tar_source("~/tmp/rtmle/R/")

x <- rtmle_init(intervals = 3,name_id = "id",name_outcome = "Y",name_competing = "Dead",name_censoring = "Censored",censored_label = "censored")
treatment_data(x) <- w_treatment
outcome_data(x) <- w_outcome
timevar_data(x) <- w_como
baseline_data(x) <- pop[,.(id,age,sex)]
prepare_data(x) <- "A"
protocol(x,name = "always_A",variable = "A") <- 1
model(x) <- "additive"
x$estimator <- "tmle"
system.time(x <- run_rtmle(x))
x$estimates
all.equal(summary(tfit)$estimate,x$estimate)



