### simulate_long_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 11 2024 (13:24) 
## Version: 
## Last-Updated: Jul 19 2024 (11:29) 
##           By: Thomas Alexander Gerds
##     Update #: 136
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#'
#' library(data.table)
#' n <- 2000
#' ld <- simulate_long_data(n,time = d$time,number_epochs=20,beta=0.05)
#' ld2 <- simulate_long_data(n,time = d$time,number_epochs=20,beta=-0.05)
#' ld[,table(table(id))]
#' ld[,num:=.N,by="id"]
#' ld[num==1]
#' ld
#' @export
simulate_long_data <- function(n,
                               number_epochs,
                               beta_A0_on_Y,
                               beta_sum_A_on_Y,
                               ...) {
    # mimick a 10-year time frame
    time = seq(1,10*365,1)
    nt <- length(time)
    baseline_hazard = list(A = cumsum(rep(0.0001,nt)),
                           L = cumsum(rep(0.002,nt)),
                           Y = cumsum(rep(0.001,nt)),
                           D = cumsum(rep(0.001,nt)),
                           C = cumsum(rep(0.001,nt)))
    # baseline variables
    pop <- data.table(
        id = 1:n,
        sex=rbinom(n,1,.4),
        age=runif(n,40,90),
        L_0=rbinom(n,1,.17),
        time = numeric(n),
        event = rep("0",n))
    # baseline treatment depends on baseline variables
    ## pop[, A_0:=rbinom(.N,1,lava::expit(0.35+0.1*L_0-0.3*sex-0.01*age))]
    pop[, A_0:=rbinom(.N,1,0.5)]
    people_atrisk <- pop[,.(id,entrytime = time,age,sex,L_0,A_0)]
    people_atrisk[,hazard_ratio_L := exp(0.01*age+0.3*sex)]
    people_atrisk[,hazard_ratio_A := exp(0.0*age-0.0*sex+0.0*L_0+6*A_0)]
    people_atrisk[,hazard_ratio_Y := exp(0.01*age+0.2*sex+0.1*L_0+beta_A0_on_Y*A_0)]
    people_atrisk[,hazard_ratio_D := exp(0.02*age-0.2*sex+0.1*L_0)]
    people_atrisk[,hazard_ratio_C := rep(1,.N)]
    fup <- NULL
    has_terminal <- NULL
    # time loop
    j <- 1
    while (j < number_epochs && nrow(people_atrisk)>0){
        # calculate the time and type of the minimum of latent times to L,A,C,Y,D
        # matrix with latent times
        # FIXME: when both time_A and time_L are smaller than min(time_D,time_Y,time_C)
        #        both could be added but then the first would not have an effect on the second
        ttt = do.call("cbind",list(
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = baseline_hazard[["L"]],hazardratio = people_atrisk$hazard_ratio_L,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = baseline_hazard[["A"]],hazardratio = people_atrisk$hazard_ratio_A,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = baseline_hazard[["C"]],hazardratio = people_atrisk$hazard_ratio_C,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = baseline_hazard[["Y"]],hazardratio = people_atrisk$hazard_ratio_Y,entrytime = people_atrisk$entrytime,decimals = 2),
                                  reventtime(n = nrow(people_atrisk),breaks = time,cumhazard = baseline_hazard[["D"]],hazardratio = people_atrisk$hazard_ratio_D,entrytime = people_atrisk$entrytime,decimals = 2))
                      )
        mins = Rfast::rowMins(ttt,value = FALSE)
        people_atrisk[,event := factor(mins,levels = 1:5,labels = c("L","A","C","Y","D"))]
        people_atrisk[,time := Rfast::rowMins(ttt,value = TRUE)]
        is_terminal <- !(people_atrisk$event%in%c("A","L"))
        # collect followup information
        fup <- rbind(fup,people_atrisk[!is_terminal],fill = TRUE)
        # collect terminal information
        has_terminal <- rbind(has_terminal,people_atrisk[is_terminal,.(id,terminal_time = time,terminal_event = event)])
        # prepare next epoch
        people_atrisk = people_atrisk[!is_terminal]
        # set the entry time to last event time for next epoch
        people_atrisk[,entrytime := time]
        #
        # calculate summaries of the current L and A histories
        #
        current_history = fup[,.(id[1],sum_L = sum(event == "L"),sum_A = sum(event == "A")),by = "id"]
        people_atrisk = current_history[people_atrisk,on = "id"]
        people_atrisk[is.na(sum_A),sum_A := 0]
        people_atrisk[is.na(sum_L),sum_L := 0]
        # the treatment rate is increased the more L
        people_atrisk[,hazard_ratio_A := hazard_ratio_A*exp((0.1*sum_L))]
        people_atrisk[,hazard_ratio_L := hazard_ratio_L]
        # the outcome rate is increased the more L and decreased the more A 
        people_atrisk[,hazard_ratio_Y := hazard_ratio_Y*exp(((0.05*sum_L)+(beta_sum_A_on_Y*sum_A)))]
        ## if (any(people_atrisk$hazard_ratio_Y<0)) browser()
        people_atrisk[,hazard_ratio_D := hazard_ratio_D]
        people_atrisk[,hazard_ratio_C := hazard_ratio_C]
        j = j+1
    }
    # the first epochs (time, event) is not interesting 
    pop[,time := NULL]
    pop[,event := NULL]
    out <- fup[,.(id,event,time)][pop,on = "id"]
    out <- has_terminal[out,on = "id"]
    # if any missing terminal event censor at last event time
    if (any(is.na(out$terminal_event))){
        out[,maxtime := max(time),by = "id"]
        out[is.na(terminal_event),`:=`(terminal_time = maxtime,terminal_event = "C")]
        out[,maxtime := NULL]
    }
    setkey(out,id,time,event)
    # clean up for those who directly have a terminal event
    out[is.na(event),event := terminal_event]
    out[is.na(time),time := terminal_time]
    out[]
}



######################################################################
### simulate_long_data.R ends here
