### run_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  1 2024 (09:11) 
## Version: 
## Last-Updated: Jul 19 2024 (11:44) 
##           By: Thomas Alexander Gerds
##     Update #: 114
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
run_rtmle <- function(x,speed = FALSE,...){
    #
    # g-part: fit nuisance parameter models for propensity and censoring
    #
    x$IC <- numeric(nrow(x$data))
    intervention_probs <- data.table(ID = x$data[[x$name_id]])
    intervention_probs <- cbind(intervention_probs,matrix(1,nrow = NROW(x$data),ncol = length(x$Anodes)+length(x$Cnodes)))
    setnames(intervention_probs,c(x$name_id,c(rbind(x$Anodes,x$Cnodes))))
    for (m in c("propensity","censoring")){
        for (j in 1:length(x$models[[m]])){
            # data at-risk at the beginning of the interval
            if (j > 1){
                outcome_free <- x$data[[paste0(x$name_outcome,"_",j-1)]]%in%0
                # competing risks
                if (length(x$Dnodes)>0) outcome_free <- outcome_free&x$data[[paste0(x$name_competing,"_",j-1)]]%in%0
                uncensored <- x$data[[paste0(x$name_censoring,"_",j-1)]]%in%"uncensored"
            }else{
                outcome_free <- rep(TRUE,nrow(x$data))
                uncensored <- rep(TRUE,nrow(x$data))
            }
            if (any(outcome_free&uncensored)){
                # fit the propensity/censoring regression model
                if (!is.null(ff <- x$models[[m]][[j]]$formula)){
                    if (speed & !inherits(try(
                                     x$models[[m]][[j]]$fit <- speedglm::speedglm(formula = ff,data = x$data[outcome_free&uncensored],family = binomial(),maxit = 100),silent = TRUE),
                                     "try-error")){
                    } else{
                        if (inherits(try(x$models[[m]][[j]]$fit <- glm(formula = ff,data = x$data[outcome_free&uncensored],family = binomial()),
                                         silent = TRUE),"try-error"))
                            browser()
                    }
                    # predict the propensity score/1-probability of censored
                    # intervene according to protocol for targets
                    # FIXME: predict also those censored?
                    intervened_data = intervene(formula = ff,
                                                data = x$data[outcome_free],
                                                targets = x$targets,
                                                time = j)
                    if (m == "propensity")
                        current_node <- x$Anodes[[j]]
                    else
                        current_node <- x$Cnodes[[j]]
                    intervention_probs[outcome_free][[current_node]] <- predict(x$models[[m]][[j]]$fit,newdata = intervened_data,type = "response")
                }
            }
        }
    }
    x$cumulative_intervention_probs <- matrixStats::rowCumprods(as.matrix(intervention_probs[,-1,with = FALSE]))
    # fixme for target in targets
    target <- 1
    intervention_match <- matrix(0,ncol = length(x$Anodes),nrow = nrow(x$data))
    for(j in x$times[-c(length(x$times))]){
        target_j <- x$targets[[target]][[paste0("time_",j)]]
        if (j == 0)
            intervention_match[,j+1] <- previous <- (x$data[[target_j$variable]] %in% c(target_j$value,NA))
        else
            intervention_match[,j+1] <- previous <- previous*(x$data[[target_j$variable]] %in% c(target_j$value,NA))
    }
    # fixme: instead of Anodes should play safe and use sapply(x$targets[[target]],function(t)t$variable)?
    colnames(intervention_match) <- x$Anodes
    x$intervention_match = intervention_match
    # 
    # Q-part: loop backwards in time through iterative condtional expectations
    #
    # the first step is the outcome regression of the last time interval
    for (j in rev(x$times)[-length(x$times)]){
        # formula
        interval_outcome_formula = formula(x$models[["outcome"]][[paste0(x$name_outcome,"_",j)]]$formula)
        #  at the last time interval the observed outcome is used
        if (j != max(x$times)) {
            interval_outcome_formula <- update(interval_outcome_formula,"predicted_outcome~.")
            Yhat <- x$data[["predicted_outcome"]]
        }else{
            Yhat <- x$data[[all.vars(interval_outcome_formula)[[1]]]]
        }
        ## print(interval_outcome_formula)
        ## print(head(Yhat))
        # data at-risk at the beginning of the interval
        if (j > 1){
            outcome_free <- x$data[[paste0(x$name_outcome,"_",j-1)]]%in%0
            # competing risks
            if (length(x$Dnodes)>0) outcome_free <- outcome_free&x$data[[paste0(x$name_competing,"_",j-1)]]%in%0
            uncensored <- x$data[[paste0(x$name_censoring,"_",j-1)]]%in%"uncensored"
        }else{
            outcome_free <- rep(TRUE,nrow(x$data))
            uncensored <- rep(TRUE,nrow(x$data))
        }
        # fit outcome regression
        if (speed & !inherits(try(fit_last <- speedglm::speedglm(formula = interval_outcome_formula,
                                                                 data = x$data[outcome_free&uncensored],
                                                                 family = binomial(),
                                                                 maxit = 100),silent = TRUE),
                              "try-error")){
        } else
            fit_last <- glm(formula = interval_outcome_formula,
                            data = x$data[outcome_free&uncensored],
                            family = "binomial")
        # save fitted object
        x$models[["outcome"]][[j]]$fit <- fit_last
        # intervene according to protocol for targets
        intervened_data = intervene(formula = interval_outcome_formula,
                                    ## data = x$data[outcome_free&uncensored],
                                    data = x$data,
                                    targets = x$targets,
                                    time = j)
        # set predicted value as outcome for next regression
        ## FIXME: if (target$estimator == "tmle")
        y <- predict(fit_last,newdata = intervened_data ,type = "response")
        current_cnode = x$data[[paste0(x$name_censoring,"_",j)]]
        ## print("before k+1")
        ## print(head(Yhat))
        ## print("before k")
        ## print(head(y))
        if (length(x$estimator) == 0 || x$estimator == "tmle"){
            ## browser()
            W = update_Q(Y = Yhat,
                         logitQ = lava::logit(y),
                         cum.g = x$cumulative_intervention_probs[,match(paste0("Censored_",j),colnames(x$cumulative_intervention_probs))], 
                         uncensored_undeterministic = outcome_free & (current_cnode%in%"uncensored"),
                         intervention.match = intervention_match[,x$targets[[target]][[paste0("time_",j-1)]]$variable])
        }else{
            W <- predict(fit_last,
                         newdata = intervened_data ,
                         type = "response")
        }
        ## print("after")
        ## print(head(W))
        # calculate contribution to influence function
        h.g.ratio <- 1/x$cumulative_intervention_probs[,match(paste0("Censored_",j),colnames(x$cumulative_intervention_probs))]
        index <- (current_cnode%in%"uncensored") & intervention_match[,x$targets[[target]][[paste0("time_",j-1)]]$variable]
        if (any(h.g.ratio[index] != 0)) {
            x$IC[index] <- x$IC[index] + (Yhat[index] - W[index]) * h.g.ratio[index]
        }
        ## curIC <- CalcIC(Qstar.kplus1, Qstar, update.list$h.g.ratio,
        ## uncensored, intervention.match, regimes.with.positive.weight)
        # prepare next iteration
        set(x = x$data,
            i = which(outcome_free&uncensored), # atrisk: not censored, event-free, no competing risk
            j = "predicted_outcome",
            value = W[which(outcome_free&uncensored)])
        # for those who have had an event or died or censored earlier
        set(x = x$data,
            i = which(!outcome_free|!uncensored), # not atrisk
            j = "predicted_outcome",
            # fixme j or j-1?
            value = x$data[[paste0(x$name_outcome,"_",j)]][which(!outcome_free|!uncensored)])
    }
    # g-formula and tmle estimator
    x$estimates <- mean(x$data$predicted_outcome)
    x$IC <- x$IC + x$data$predicted_outcome - mean(x$data$predicted_outcome)
    x
}
######################################################################
### run_rtmle.R ends here
