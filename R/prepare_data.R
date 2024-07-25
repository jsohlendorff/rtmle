### prepare_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 19 2024 (10:07) 
## Version: 
## Last-Updated: Jul 19 2024 (11:43) 
##           By: Thomas Alexander Gerds
##     Update #: 3
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
"prepare_data<-" <- function(x,subset_id = NULL,...,value){
    x$name_regimen = value
    time_horizon = max(x$times)
    K = length(x$times)
    stopifnot(inherits(x$treatment_data,"data.table"))
    stopifnot(inherits(x$outcome_data,"data.table"))
    stopifnot(match(x$name_id,names(x$treatment_data),nomatch = 0)>0)
    stopifnot(match(x$name_id,names(x$outcome_data),nomatch = 0)>0)
    data.table::setkeyv(x$treatment_data,x$name_id)
    data.table::setkeyv(x$outcome_data,x$name_id)
    work_data=x$outcome_data[x$treatment_data,on = x$name_id]
    # deal with outcome/death/censored at index
    Y_0 = match(paste0(x$name_outcome,"_",0),names(work_data))
    D_0 = match(paste0(x$name_competing,"_",0),names(work_data))
    C_0 = match(paste0(x$name_censoring,"_",0),names(work_data))
    excluded = FALSE
    if (!is.na(Y_0)){
        if(!is.na(D_0)&!is.na(C_0)){
            excluded <- (work_data[[Y_0]]%in%1)|(work_data[[D_0]]%in%1)|(work_data[[C_0]]%in%x$censored_label)
        }
        if(!is.na(D_0)){
            excluded <- (work_data[[Y_0]]%in%1)|(work_data[[D_0]]%in%1)
        }
        if(!is.na(C_0)){
            excluded <- (work_data[[Y_0]]%in%1)|(work_data[[C_0]]%in%x$censored_label)
        }
    }else{
        if(!is.na(D_0)&!is.na(C_0)){
            excluded <- (work_data[[D_0]]%in%1)|(work_data[[C_0]]%in%x$censored_label)
        }
        if(!is.na(D_0)){
            excluded <- (work_data[[D_0]]%in%1)
        }
        if(!is.na(C_0)){excluded <- (work_data[[C_0]]%in%x$censored_label)
        }
    }
    if (any(excluded)){
        work_data = work_data[!excluded]
        stopifnot(nrow(work_data)>0)
    }
    # adding the baseline covariates
    if (!is.null(x$baseline_data)){
        stopifnot(inherits(x$baseline_data,"data.table"))
        stopifnot(match(x$name_id,names(x$baseline_data),nomatch = 0)>0)
        work_data=x$baseline_data[work_data,on = x$name_id]
    }
    # add time covariates
    # first remove outcome if overlap
    if (length(x$timevar_data)>0){
        stopifnot(inherits(x$timevar_data,"data.table"))
        stopifnot(match(x$name_id,names(x$timevar_data),nomatch = 0)>0)
        if (length((outcome_overlap <- grep(paste0(x$name_outcome,"_"),names(x$timevar_data)))>0)){
            timevar_data <- x$timevar_data[,-outcome_overlap, with=FALSE]}
        else
            timevar_data = x$timevar_data
        data.table::setkeyv(timevar_data,x$name_id)
        work_data=timevar_data[work_data, on = x$name_id]
        name_time_covariates = unlist(lapply(grep("_0",names(timevar_data),value=TRUE),
                                             function(x){substring(x,0,nchar(x)-2)}))
    }else{
        name_time_covariates <- NULL
    }
    name_baseline_covariates = setdiff(names(x$baseline_data),x$name_id)
    # sorting the variables for LTMLE
    work_data = work_data[,c(x$name_id, intersect(c(name_baseline_covariates,unlist(sapply(x$times, function(timepoint){
        if(timepoint == 0){
            paste0(c(name_time_covariates, x$name_regimen),"_",timepoint)
        } else{
            if(timepoint != x$times[K]){
                paste0(c(x$name_censoring, x$name_outcome, x$name_competing, name_time_covariates, x$name_regimen),"_",timepoint)
            } else {
                paste0(c(x$name_censoring, x$name_outcome),"_",timepoint)
            }
        }
    }))), names(work_data))), with = FALSE]
    # subset data
    if(length(subset_id)>0){
        subset_dt = data.table(ID = subset_id)
        setnames(subset_dt,"ID",x$name_id)
        work_data = work_data$data[subset_dt,on = x$name_id]
    }
    # label the variables that are constant in the subset data
    same = sapply(work_data, function(x){length(unique(x))==1})
    if(sum(same)>0){
        constant_variables <- names(work_data)[same]
    } else{
        constant_variables <- NULL}
    name_baseline_covariates <- intersect(name_baseline_covariates,names(work_data))
    ## check
    if(length(x$name_censoring)>0){
        for(col in sapply(x$times[-1], function(timepoint){paste0(x$name_censoring,"_",timepoint)})){
            set(work_data, j = col, value=as.factor(ifelse(work_data[[col]]%in%x$censored_label,"censored","uncensored")))
        }
    }
    ## Manipulation of the event nodes
    A_nodes = unlist(lapply(x$times[-K], function(time){paste0(x$name_regimen, "_", time)}))
    Y_nodes = unlist(lapply(x$times[-1], function(time){paste0(x$name_outcome, "_", time)}))
    D_nodes = unlist(lapply(x$times[-c(1,K)], function(time){paste0(x$name_competing, "_", time)}))
    C_nodes = unlist(lapply(x$times[-1], function(time){paste0(x$name_censoring, "_", time)}))
    # if A,B then B_0 is obsolete because A0=1-B0
    if (length(x$name_regimen) == 2)
        A_nodes <- A_nodes[A_nodes != paste0(x$name_regimen[[2]],"_0")]
    A_nodes_position = match(A_nodes, names(work_data))
    Y_nodes_position = match(Y_nodes, names(work_data))
    D_nodes_position = match(D_nodes, names(work_data))
    C_nodes_position = match(C_nodes, names(work_data))
    ## Adjust data depending on censoring/event/competing event with NA
    for(q in 1:(K-1)){
        if(q<(K-1)){
            has_outcome_or_death_and_censored = (((work_data[[Y_nodes_position[[q]]]]%in%"1")|(work_data[[D_nodes_position[[q]]]]%in%"1"))&
                                                 (work_data[[C_nodes_position[[q]]]]%in%"censored"))
        } else{
            has_outcome_or_death_and_censored = ((work_data[[Y_nodes_position[[q]]]]%in%1)&(work_data[[C_nodes_position[[q]]]]%in%"censored"))
        }
        if(any(has_outcome_or_death_and_censored)){
            set(work_data,j=C_nodes_position[[q]],i=which(has_outcome_or_death_and_censored),value="uncensored")
        }
    }
    ## All nodes (except outcome and competing risk) should be NA after an event (outcome or death)
    if(time_horizon!= 1){
        for(k in Y_nodes_position[-(K-1)]){
            later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
            later_Y_nodes=intersect((k+1):Y_nodes_position[length(Y_nodes_position)],Y_nodes_position)
            # later_C_nodes=intersect((k+1):NCOL(work_data),C_nodes_position)
            # later_D_nodes=intersect((k+1):NCOL(work_data),D_nodes_position)
            if(any(has_outcome <- (work_data[[k]]%in%1))){
                for(l in later_nodes) {set(work_data,j=l,i=which(has_outcome),value=NA)}
                for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_outcome),value=1)}
                # for(l in later_C_nodes) {set(work_data,j=l,i=which(has_outcome),value="uncensored")}
                # for(l in later_D_nodes) {set(work_data,j=l,i=which(has_outcome),value=0)}
            }
        }
        if(length(x$name_competing)>0){
            for(k in D_nodes_position){
                later_nodes=setdiff((k+1):NCOL(work_data),Y_nodes_position)
                # Later outcome event nodes are set to 0
                later_Y_nodes=intersect((k+1):NCOL(work_data),Y_nodes_position)
                # later_C_nodes=intersect((k+1):NCOL(work_data),C_nodes_position)
                # later_D_nodes=intersect((k+1):NCOL(work_data),D_nodes_position)
                if(any(has_died <- (work_data[[k]]%in%1))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_died),value=NA)}
                    for(l in later_Y_nodes) {set(work_data,j=l,i=which(has_died),value=0)}
                    # for(l in later_C_nodes) {set(work_data,j=l,i=which(has_died),value="uncensored")}
                    # for(l in later_D_nodes) {set(work_data,j=l,i=which(has_died),value=1)}
                }
            }
        }
        ## All nodes should be NA as soon as censoring has occurred
        if(length(x$name_censoring)>0){
            for(k in C_nodes_position){
                later_nodes=(k+1):NCOL(work_data)
                ## print(k)
                ## print(names(work_data)[later_nodes])
                if(any(has_censored <- (work_data[[k]]%in%"censored"))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }else{
        # 
        # time_horizon = 1 we set the outcome to NA in case of censored
        #                  and same for competing risks
        #
        if(length(x$name_censoring)>0){
            for(k in C_nodes_position){
                later_nodes=(k+1):NCOL(work_data)
                ## print(names(work_data)[later_nodes])
                if(any(has_censored <- (work_data[[k]]%in%"censored"))){
                    for(l in later_nodes) {set(work_data,j=l,i=which(has_censored),value=NA)}
                }
            }
        }
    }
    L_nodes <- c(sapply(x$times, function(k) {paste0(c(name_time_covariates), "_", k)}))
    L_nodes <- L_nodes[match(L_nodes, names(work_data),nomatch = 0)!=0]
    if(length(x$name_censoring)==0) {C_nodes = NULL}
    x$data <- work_data[]
    x$name_time_covariates <- name_time_covariates
    x$name_baseline_covariates <- name_baseline_covariates
    x$name_constant_variables <- constant_variables
    x$Anodes = intersect(A_nodes, names(work_data))
    x$Cnodes = intersect(C_nodes, names(work_data))
    x$Dnodes = intersect(D_nodes, names(work_data))
    x$Lnodes = intersect(L_nodes, names(work_data)) 
    x$Ynodes = intersect(Y_nodes, names(work_data))
    x
}
######################################################################
### prepare_data.R ends here
