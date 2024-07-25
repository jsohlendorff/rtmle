### intervene.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Jul 19 2024 (11:26) 
##           By: Thomas Alexander Gerds
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
intervene <- function(formula,
                      data,
                      targets,
                      time){
    ## interdata <- copy(model.frame(delete.response(terms(formula(formula))),data,na.action = "na.fail"))
    av <- all.vars(delete.response(terms(formula(formula))))
    if (length(av)>0){
        interdata <- data[,av,with = FALSE]
        ## FIXME: for (tar in targets)
        for (k in 1:time){
            set(interdata,
                j = x$targets[[1]][[paste0("time_",k-1)]]$variable,
                value = x$targets[[1]][[paste0("time_",k-1)]]$value)
        }
    }else{
        interdata <- NULL
    }
    interdata
}
######################################################################
### intervene.R ends here
