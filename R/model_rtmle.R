### model_rtmle.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Jul 19 2024 (11:38) 
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
#' @export
"model<-" <- function(x,time,model,...,value) {
    if (is.character(value[[1]]) && value[[1]] == "additive"){
        x$models = additive_formalizer(x = x,Markov = NULL)
    }else{
        x$models[[model]][[paste0("time_",time)]]$formula <- value
    }
    x
}
######################################################################
### model_rtmle.R ends here
