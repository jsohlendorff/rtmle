### protocol.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul  3 2024 (13:46) 
## Version: 
## Last-Updated: Jul 19 2024 (11:38) 
##           By: Thomas Alexander Gerds
##     Update #: 7
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' @export
"protocol<-" <- function(x,name,variable,...,value) {
    x$targets[[name]] <- vector(length(x$times),mode = "list")
    names(x$targets[[name]]) <- paste0("time_",x$times)
    for (j in 1:(length(names(x$targets[[name]]))-1)){
        x$targets[[name]][[j]]$variable <- paste0(variable,"_",(j-1))
        x$targets[[name]][[j]]$value <- value
    }
    x
}
######################################################################
### protocol.R ends here
