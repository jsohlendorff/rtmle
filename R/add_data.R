### add_data.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 25 2024 (11:24) 
## Version: 
## Last-Updated: Jul 25 2024 (11:51) 
##           By: Thomas Alexander Gerds
##     Update #: 5
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' @export
"add_data<-" <- function(x,...,value){
    x$long_data <- value
    ## x$data <- value
    x
}


######################################################################
### add_data.R ends here
