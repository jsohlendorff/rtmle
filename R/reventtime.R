### reventtime.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jul 12 2024 (09:38) 
## Version: 
## Last-Updated: Jul 18 2024 (13:29) 
##           By: Thomas Alexander Gerds
##     Update #: 59
#----------------------------------------------------------------------
## 
### Comaxtimeentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
reventtime <- function(n,
                       breaks,
                       cumhazard,
                       hazardratio,
                       entrytime = NULL,
                       decimals = NULL){
    if (missing(n)) stop("Please specify sample size 'n'")
    if (any(hazardratio == 0)) stop("Hazardratios cannot contain zeros")
    if (all(breaks>0)) {
        breaks <- c(0,breaks)
        cumhazard <- c(0,cumhazard)
    }
    ## linear approximation
    lapprox <- function(values, breaks, fx){
        pos <- sindex(jump.times = breaks, eval.times = values)
        if (any(pos == 0)) {
            warning("Aproximating outside of the range")
            ## browser()
        }
        maxindex <- which(pos == length(breaks))
        next_pos <- pos + 1
        pos <- pmax(pos,1)
        next_pos[maxindex] <- length(breaks)
        approx_value <- (values - breaks[pos])/(breaks[next_pos] - breaks[pos])
        approx_value[maxindex] <- 0
        res <- approx_value * (fx[next_pos] - fx[pos]) + fx[pos]
        res[is.na(res)] <- tail(fx, 1)
        return(res)
    }
    if (missing(hazardratio)) {
        hazardratio <- rep(1,n)
    }
    maxtime <- tail(breaks, 1)
    if (!is.null(entrytime) && any(entrytime>0)) {
        if (length(entrytime) == 1) entrytime <- rep(entrytime, n)
        entry_cumhazard <- lapprox(values = entrytime,breaks = breaks,fx = cumhazard)
    }
    else {
        entry_cumhazard <- rep(0, n)
    }
    # The distribution of T|X is the same as Lambda^-1(E/HR)
    # where E is expontential with rate 1
    erate <- rexp(n)/hazardratio + entry_cumhazard
    etime <- pmin(lapprox(values = erate,fx = breaks,breaks = cumhazard),maxtime)
    if (!is.null(decimals))
        etime = round(etime,decimals)
    return(etime)
}


######################################################################
### reventtime.R ends here
