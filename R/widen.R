widen <- function(data,
                  name,
                  pop,
                  intervals,
                  fun.aggregate = NULL,
                  id){
    ## stopifnot(all(c(id,"start","end")%in%names(pop)))
    stopifnot(all(c(id,"date")%in%names(data)))
    # better safe then sorry
    data <- copy(data)
    # if no start of followup variable is given
    # we assume all start at zero
    if (("start" %in% names(pop))){
        grid <- pop[,.(date=start+intervals,end = pmin(death_date,censored_date,na.rm = TRUE)),by=id]
    }else{
        grid <- pop[,.(date=intervals,end = pmin(death_date,censored_date,na.rm = TRUE)),by=id]
    }
    grid[,interval:=0:(length(intervals)-1),by=id]
    grid <- pop[,.SD,.SDcols = c(id)][grid,on = id]
    length_interval=unique(round(diff(intervals),0)) 
    grid <- grid[date<=end+length_interval]
    # find out if patient used data in personalized intervals
    w=map_intervals(grid=grid,
                    data=data,
                    name=name,
                    fun.aggregate = fun.aggregate,
                    rollforward=(length_interval - 1),id = id)
    # We subtract 1 to get half-closed intervals incl. left endpoint
    # this way no observation will count for two intervals
    w[]
}



