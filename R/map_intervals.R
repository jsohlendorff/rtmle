map_intervals <- function(grid,
                          data,
                          name,
                          rollforward,
                          values=c(1,0),
                          fill=NA,
                          X.factor=FALSE,
                          id){
    # make sure that only patients who have a date are rolled
    data=data[!is.na(date)]
    if (length(data)==0) return(NULL)
    setkeyv(grid,c(id,"date"))
    setkeyv(data,c(id,"date"))
    data[,X:=values[[1]]]
    grid <- data[grid,roll=rollforward]
    # missing value means no event in this interval
    grid[is.na(grid$X),X:=values[[2]]]
    setkeyv(grid,c(id,"interval"))
    wide <- dcast(grid,
                  formula(paste0(id,"~interval")),
                  value.var="X",
                  sep="_",
                  # FIXME: if more than one treatment event happens
                  fun.aggregate = function(x){1*sum(x)>0},
                  fill=fill)
    if (X.factor) {
        # this is for ltmle censored/uncensored
        for (cc in names(wide)[-1]){
            set(wide,j=cc,value=factor(wide[[cc]],levels=values))
        }
    }
    grid[,X:=NULL]
    data[,X:=NULL]
    # dcast assigns numeric column names when value.var
    # has length one
    setnames(wide,c(id,paste0(name,"_",names(wide)[-1])))
    setkeyv(wide,id)
    wide[]
}
