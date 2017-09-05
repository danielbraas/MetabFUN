NA_function <- function(x) if (sum(is.na(x)) < length(x)) return(x=x) else return(x=rep(0, length(x)))
