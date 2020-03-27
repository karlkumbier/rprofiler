# This functions is a modified version of the R-base function ks.test. It
# computes a signed KS statistic measuring the maximum difference between a
# distribution x and some reference distribution y.
ksTest <- function (x, y) {
    x <- x[!is.na(x)]
    n.x <- length(x)
    if (n.x < 1L) return(NA) 
        
    y <- y[!is.na(y)]
    n.y <- length(y)
    if (n.y < 1L) return(NA) 
    
    w <- c(x, y)
    z <- cumsum(ifelse(order(w) <= n.x, 1 / n.x, -1 / n.y))
    
    if (length(unique(w)) < (n.x + n.y)) {
      z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
    }
    
    statistic <- c(max(z), min(z))
    statistic <- statistic[which.max(abs(statistic))]
    return(statistic)
}
