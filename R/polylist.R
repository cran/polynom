polylist <- function(...) {
    P <- lapply(list(...), as.polynomial)
    structure(P, class = "polylist")
}

deriv.polylist <- function(p, ...) 
    structure(lapply(p, deriv), class = class(p))

integral.polylist <- function(p, ...) {
    result <- lapply(p, integral, ...)
    if (length(result) > 0 && is.polynomial(result[[1]]))
        class(result) <- class(p)
    result
}

plot.polylist <- function(p, xlim = 0:1, ylim = range(Px), type = "l",
                          len = 100, ...) {
    if(missing(xlim)) {
        ## try to cover the "interesting" region
        xlim <- range(Re(unlist(lapply(p, summary.polynomial))))
    }
    if(any(is.na(xlim))) {
        warning("summary of polynomial fails. Using nominal xlim")
        xlim <- 0:1
    }
    if(diff(xlim) == 0)
        xlim <- xlim + c(-1, 1)/2
    if(length(xlim) > 2)
        x <- xlim
    else {
        eps <- diff(xlim)/100
        xlim <- xlim + c( - eps, eps)
        x <- seq(xlim[1], xlim[2], len = len)
    }
    Px <- unlist(lapply(p, predict.polynomial, x))
    if(!missing(ylim))
        Px[Px < ylim[1]] <- Px[Px > ylim[2]] <- NA
    plot(cbind(x, Px), xlab = "x", ylab = "P(x)", type = "n",
         xlim = xlim, ylim = ylim, ...)
    for(i in seq(along = p))
        lines(p[[i]], lty = i)
    invisible()
}

print.polylist <- function(x, ...) {
    cat("List of polynomials:\n")
    y <- x
    x <- unclass(x)
    NextMethod()
    invisible(y)
}
