#' @exportS3Method
#'
print.qfrm <- function(x, digits = getOption("digits"),
                       show_range = !is.null(x$error_bound),
                       prefix = "\t", ...) {
    stat <- x$statistic
    errorb <- x$error_bound
    exact <- isTRUE(attr(errorb, "exact"))
    cat("\n")
    cat(strwrap("Moment of ratio of quadratic form", prefix = prefix), sep = "\n")
    cat("\n")
    out <- character()
    if(!is.null(stat)) {
        out <- c(out, paste("Moment =", format(stat,
            digits = max(1L, digits - 2L))))
    }
    if(length(errorb) > 0 && !xor(all(is.na(errorb)), all(is.nan(errorb))) && !exact) {
        out <- c(out, paste("Error =", format(errorb,
            digits = max(1L, digits - 2L))))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if(exact) {
        cat("This value is exact\n")
    } else if(length(errorb) > 0 && xor(all(is.na(errorb)), all(is.nan(errorb)))) {
        cat("Error bound unavailable; recommended to inspect plot() of this xect\n")
    } else if(show_range) {
        if(isTRUE(attr(errorb, "twosided"))) {
            ra <- sort(c(stat - errorb, stat + errorb))
        } else {
            ra <- sort(c(stat, stat + errorb))
        }
        cat("Possible range:\n ",
        paste(format(ra, digits = digits), collapse = " "), "\n", sep = "")
    }
    if(isTRUE(attr(errorb, "singular"))) {
        cat(paste("Note: Argument matrix (numerically) singular, so error bound is unreliable\n"))
    }
    if(isTRUE(attr(errorb, "alphaout"))) {
        cat(paste("Note: Adjustment parameter(s) alpha above 1, so error bound is unreliable\n"))
    }
    cat("\n")
    invisible(x)
}

#' @exportS3Method
#'
plot.qfrm <- function(x, add_error = length(errseq) > 0,
                      add_legend = add_error,
                      ylim = sum(ansseq) * ylim_f, ylim_f = c(0.9, 1.1),
                      xlab = "Order of polynomials", ylab = "Moment of ratio",
                      col_m = "royalblue4", col_e = "tomato",
                      lwd_m = 1, lwd_e = 1, lty_m = 1, lty_e = 2,
                      pos_leg = "topright", ...) {
    if(!requireNamespace("graphics", quietly = TRUE)) {
        stop("Package \"graphics\" is required for plot.qfrm")
    }
    ansseq <- x$res_seq
    errseq <- x$err_seq
    cumseq <- cumsum(ansseq)
    try(plot(seq_along(ansseq) - 1L, cumseq, type = "l", col = col_m,
         ylim = ylim, xlab = xlab,
         ylab = ylab, lwd = lwd_m, lty = lty_m, ...))
    if(add_error) {
        graphics::lines(seq_along(ansseq) - 1L, errseq + cumseq, col = col_e, lwd = lwd_e, lty = lty_e)
        if(isTRUE(attr(errseq, "twosided"))) {
            graphics::lines(-errseq + cumseq, col = col_e, lwd = lwd_e, lty = lty_e)
        }
    }
    if(add_legend) {
        graphics::legend(pos_leg, legend = c("Moment", "Error bound"),
               col = c(col_m, col_e), lwd = c(lwd_m, lwd_e),
               lty = c(lty_m, lty_e))
    }
}
