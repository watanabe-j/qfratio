##### methods.qfrm (documentation) #####
#' Methods for qfrm and qfpm objects
#'
#' This package defines straightforward \code{print} and \code{plot} methods
#' for \code{qfrm} and \code{qfpm} objects which result from the
#' \code{\link{qfrm}}, \code{\link{qfmrm}}, and \code{\link{qfpm}} functions.
#'
#' The \code{print} methods simply display the moment (typically
#' a truncated sum), along with its error bound (when available).
#'
#' The \code{plot} method is designed for quick inspection of the profile of
#' a series expression along varying polynomial orders.
#' When the object has a sequence for error bounds, this is also shown
#' with a broken line (by default).
#' When the object has an exact moment (i.e., resulting from
#' \code{\link{qfrm_ApIq_int}()} or the \code{\link{qfpm}} functions), a warning
#' is thrown because inspection of the plot will not be required in this case.
#'
#' @param x
#'   \code{qfrm} or \code{qfpm} object
#' @param digits
#'   Number of significant digits to be printed.
#'   By default 2 digits larger than the setting in options.
#' @param show_range
#'   Logical to specify whether the possible range for the moment
#'   is printed (when available).  Default \code{TRUE} when available.
#' @param prefix
#'   String passed to \code{\link{strwrap}}, as in
#'   \code{\link[stats]{print.power.htest}}
#' @param add_error
#'   Logical to specify whether the sequence of error bounds is plotted
#'   (when available).  Default \code{TRUE} when available.
#' @param add_legend
#'   Logical to specify whether a legend is added.  Turned on by default
#'   when \code{add_error = TRUE}.
#' @param ylim,ylim_f
#'   \code{ylim} is passed to \code{\link[graphics]{plot.default}};
#'   By default, this is automatically set to \code{ylim_f} times
#'   the terminal value of the sequence expression (\code{sum(x$res_seq)}).
#'   \code{ylim_f} is by default \code{c(0.9, 1.1)}.
#' @param xlab,ylab
#'   Passed to \code{\link[graphics]{plot.default}}
#' @param col_m,col_e,lwd_m,lwd_e,lty_m,lty_e
#'   \code{col}, \code{lwd}, and \code{lty} to plot the sequences of
#'   the moment (\code{***_m}) and its error bound (\code{***_e})
#' @param pos_leg
#'   Position of the legend, e.g., \code{"topright"}, \code{"bottomright"},
#'   passed as the first argument for \code{\link[graphics]{legend}}
#' @param ...
#'   In the \code{plot} methods, passed to \code{\link[graphics]{plot.default}}.
#'   In the \code{print} methods, ignored (retained for the compatibility
#'   with the generic method).
#'
#' @name methods.qfrm
#'
#' @examples
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(1:nv)
#' mu <- rep.int(1, nv)
#'
#' res1 <- qfrm(A, B, p = 3, mu = mu)
#' print(res1)
#' print(res1, digits = 5)
#' print(res1, digits = 10)
#'
#' ## Default plot: ylim too narrow to see the error bound at this m
#' plot(res1)
#'
#' ## With extended ylim
#' plot(res1, ylim_f = c(0.8, 1.2), pos_leg = "topleft")
#'
#' ## In this case, it is easy to increase m
#' (res2 <- qfrm(A, B, p = 3, mu = mu, m = 200))
#' plot(res2)
#'
NULL

#' @rdname methods.qfrm
#' @order 1
#'
#' @exportS3Method
#'
print.qfrm <- function(x, digits = getOption("digits") + 2,
                       show_range = !is.null(x$error_bound),
                       prefix = "\t", ...) {
    stat <- x$statistic
    errorb <- x$error_bound
    exact <- isTRUE(attr(errorb, "exact"))
    cat("\n")
    cat(strwrap("Moment of ratio of quadratic forms", prefix = prefix), sep = "\n")
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
        cat("Error bound unavailable; recommended to inspect plot() of this object\n")
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

#' @rdname methods.qfrm
#' @order 3
#'
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
    if(isTRUE(attr(errseq, "exact"))) {
        warning("plot method for this class is for inspection of truncated sums.\n  ",
                "This object has an exact moment, in which case plot is moot.")
    }
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

#' @rdname methods.qfrm
#' @order 2
#'
#' @exportS3Method
#'
print.qfpm <- function(x, digits = getOption("digits") + 2,
                       prefix = "\t", ...) {
    stat <- x$statistic
    cat("\n")
    cat(strwrap("Moment of (product of) quadratic form(s)", prefix = prefix), sep = "\n")
    cat("\n")
    out <- paste("Moment =", format(stat, digits = max(1L, digits - 2L)))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    cat("This value is exact\n")
    cat("\n")
    invisible(x)
}
