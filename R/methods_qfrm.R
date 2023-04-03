##### new_qfrm #####
#' Construct qfrm object
#'
#' These are internal \dQuote{constructor} functions used to make \code{qfrm}
#' and \code{qfpm} objects, which are used as a return value from the
#' \code{\link{qfrm}}, \code{\link{qfmrm}}, and \code{\link{qfpm}} functions.
#'
#' @param statistic
#'   Terminal value (partial sum) for the moment.  When missing,
#'   obtained as \code{sum(terms)}.
#' @param error_bound
#'   Terminal error bound.  When missing,
#'   obtained as \code{seq_error[length(seq_error)]}.
#' @param terms
#'   Terms in series expression for the moment along varying polynomial degrees
#' @param seq_error
#'   Vector of error bounds corresponding to \code{cumsum(terms)}
#' @param exact,twosided,alphaout,singular_arg
#'   Logicals used to append attributes to the resultant error bound
#'   (see \dQuote{Value})
#' @param diminished
#'   Logical used to append attribute to the resultant statistic and terms
#'   (see \dQuote{Value})
#' @param class
#'   Character vector to (pre-)append classes to the return value
#' @param ...
#'   Additional arguments for accommodating subclasses
#'
#' @return
#' \code{new_qfrm()} and \code{new_qfpm()} return a list of class \code{qfrm}
#' and \code{c(qfpm, qfrm)}, respectively.  These classes are defined for
#' the \code{print} and \code{plot} methods.
#'
#' The return object is a list of 4 elements which are intended to be:
#' \itemize{
#'   \item{\code{$statistic}: }{evaluation result (\code{sum(terms)})}
#'   \item{\code{$terms}: }{vector of \eqn{0}th to \eqn{m}th order terms}
#'   \item{\code{$error_bound}: }{error bound of \code{statistic}}
#'   \item{\code{$seq_error}: }{vector of error bounds corresponding to
#'                              partial sums (\code{cumsum(terms)})}
#' }
#' When the result is exact, \code{$terms} can be of length 1 and equal to
#' \code{$statistic}.  This is always the case for the \code{qfpm} class.
#'
#' When the relevant flags are provided in the constructor, \code{$error_bound}
#' and \code{$seq_error} have the following attributes which control behaviors
#' of the \code{print} and \code{plot} methods:
#' \itemize{
#'   \item{\code{"exact"}: }{indicates whether the moment is exact}
#'   \item{\code{"twosided"}: }{indicates whether the error bounds are
#'                two-sided}
#'   \item{\code{"alphaout"}: }{indicates whether any of the scaling factors
#'                (\code{alphaA}, \code{alphaB}, \code{alphaD}) is outside
#'                \eqn{(0, 1]}, when error bound does not strictly hold}
#'   \item{\code{"singular"}: }{indicates whether the relevant argument matrix
#'                is (numerically) singular, in which case the error bound is
#'                invalid}
#' }
#' Similarly, when \code{diminished = TRUE}, \code{$statistic} and \code{$terms}
#' have the attribute \code{"diminished"} being \code{TRUE}, which indicates
#' that numerical underflow/diminishing happened during scaling
#' (see \dQuote{Scaling} in \code{\link{d1_i}}).
#'
#' @seealso
#' \code{\link{qfrm}}, \code{\link{qfmrm}}, \code{\link{qfpm}}: functions
#' that return objects of these classes
#'
#' \code{\link{methods.qfrm}}: the \code{print} and \code{plot} methods
#'
#' @name new_qfrm
#'
new_qfrm <- function(statistic, error_bound = NULL,
                     terms = statistic, seq_error = NULL,
                     exact = FALSE, twosided = FALSE, alphaout = FALSE,
                     singular_arg = FALSE, diminished = FALSE, ...,
                     class = character()) {
    if(missing(statistic) && !missing(terms)) {
        statistic <- sum(terms)
    }
    if(missing(error_bound) && !missing(seq_error)) {
        error_bound <- seq_error[length(seq_error)]
    }
    if(!is.null(seq_error)) {
        if(is.na(error_bound) && all(is.na(error_bound))
           && !all(is.nan(error_bound))) {
            seq_error <- NULL
        }
    }
    if(isTRUE(exact)) {
        if(!is.null(error_bound)) attr(error_bound, "exact") <- TRUE
        if(!is.null(seq_error)) attr(seq_error, "exact") <- TRUE
    }
    if(isTRUE(twosided)) {
        if(!is.null(error_bound)) attr(error_bound, "twosided") <- TRUE
        if(!is.null(seq_error)) attr(seq_error, "twosided") <- TRUE
    }
    if(isTRUE(alphaout)) {
        if(!is.null(error_bound)) attr(error_bound, "alphaout") <- TRUE
        if(!is.null(seq_error)) attr(seq_error, "alphaout") <- TRUE
    }
    if(isTRUE(singular_arg)) {
        if(!is.null(error_bound)) attr(error_bound, "singular") <- TRUE
        if(!is.null(seq_error)) attr(seq_error, "singular") <- TRUE
    }
    if(isTRUE(diminished)) {
        if(!is.null(terms)) attr(terms, "diminished") <- TRUE
        if(!is.null(statistic)) attr(statistic, "diminished") <- TRUE
    }
    structure(list(statistic = statistic, error_bound = error_bound,
                   terms = terms, seq_error = seq_error),
              class = c(class, "qfrm"))
}

#' @rdname new_qfrm
#'
new_qfpm <- function(statistic, exact = TRUE, ..., class = character()) {
    error <- 0
    new_qfrm(statistic, error, statistic, error, exact = exact, ...,
             class = c(class, "qfpm"))
}

##### methods.qfrm (documentation) #####
#' Methods for qfrm and qfpm objects
#'
#' Straightforward \code{print} and \code{plot} methods are defined
#' for \code{qfrm} and \code{qfpm} objects which result from the
#' \code{\link{qfrm}}, \code{\link{qfmrm}}, and \code{\link{qfpm}} functions.
#'
#' The \code{print} methods simply display the moment \code{x$statistic}
#' (typically a partial sum), its error bound \code{x$error_bound}
#' (when available), and the possible range of the moment
#' (\code{x$statistic} to \code{x$statistic + x$error_bound} in case of
#' one-sided error bound; \code{x$statistic - x$error_bound} to
#' \code{x$statistic + x$error_bound} in case of two-sided).
#'
#' The \code{plot} method is designed for quick inspection of the profile of
#' the partial sum of the series along varying orders
#' \code{cumsum(x$terms)}.  When the object has a sequence for error bounds
#' \code{x$seq_error}, this is also shown with a broken line (by default).  When
#' the object has an exact moment (i.e., resulting from
#' \code{\link{qfrm_ApIq_int}()} or the \code{\link{qfpm}} functions), a message
#' is thrown to tell inspection of the plot will not be required in this case.
#'
#' @param x
#'   \code{qfrm} or \code{qfpm} object
#' @param digits
#'   Number of significant digits to be printed.
#' @param show_range
#'   Logical to specify whether the possible range for the moment
#'   is printed (when available).  Default \code{TRUE} when available.
#' @param add_error
#'   Logical to specify whether the sequence of error bounds is plotted
#'   (when available).  Default \code{TRUE} when available.
#' @param add_legend
#'   Logical to specify whether a legend is added.  Turned on by default
#'   when \code{add_error = TRUE}.
#' @param ylim,ylim_f
#'   \code{ylim} is passed to \code{\link[graphics]{plot.default}}.  By default,
#'   this is automatically set to \code{ylim_f} times the terminal value of
#'   the seies expression (\code{x$statistic}).  \code{ylim_f} is
#'   by default \code{c(0.9, 1.1)}.
#' @param xlab,ylab
#'   Passed to \code{\link[graphics]{plot.default}}
#' @param col_m,col_e,lwd_m,lwd_e,lty_m,lty_e
#'   \code{col}, \code{lwd}, and \code{lty} to plot the sequences of
#'   the moment (\code{***_m}) and its error bound (\code{***_e})
#' @param pos_leg
#'   Position of the legend, e.g., \code{"topright"}, \code{"bottomright"},
#'   passed as the first argument for \code{\link[graphics]{legend}}
#' @param ...
#'   In the \code{plot} methods, passed to
#'   \code{\link[graphics]{plot.default}}.  In the \code{print} methods,
#'   ignored (retained for the compatibility with the generic method).
#'
#' @return
#' The \code{print} method invisibly returns the input.
#' 
#' The \code{plot} method is used for the side effect (and invisibly returns
#' \code{NULL}).
#'
#' @name methods.qfrm
#'
#' @seealso
#' \code{\link{new_qfrm}}: descriptions of the classes and their
#' \dQuote{constructors}
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
print.qfrm <- function(x, digits = getOption("digits"),
                       show_range = !is.null(x$error_bound), ...) {
    stat <- x$statistic
    errorb <- x$error_bound
    exact <- isTRUE(attr(errorb, "exact"))
    twosided <- isTRUE(attr(errorb, "twosided"))
    cat("\n")
    cat(strwrap("Moment of ratio of quadratic forms", prefix = "\t"),
        sep = "\n")
    cat("\n")
    out <- character()
    if(!is.null(stat)) {
        out <- c(out, paste("Moment =", format(stat, digits = max(1L, digits))))
    }
    # Is errorb all NA? (NaN should be excluded as it returns TRUE for is.na())
    all_na_errorb <- xor(all(is.na(errorb)), all(is.nan(errorb)))
    if(length(errorb) > 0 && !all_na_errorb && !exact) {
        out <- c(out,
                 paste("Error =", format(errorb, digits = max(1L, digits)),
                       if(twosided) " (two-sided)" else " (one-sided)"))
    }
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if(exact) {
        cat("This value is exact\n")
    } else if(length(errorb) > 0 && all_na_errorb) {
        cat("Error bound unavailable;",
            "recommended to inspect plot() of this object\n")
    } else if(show_range) {
        if(twosided) {
            ra <- sort(c(stat - errorb, stat + errorb))
        } else {
            ra <- sort(c(stat, stat + errorb))
        }
        cat("Possible range:\n ",
        paste(format(ra, digits = digits + 2L), collapse = " "), "\n", sep = "")
    }
    if(isTRUE(attr(errorb, "singular")) && !all_na_errorb) {
        cat("Note: Argument matrix numerically singular;",
            "error bound is unreliable\n")
    }
    if(isTRUE(attr(errorb, "alphaout")) && !all_na_errorb) {
        cat("Note: Adjustment parameter(s) alpha > 1;",
            "error bound is unreliable\n")
    }
    if(isTRUE(attr(stat, "diminished"))) {
        cat("Note: Numerical underflow encountered;",
            "result is likely inaccurate\n")
    }
    cat("\n")
    invisible(x)
}

#' @rdname methods.qfrm
#' @order 3
#'
#' @exportS3Method
#'
plot.qfrm <- function(x, add_error = length(x$seq_error) > 0,
                      add_legend = add_error,
                      ylim = x$statistic * ylim_f, ylim_f = c(0.9, 1.1),
                      xlab = "Order of evaluation", ylab = "Moment of ratio",
                      col_m = "royalblue4", col_e = "tomato",
                      lwd_m = 1, lwd_e = 1, lty_m = 1, lty_e = 2,
                      pos_leg = "topright", ...) {
    if(!requireNamespace("graphics", quietly = TRUE)) {
        message("Package 'graphics' not found. Exit from plot.qfrm")
        invisible()
    }
    terms <- x$terms
    seq_error <- x$seq_error
    seq_partial <- cumsum(terms)
    if(isTRUE(attr(seq_error, "exact"))) {
        message("plot method for this class is for inspecting ",
                "partial sums.\n  This particular object has an exact moment, ",
                "so plot is moot.")
    }
    try(plot(seq_along(terms) - 1L, seq_partial, type = "l", col = col_m,
         ylim = ylim, xlab = xlab,
         ylab = ylab, lwd = lwd_m, lty = lty_m, ...))
    if(add_error) {
        graphics::lines(seq_along(terms) - 1L, seq_error + seq_partial,
                        col = col_e, lwd = lwd_e, lty = lty_e)
        if(isTRUE(attr(seq_error, "twosided"))) {
            graphics::lines(-seq_error + seq_partial,
                            col = col_e, lwd = lwd_e, lty = lty_e)
        }
    }
    if(add_legend) {
        graphics::legend(pos_leg, legend = c("Partial sum", "Error bound"),
               col = c(col_m, col_e), lwd = c(lwd_m, lwd_e),
               lty = c(lty_m, lty_e))
    }
    invisible()
}

#' @rdname methods.qfrm
#' @order 2
#'
#' @exportS3Method
#'
print.qfpm <- function(x, digits = getOption("digits"), ...) {
    stat <- x$statistic
    errorb <- x$error_bound
    exact <- isTRUE(attr(errorb, "exact"))
    cat("\n")
    cat(strwrap("Moment of (product of) quadratic form(s)", prefix = "\t"),
        sep = "\n")
    cat("\n")
    out <- paste("Moment =", format(stat, digits = max(1L, digits)))
    cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
    if(exact) cat("This value is exact\n")
    cat("\n")
    invisible(x)
}
