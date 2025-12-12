#' Calculate hypergeometric series
#'
#' These internal functions calculate (summands of) hypergeometric series.
#'
#' The denominator parameter \code{b} is assumed positive,
#' whereas the numerator parameters can be positive or negative.  The signs
#' of the latter will be reflected in the result.
#'
#' @name hgs
#'
#' @param dks
#'   \code{(m + 1)} vector for \eqn{d_{i}},
#'   \code{(m + 1) * (m + 1)} square matrix for \eqn{d_{i,j}}, or
#'   \code{(m + 1) * (m + 1) * (m + 1)} array for \eqn{d_{i,j,k}}
#'   (\eqn{i, j, k = 0, 1, \dots m})
#' @param a1,a2,a3
#'   Numerator parameters
#' @param b
#'   Denominator parameter
#' @param lconst
#'   Scalar \eqn{\log c}{log c}
#'
#' @return
#' Numeric with the same dimension with \code{dks}
#'
NULL

# Internal function to calculate sequence of log rising factorial
.lrfseq <- function(a, m) {
    if(a <= 0 && (a %% 1) == 0) {
        return(cumsum(suppressWarnings(log(c(1, pmax(seq.int(-a, -a - (m - 1)), 0))))))
    } else {
        return(lgamma(seq.int(a, a + m)) - lgamma(a))
    }
}

#' Calculate 1D hypergeometric series
#'
#' \code{hgs_1d()} calculates the hypergeometric series
#' \eqn{c \frac{(a_1)_i}{(b)_i} d_{i}}{c ( (a_1)_i / (b)_i ) d_{i}}
#'
#' @rdname hgs
#'
hgs_1d <- function(dks, a1, b, lconst = 0) {
    m <- length(dks) - 1
    lnum_i <- .lrfseq(a1, m)
    lden_i <- .lrfseq(b, m)
    lcoefm <- lnum_i - lden_i
    ansseq <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (a1)_i, where i = 0:m is along the rows
    sgnseq <- cumprod(c(1, sign(seq.int(a1, a1 + m - 1)))) * sign(dks)
    ansseq * sgnseq
}

#' Calculate 2D hypergeometric series
#'
#' \code{hgs_2d()} calculates the hypergeometric series
#' \eqn{c \frac{(a_1)_i (a_2)_j}{(b)_{i+j}} d_{i, j}}{c
#'      ( (a_1)_i (a_2)_j / (b)_{i+j} ) d_{i,j}}
#'
#' @rdname hgs
#'
hgs_2d <- function(dks, a1, a2, b, lconst = 0) {
    m <- ncol(dks) - 1
    lnum_i <- .lrfseq(a1, m)
    lnum_j <- .lrfseq(a2, m)
    lden_ij <- lgamma(b + outer(0:m, 0:m, "+")) - lgamma(b)
    lcoefm <- (outer(lnum_i, lnum_j, "+") - lden_ij)
    ansmat <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (a1)_i and (a2)_j, where i = 0:m is along the rows
    sgnseq <- outer(cumprod(c(1, sign(seq.int(a1, a1 + m - 1)))),
                    cumprod(c(1, sign(seq.int(a2, a2 + m - 1))))) * sign(dks)
    ansmat * sgnseq
}

#' Calculate 3D hypergeometric series
#'
#' \code{hgs_3d()} calculates the hypergeometric series
#' \eqn{c \frac{(a_1)_i (a_2)_j (a_3)_k}{(b)_{i+j+k}} d_{i, j, k}}{c
#'      ( (a_1)_i (a_2)_j (a_3)_k / (b)_{i+j+k} ) d_{i,j,k}}
#'
#' @rdname hgs
#'
hgs_3d <- function(dks, a1, a2, a3, b, lconst = 0) {
    m <- dim(dks)[1] - 1
    lnum_i <- .lrfseq(a1, m)
    lnum_j <- .lrfseq(a2, m)
    lnum_k <- .lrfseq(a3, m)
    lden_ijk <- lgamma(b + outer(outer(0:m, 0:m, "+"), 0:m, "+")) -
                lgamma(b)
    lcoefm <- (outer(outer(lnum_i, lnum_j, "+"), lnum_k, "+") - lden_ijk)
    ansmat <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (a1)_i and (a2)_j, where i = 0:m is along the rows
    sgnseq <- outer(outer(cumprod(c(1, sign(seq.int(a1, a1 + m - 1)))),
                          cumprod(c(1, sign(seq.int(a2, a2 + m - 1))))),
                    cumprod(c(1, sign(seq.int(a3, a3 + m - 1))))) * sign(dks)
    ansmat * sgnseq
}
