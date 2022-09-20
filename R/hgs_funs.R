#' Calculate hypergeometric series
#'
#' These internal functions calculate (summands of) hypergeometric series.
#'
#' The denominator paramter \code{dpi} (or \eqn{q}) is assumed positive,
#' whereas the numerator parameters can be positive or negative.
#' The signs of the latter will be reflected in the result.
#'
#' @name hgs
#'
#' @param dks
#'   \eqn{(m + 1)} vector for \eqn{d_{i}},
#'   \eqn{(m + 1) \times (m + 1)} square matrix for \eqn{d_{i, j}}, or
#'   \eqn{(m + 1) \times (m + 1) \times (m + 1)} array for \eqn{d_{i, j, k}}
#'   (\eqn{i, j, k = 0, 1, \dots m})
#' @param np1,np2,np3
#'   Numerator parameters. Scalars \eqn{p_1}, \eqn{p_2}, and \eqn{p_3}
#' @param dpi,dpij,dpijk
#'   Denominator parameter. Scalar \eqn{q}
#' @param lconst
#'   Scalar \eqn{\log c}
#'
#' @return
#' Numeric with the same dimension with \code{dks}
#'
NULL

#' Calculate 1D hypergeometric series
#'
#' \code{hgs_1d()} calculates the hypergeometric series
#' \eqn{c \frac{(p_1)_i}{(q)_i} d_{i}}
#'
#' @rdname hgs
#'
hgs_1d <- function (dks, np1, dpi, lconst = 0) {
    m <- length(dks) - 1
    if(np1 < 0 && (np1 %% 1) == 0) {
        lnum_i <- cumsum(log(c(1, -np1 - 0:(m - 1))))
    } else {
        lnum_i <- lgamma(np1 + 0:m) - lgamma(np1)
    }
    lden_i <- lgamma(dpi + 0:m) - lgamma(dpi)
    lcoefm <- lnum_i - lden_i
    ansseq <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (np1)_i, where i = 0:m is along the rows
    sgnseq <- cumprod(c(1, sign(np1 + 0:(m - 1)))) * sign(dks)
    ansseq * sgnseq
}

#' Calculate 2D hypergeometric series
#'
#' \code{hgs_2d()} calculates the hypergeometric series
#' \eqn{c \frac{(p_1)_i (p_2)_j}{(q)_{i+j}} d_{i, j}}
#'
#' @rdname hgs
#'
hgs_2d <- function (dks, np1, np2, dpij, lconst = 0) {
    m <- ncol(dks) - 1
    if(np1 < 0 && (np1 %% 1) == 0) {
        lnum_i <- cumsum(log(c(1, -np1 - 0:(m - 1))))
    } else {
        lnum_i <- lgamma(np1 + 0:m) - lgamma(np1)
    }
    if(np2 < 0 && (np2 %% 1) == 0) {
        lnum_j <- cumsum(log(c(1, -np2 - 0:(m - 1))))
    } else {
        lnum_j <- lgamma(np2 + 0:m) - lgamma(np2)
    }
    lden_ij <- lgamma(dpij + outer(0:m, 0:m, "+")) - lgamma(dpij)
    lcoefm <- (outer(lnum_i, lnum_j, "+") - lden_ij)
    ansmat <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (np1)_i and (np2)_j, where i = 0:m is along the rows
    sgnseq <- outer(cumprod(c(1, sign(np1 + 0:(m - 1)))),
                    cumprod(c(1, sign(np2 + 0:(m - 1))))) * sign(dks)
    ansmat * sgnseq
}

#' Calculate 3D hypergeometric series
#'
#' \code{hgs_3d()} calculates the hypergeometric series
#' \eqn{c \frac{(p_1)_i (p_2)_j (p_3)_k}{(q)_{i+j+k}} d_{i, j, k}}
#'
#' @rdname hgs
#'
hgs_3d <- function (dks, np1, np2, np3, dpijk, lconst = 0) {
    m <- dim(dks)[1] - 1
    if(np1 < 0 && (np1 %% 1) == 0) {
        lnum_i <- cumsum(log(c(1, -np1 - 0:(m - 1))))
    } else {
        lnum_i <- lgamma(np1 + 0:m) - lgamma(np1)
    }
    if(np2 < 0 && (np2 %% 1) == 0) {
        lnum_j <- cumsum(log(c(1, -np2 - 0:(m - 1))))
    } else {
        lnum_j <- lgamma(np2 + 0:m) - lgamma(np2)
    }
    if(np3 < 0 && (np3 %% 1) == 0) {
        lnum_k <- cumsum(log(c(1, -np3 - 0:(m - 1))))
    } else {
        lnum_k <- lgamma(np3 + 0:m) - lgamma(np3)
    }
    lden_ijk <- lgamma(dpijk + outer(outer(0:m, 0:m, "+"), 0:m, "+")) -
                lgamma(dpijk)
    lcoefm <- (outer(outer(lnum_i, lnum_j, "+"), lnum_k, "+") - lden_ijk)
    ansmat <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (np1)_i and (np2)_j, where i = 0:m is along the rows
    sgnseq <- outer(outer(cumprod(c(1, sign(np1 + 0:(m - 1)))),
                          cumprod(c(1, sign(np2 + 0:(m - 1))))),
                    cumprod(c(1, sign(np3 + 0:(m - 1))))) * sign(dks)
    ansmat * sgnseq
}

#' Calculate 2D hypergeometric series
#'
#' \code{hgs_2d()} calculates the hypergeometric series
#' \eqn{c \frac{(p_1)_i (p_2)_j}{(q)_{i+j}} d_{i, j}}
#'
#' @rdname hgs
#'
hgs_dmu_2d <- function (dks, np1, np2, dpij, lconst = 0) {
    m <- ncol(dks) - 1
    if(np1 < 0 && (np1 %% 1) == 0) {
        lnum_i <- cumsum(log(c(1, -np1 - 0:(m - 1))))
    } else {
        lnum_i <- lgamma(np1 + 0:m) - lgamma(np1)
    }
    if(np2 < 0 && (np2 %% 1) == 0) {
        lnum_j <- cumsum(log(c(1, -np2 - 0:(m - 1))))
    } else {
        lnum_j <- lgamma(np2 + 0:m) - lgamma(np2)
    }
    lden_ij <- lgamma(dpij + outer(0:m, 0:m, "+")) - lgamma(dpij)
    lden_ij <- t(t(lden_ij) + 0:m * log(2) + lgamma(1 / 2 + 0:m) - lgamma(1 / 2))
    lcoefm <- (outer(lnum_i, lnum_j, "+") - lden_ij)
    ansmat <- exp(lcoefm + log(abs(dks)) + lconst)
    ## Signs for (np1)_i and (np2)_j, where i = 0:m is along the rows
    sgnseq <- outer(cumprod(c(1, sign(np1 + 0:(m - 1)))),
                    cumprod(c(1, sign(np2 + 0:(m - 1))))) * sign(dks)
    ansmat * sgnseq
}
