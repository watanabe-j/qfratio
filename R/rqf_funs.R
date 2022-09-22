##### rqfr #####
#' Monte Carlo sampling of ratio/product of quadratic forms
#'
#' \code{rqfr()}, \code{rqfmr()}, \code{rqfp()} calculates a random sample of
#' a simple ratio, multiple ratio (of special form), and product, respectively,
#' of quadratic forms in normal variables of specified mean and covariance
#' (standard multivariate normal by default).
#' These functions are primarily for empirical verification of the analytic
#' results provided in this package.
#'
#' These functions generate a random smaple of
#' \eqn{ \frac{(\mathbf{x^\mathit{T} A x})^p }{(\mathbf{x^\mathit{T} B x})^q} } (\code{rqfr()}),
#' \eqn{ \frac{(\mathbf{x^\mathit{T} A x})^p}{(\mathbf{x^\mathit{T} B x})^q (\mathbf{x^\mathit{T} Dx})^r} } (\code{rqfmr()}), and
#' \eqn{ (\mathbf{x^\mathit{T} A x})^p (\mathbf{x^\mathit{T} B x})^q (\mathbf{x^\mathit{T} D x})^r } (\code{rqfp()}),
#' where \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{\Sigma})}.
#'
#' When only one of \code{p} and \code{q} are provided in \code{rqfr()},
#' the other (missing) one is set to the same value.
#'
#' In \code{rqfmr()}, \code{q} and \code{r} are set to \code{p/2}
#' when both missing, and set to the same value when only one is missing.
#' When \code{p} is missing, this is set to be \code{q + r}.
#' If unsure, specify all these explicitly.
#'
#' In \code{rqfp()}, \code{p}, \code{q} and \code{r} are \code{1} by default.
#'
#' @param nit
#'   Number of iteration or sample size.  Should be an integer-alike of
#'   length 1.
#' @param A,B,D
#'   Argument matrices (see "Details").  Assumed to be square matrices of
#'   the same order.  When missing, set to the identity matrix.  At least
#'   one of these should be specified.
#' @param p,q,r
#'   Exponents for A, B, D, respectively (see "Details").  Assumed to be
#'   numeric of length 1 each.  See "Details" for default values.
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}. Default zero vector.
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}} for \eqn{\mathbf{x}}.
#'   Default identity matrix.
#'   \code{mu} and \code{Sigma} are assumed to be of the same order as
#'   the argument matrices.
#' @param use_cpp
#'   Logical to specify whether an \code{C++} version is called or not.
#'   Default \code{FALSE}.
#'
#' @return Numeric vector of length \code{nit}.
#'
#' @seealso \code{\link{qfrm}} and \code{\link{qfpm}} for analytic moments
#'
#' @name rqfr
#'
#' @export
#'
#' @examples
#' p <- 4
#' A <- diag(1:p)
#' B <- diag(p:1)
#' D <- diag(sqrt(1:p))
#' rqfr(5, A) ## By default B = I, p = q = 1, x ~ N(0, I)
#'
#' rqfmr(5, A, B, D, 1, 1/2, 1/2)
#'
#' rqfp(5, A, B, D, 1, 1, 1)
#'
#' ## Example with non-standard normal
#' mu <- 1:p / p
#' Sigma <- matrix(0.5, p, p)
#' diag(Sigma) <- 1
#' rqfr(5, A, mu = 1:p / p, Sigma = Sigma)
#'
#' ## Compare Monte Carlo mean and analytic expression
#' mcres <- rqfr(1000, A, p = 2)
#' mean(mcres)
#' sd(mcres)
#' qfrm(A, p = 2)
#'
rqfr <- function(nit = 1000L, A, B, p = 1, q = p, mu, Sigma, use_cpp = FALSE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package \"mvtnorm\" is required to use this function.")
    }
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
        A <- (A + t(A)) / 2
    }
    if(missing(B)) {
        B <- In
    } else {
        B <- (B + t(B)) / 2
    }
    if(missing(p) && !missing(q)) p <- q
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    # X <- eigvaldisp::rmvn(nit, mean = mu, Sigma = Sigma)
    if(use_cpp) {
        return(rqfrE(nit, A, B, p, q, mu, Sigma))
    }
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    diag(tcrossprod(tcrossprod(X, A), X)) ^ p / diag(tcrossprod(tcrossprod(X, B), X)) ^ q
}

##### rqfmr #####
#' @rdname rqfr
#' @export
#'
rqfmr <- function(nit = 1000L, A, B, D, p = 1, q = p / 2, r = q, mu, Sigma, use_cpp = FALSE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package \"mvtnorm\" is required to use this function.")
    }
    if(missing(A)) {
        if(missing(B)) {
            if(missing(D)) {
                stop("Provide at least one of A, B and D")
            } else {
                n <- dim(D)[1L]
            }
        } else {
            n <- dim(B)[1L]
        }
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
        A <- (A + t(A)) / 2
    }
    if(missing(B)) {
        B <- In
    } else {
        B <- (B + t(B)) / 2
    }
    if(missing(D)) {
        D <- In
    } else {
        D <- (D + t(D)) / 2
    }
    if(missing(q) && !missing(r)) q <- r
    if(missing(p) && !missing(q)) p <- q + r
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    if(use_cpp) {
        return(rqfmrE(nit, A, B, D, p, q, r, mu, Sigma))
    }
    # X <- eigvaldisp::rmvn(nit, mean = mu, Sigma = Sigma)
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    diag(tcrossprod(tcrossprod(X, A), X)) ^ p / diag(tcrossprod(tcrossprod(X, B), X)) ^ q / diag(tcrossprod(tcrossprod(X, D), X)) ^ r
}

##### rqfp #####
#' @rdname rqfr
#' @export
#'
rqfp <- function(nit = 1000L, A, B, D, p = 1, q = 1, r = 1, mu, Sigma, use_cpp = FALSE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package \"mvtnorm\" is required to use this function.")
    }
    if(missing(A)) {
        if(missing(B)) {
            if(missing(D)) {
                stop("Provide at least one of A, B and D")
            } else {
                n <- dim(D)[1L]
            }
        } else {
            n <- dim(B)[1L]
        }
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
        A <- (A + t(A)) / 2
    }
    if(missing(B)) {
        B <- In
    } else {
        B <- (B + t(B)) / 2
    }
    if(missing(D)) {
        D <- In
    } else {
        D <- (D + t(D)) / 2
    }
    # if(missing(q) && !missing(r)) q <- r
    # if(missing(p) && !missing(q)) p <- q
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    if(use_cpp) {
        return(rqfpE(nit, A, B, D, p, q, r, mu, Sigma))
    }
    # X <- eigvaldisp::rmvn(nit, mean = mu, Sigma = Sigma)
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    diag(tcrossprod(tcrossprod(X, A), X)) ^ p * diag(tcrossprod(tcrossprod(X, B), X)) ^ q * diag(tcrossprod(tcrossprod(X, D), X)) ^ r
}
