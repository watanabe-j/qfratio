##### rqfr #####
#' Monte Carlo sampling of ratio/product of quadratic forms
#'
#' \code{rqfr()}, \code{rqfmr()}, and \code{rqfp()} calculate a random sample of
#' a simple ratio, multiple ratio (of special form), and product, respectively,
#' of quadratic forms in normal variables of specified mean and covariance
#' (standard multivariate normal by default).  These functions are primarily for
#' empirical verification of the analytic results provided in this package.
#'
#' These functions generate a random sample of
#' \eqn{ \frac{(\mathbf{x^\mathit{T} A x})^p}{(\mathbf{x^\mathit{T} B x})^q}
#'      }{(x^T A x)^p / (x^T B x)^q}
#' (\code{rqfr()}),
#' \eqn{ \frac{(\mathbf{x^\mathit{T} A x})^p}
#'            {(\mathbf{x^\mathit{T} B x})^q (\mathbf{x^\mathit{T} Dx})^r}
#'      }{(x^T A x)^p / ( (x^T B x)^q (x^T D x)^r )}
#' (\code{rqfmr()}), and
#' \eqn{ (\mathbf{x^\mathit{T} A x})^p (\mathbf{x^\mathit{T} B x})^q
#'       (\mathbf{x^\mathit{T} D x})^r }{(x^T A x)^p (x^T B x)^q (x^T D x)^r}
#' (\code{rqfp()}), where
#' \eqn{\mathbf{x} \sim N_n(\bm{\mu}, \mathbf{\Sigma})
#'      }{x ~ N_n(\mu, \Sigma)}.  (Internally, \code{rqfr()} and \code{rqfmr()}
#' just call \code{rqfp()} with negative exponents.)
#'
#' When only one of \code{p} and \code{q} is provided in \code{rqfr()},
#' the other (missing) one is set to the same value.
#'
#' In \code{rqfmr()}, \code{q} and \code{r} are set to \code{p/2}
#' when both missing, and set to the same value when only one is missing.  When
#' \code{p} is missing, this is set to be \code{q + r}.  If unsure,
#' specify all these explicitly.
#'
#' In \code{rqfp()}, \code{p}, \code{q} and \code{r} are \code{1} by default,
#' provided that the corresponding argument matrices are given.  If both
#' an argument matrix and its exponent (e.g., \code{D} and \code{r})
#' are missing, the exponent is set to \code{0} so that the factor be unity.
#'
#' @param nit
#'   Number of iteration or sample size.  Should be an integer-alike of
#'   length 1.
#' @param A,B,D
#'   Argument matrices (see \dQuote{Details}).  Assumed to be square matrices of
#'   the same order.  When missing, set to the identity matrix.  At least
#'   one of these must be specified.
#' @param p,q,r
#'   Exponents for A, B, D, respectively (see \dQuote{Details}).  Assumed to be
#'   numeric of length 1 each.  See \dQuote{Details} for default values.
#' @param mu
#'   Mean vector \eqn{\bm{\mu}}{\mu} for \eqn{\mathbf{x}}{x}.  Default
#'   zero vector.
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}}{\Sigma} for
#'   \eqn{\mathbf{x}}{x}.  Default identity matrix.  \code{mu} and
#'   \code{Sigma} are assumed to be of the same order as the argument matrices.
#' @param use_cpp
#'   Logical to specify whether an \proglang{C++} version is called or
#'   not.  \code{TRUE} by default.
#'
#' @return Numeric vector of length \code{nit}.
#'
#' @seealso \code{\link{qfrm}} and \code{\link{qfpm}} for analytic moments
#'
#' \code{\link{dqfr}} for analytic distribution-related functions for
#' simple ratios
#'
#' @name rqfr
#'
#' @examples
#' p <- 4
#' A <- diag(1:p)
#' B <- diag(p:1)
#' D <- diag(sqrt(1:p))
#'
#' ## By default B = I, p = q = 1;
#' ## i.e., (x^T A x) / (x^T x), x ~ N(0, I)
#' rqfr(5, A)
#'
#' ## (x^T A x) / ((x^T B x)(x^T D x))^(1/2), x ~ N(0, I)
#' rqfmr(5, A, B, D, 1, 1/2, 1/2)
#'
#' ## (x^T A x), x ~ N(0, I)
#' rqfp(5, A)
#'
#' ## (x^T A x) (x^T B x), x ~ N(0, I)
#' rqfp(5, A, B)
#'
#' ## (x^T A x) (x^T B x) (x^T D x), x ~ N(0, I)
#' rqfp(5, A, B, D)
#'
#' ## Example with non-standard normal
#' mu <- 1:p / p
#' Sigma <- matrix(0.5, p, p)
#' diag(Sigma) <- 1
#' rqfr(5, A, mu = 1:p / p, Sigma = Sigma)
#'
#' ## Compare Monte Carlo sample and analytic expression
#' set.seed(3)
#' mcres <- rqfr(1000, A, p = 2)
#' mean(mcres)
#' (anres <- qfrm(A, p = 2))
#' stats::t.test(mcres, mu = anres$statistic)
#'
#' @export
#'
rqfr <- function(nit = 1000L, A, B, p = 1, q = p, mu, Sigma, use_cpp = TRUE) {
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
    rqfp(nit, A, B, p = p, q = -q, r = 0,
         mu = mu, Sigma = Sigma, use_cpp = use_cpp)
}

##### rqfmr #####
#' @rdname rqfr
#' @export
#'
rqfmr <- function(nit = 1000L, A, B, D, p = 1, q = p / 2, r = q,
                  mu, Sigma, use_cpp = TRUE) {
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
    rqfp(nit, A, B, D, p = p, q = -q, r = -r,
         mu = mu, Sigma = Sigma, use_cpp = use_cpp)
}

##### rqfp #####
#' @rdname rqfr
#' @export
#'
rqfp <- function(nit = 1000L, A, B, D, p = 1, q = 1, r = 1,
                 mu, Sigma, use_cpp = TRUE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package 'mvtnorm' required to use this function")
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
        if(missing(p)) p <- 0
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
        A <- (A + t(A)) / 2
    }
    if(missing(B)) {
        B <- In
        if(missing(q)) q <- 0
    } else {
        B <- (B + t(B)) / 2
    }
    if(missing(D)) {
        D <- In
        if(missing(r)) r <- 0
    } else {
        D <- (D + t(D)) / 2
    }
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    if(use_cpp) {
        return(rqfpE(nit, A, B, D, p, q, r, mu, Sigma))
    }
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    if(p == 0) {
        qfAp <- rep.int(1, nit)
    } else {
        qfAp <- diag(tcrossprod(tcrossprod(X, A), X)) ^ p
    }
    if(q == 0) {
        qfBq <- rep.int(1, nit)
    } else {
        qfBq <- diag(tcrossprod(tcrossprod(X, B), X)) ^ q
    }
    if(r == 0) {
        qfDr <- rep.int(1, nit)
    } else {
        qfDr <- diag(tcrossprod(tcrossprod(X, D), X)) ^ r
    }
    qfAp * qfBq * qfDr
}
