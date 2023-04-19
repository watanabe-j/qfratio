##### pqfr_A1B1 #####
#' Distribution Function of Ratio of Quadratic Forms
#'
#' \sQuote{Exact} distribution function of the ratio of quadratic forms in
#' normal variables,
#' \eqn{\frac{ (\mathbf{x^\mathit{T} A x}) }{ (\mathbf{x^\mathit{T} B x}) }
#'      }{ (x^T A x) / (x^T B x) }, where
#' \eqn{\mathbf{x} \sim N_n(\bm{\mu}, \mathbf{\Sigma})}{x ~ N_n(\mu, \Sigma)},
#' based on Forchini (2002, 2005).
#'
#' This function is **experimental**.
#'
#' \code{quantile} should be of length one, because the cumulative distribution
#' is a complicated function of it and there is no practical way for
#' vectorization.
#'
#' Evaluates the (cumulative) distribution function as a partial sum of
#' infinite series involving top-order zonal or invariant polynomials.  The
#' expression for the central case is from Forchini (2002), and that for
#' the noncentral case is from Forchini (2005) (with apparent errors
#' corrected).  As in other functions of this package, this is evaluated with
#' the recursive algorithm for \eqn{d} in Hillier et al. (2009, 2014).
#'
#' The distribution function has points of nonanalyticity at the relative
#' eigenvalues of the argument matrices (Forchini 2002), and the series
#' expression tends to be *very slow* to converge around them.  In this case,
#' use a large \code{m}, or seek for alternative ways of evaluation.
#'
#' As a distribution function, the returned value is a \eqn{p}-value on the
#' lower tail; \eqn{P \left\{ \frac{ (\mathbf{x^\mathit{T} A x}) }{
#'                             (\mathbf{x^\mathit{T} B x}) } \le q \right\}
#' }{P\{(x^T A x) / (x^T B x) <= q\}}.  The reciprocal (like
#' \code{lower.tail = FALSE} in regular probability distribution functions)
#' is not supported.
#'
#' @inheritParams qfrm
#'
#' @param quantile
#'   Quantile \eqn{q}; numeric of length one
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}}{\Sigma} for
#'   \eqn{\mathbf{x}}{x}
#' @param check_convergence
#'   Specifies how numerical convergence is checked (see
#'   \code{\link{qfrm}})
#' @param cpp_method
#'   Method used in \proglang{C++} calculations to avoid numerical
#'   overflow/underflow (see \dQuote{Details} in \code{\link{qfrm}}).  Options:
#'   \itemize{
#'     \item{\code{"double"}: }{default; fastest but prone to underflow in
#'           some conditions}
#'     \item{\code{"long_double"}: }{same algorithm but using the
#'           \code{long double} variable type; robust but slow and
#'           memory-inefficient}
#'     \item{\code{"coef_wise"}: }{coefficient-wise scaling algorithm;
#'           experimental and may not work properly for this function}
#'    }
#' @param nthreads
#'   Number of threads used in \proglang{OpenMP}-enabled \proglang{C++}
#'   functions.  \code{0} or any negative value is special and means one-half of
#'   the number of processors detected.  See \dQuote{Multithreading} in
#'   \code{\link{qfrm}}.
#'
#' @return
#' List containing the following elements:
#' \itemize{
#'   \item{\code{$p}: }{\eqn{p}-value on the lower tail (\code{sum(terms)})}
#'   \item{\code{$terms}: }{vector of \eqn{0}th to \eqn{m}th order terms}
#'  }
#'
#' @references
#' Forchini, G. (2002) The exact cumulative distribution function of
#'   a ratio of quadratic forms in normal variables, with application to
#'   the AR(1) model. *Econometric Theory*, **18**, 823--852.
#'   \doi{10.1017/S0266466602184015}.
#'
#' Forchini, G. (2005) The distribution of a ratio of quadratic forms in
#'   noncentral normal variables.
#'   *Communications in Statistics---Theory and Methods*, **34**, 999--1008.
#'   \doi{10.1081/STA-200056855}.
#'
#' Hillier, G., Kan, R. and Wang, X. (2009) Computationally efficient recursions
#'   for top-order invariant polynomials with applications.
#'   *Econometric Theory*, **25**, 211--242.
#'   \doi{10.1017/S0266466608090075}.
#'
#' Hillier, G., Kan, R. and Wang, X. (2014) Generating functions and
#'   short recursions, with applications to the moments of quadratic forms
#'   in noncentral normal vectors. *Econometric Theory*, **30**, 436--473.
#'   \doi{10.1017/S0266466613000364}.
#'
#' @export
#'
#' @examples
#' ## Some symmetric matrices and parameters
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(sqrt(1:nv))
#' mu <- 0.2 * nv:1
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#'
#' ## P{ (x^T A x) / (x^T x) <= 2.5} where x ~ N(0, I)
#' pqfr_A1B1(2.5, A)$p
#'
#' ## P{ (x^T A x) / (x^T B x) <= 1.5} where x ~ N(0, I)
#' pqfr_A1B1(1.5, A, B)$p
#'
#' ## P{ (x^T A x) / (x^T B x) <= 1.5} where x ~ N(mu, I)
#' pqfr_A1B1(1.5, A, B, mu = mu)$p
#'
#' ## P{ (x^T A x) / (x^T B x) <= 1.5} where x ~ N(mu, Sigma)
#' pqfr_A1B1(1.5, A, B, mu = mu, Sigma = Sigma)$p
#'
pqfr_A1B1 <- function(quantile, A, B, m = 100L,
                      mu = rep.int(0, n),
                      Sigma = diag(n),
                      check_convergence = c("relative", "strict_relative",
                                            "absolute", "none"),
                      use_cpp = TRUE,
                      cpp_method = c("double", "long_double", "coef_wise"),
                      nthreads = 1,
                      tol_conv = .Machine$double.eps ^ (1/4),
                      tol_zero = .Machine$double.eps * 100,
                      tol_sing = .Machine$double.eps * 100,
                      thr_margin = 100) {
    if(isTRUE(check_convergence)) check_convergence <- "strict_relative"
    if(isFALSE(check_convergence)) check_convergence <- "none"
    check_convergence <- match.arg(check_convergence)
    if(!missing(cpp_method)) use_cpp <- TRUE
    cpp_method <- match.arg(cpp_method)
    ## If A or B is missing, let it be an identity matrix
    ## If they are given, symmetrize
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
    zeros <- rep.int(0, n)
    ## If Sigma is given, transform A, B, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma) && !iseq(Sigma, In, tol_zero)) {
        KiKS <- KiK(Sigma, tol_sing)
        K <- KiKS$K
        iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        KtBK <- t(K) %*% B %*% K
        iKmu <- iK %*% mu
        ## If Sigma is singular, check conditions for A, B, mu, and Sigma
        if(ncol(K) != n) {
            okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
                    (iseq(A %*% mu, zeros, tol_zero) &&
                     iseq(B %*% mu, zeros, tol_zero)) ||
                    (iseq(crossprod(iK, KtAK %*% iK), A) &&
                     iseq(crossprod(iK, KtBK %*% iK), B))
            if(!okay) {
                stop("For singular Sigma, certain condition must be met ",
                     "for A, B, mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details")
            }
        }
        return(pqfr_A1B1(quantile, KtAK, KtBK, m = m, mu = iKmu,
                         check_convergence = check_convergence,
                         use_cpp = use_cpp, cpp_method = cpp_method,
                         nthreads = nthreads,
                         tol_conv = tol_conv, tol_zero = tol_zero,
                         tol_sing = tol_sing, thr_margin = thr_margin))
    }
    ## Check basic requirements for arguments
    stopifnot(
        "A and B must be square matrices" = all(c(dim(A), dim(B)) == n),
        "quantile must be a length-one numeric" =
            is.numeric(quantile) && (length(quantile) == 1)
    )
    diminished <- FALSE
    eigA_qB <- eigen(A - quantile * B, symmetric = TRUE)
    Ds <- eigA_qB$values
    D1 <- Ds[Ds > tol_sing]
    D2 <- -Ds[Ds < -tol_sing]
    n1 <- length(D1)
    n2 <- length(D2)
    if(n2 == 0) return(list(p = 0, terms = rep.int(0, m + 1)))
    if(n1 == 0) return(list(p = 1, terms = c(1, rep.int(0, m))))
    mu <- c(crossprod(eigA_qB$vectors, c(mu)))
    mu1 <- mu[Ds > tol_sing]
    mu2 <- mu[Ds < -tol_sing]
    if(use_cpp) {
        if(cpp_method == "coef_wise") {
            cppres <- A1B1_Ec(D1, D2, mu1, mu2, m,
                              thr_margin, nthreads, tol_zero)
        } else if(cpp_method == "long_double") {
            cppres <- A1B1_El(D1, D2, mu1, mu2, m,
                              thr_margin, nthreads, tol_zero)
        } else {
            cppres <- A1B1_Ed(D1, D2, mu1, mu2, m,
                              thr_margin, nthreads, tol_zero)
        }
        ansseq <- cppres$ansseq
        diminished <- cppres$diminished
    } else {
        D1i <- 1 / D1
        D2i <- 1 / D2
        sumtrDi <- sum(D1i) + sum(D2i)
        D1d <- D1i / sumtrDi
        D2d <- D2i / sumtrDi
        central <- iseq(mu1, rep.int(0, n1), tol_zero) &&
                   iseq(mu2, rep.int(0, n2), tol_zero)
        D1h <- rep.int(sum(D1d), n1) - D1d
        D2h <- rep.int(sum(D2d), n2) - D2d
        seq0m <- seq.int(0, m)
        if(central) {
            dk1 <- d1_i(D1h, m, thr_margin)
            dk2 <- d1_i(D2h, m, thr_margin)
            lscf1 <- attr(dk1, "logscale")
            lscf2 <- attr(dk2, "logscale")
        } else {
            dkm1 <- d2_ij_m(tcrossprod(sqrt(D1d) * mu1) / 2, diag(D1h, n1), m,
                            thr_margin = thr_margin)
            dkm2 <- d2_ij_m(tcrossprod(sqrt(D2d) * mu2) / 2, diag(D2h, n2), m,
                            thr_margin = thr_margin)
            seqlrf_1_2 <- lgamma(1 / 2 + seq0m) - lgamma(1 / 2)
            dk1 <- sum_counterdiag(exp(log(dkm1) - seqlrf_1_2))
            dk2 <- sum_counterdiag(exp(log(dkm2) - seqlrf_1_2))
            lscf1 <- attr(dkm1, "logscale")[1, ]
            lscf2 <- attr(dkm2, "logscale")[1, ]
            diminished <- any(lscf1 < 0) && any(diag(dkm1[(m + 1):1, ]) == 0) ||
                          any(lscf2 < 0) && any(diag(dkm2[(m + 1):1, ]) == 0)
        }
        ordmat <- outer(seq0m, seq0m, FUN = "+") + (n1 + n2) / 2
        ansmat <- lgamma(ordmat)
        ansmat <- ansmat - lgamma(n1 / 2 + 1 + seq0m) + log(dk1)
        ansmat <- t(t(ansmat) - lgamma(n2 / 2 + seq0m) + log(dk2))
        ansmat <- ansmat + (sum(log(D1d)) + sum(log(D2d))) / 2
        ansmat <- ansmat - lscf1
        ansmat <- t(t(ansmat) - lscf2)
        ansmat <- ansmat - (sum(mu1 ^ 2) + sum(mu2 ^ 2)) / 2
        ansmat <- exp(ansmat)
        ansmat <- ansmat * gsl::hyperg_2F1(ordmat, 1, n1 / 2 + 1 + seq0m,
                                           sum(D1d))
        ansseq <- sum_counterdiag(ansmat)
    }
    if(diminished) {
        warning("Some terms in multiple series numerically diminished to 0 ",
                "as\n  scaled to avoid numerical overflow. ",
                "Result will be inaccurate",
                if(cpp_method != "coef_wise")
                    paste0(".\n  Consider using cpp_method = ",
                          if(cpp_method != "long_double") "\"long_double\" or ",
                          "\"coef_wise\"."))
    }
    if(check_convergence != "none") {
        if(check_convergence == "strict_relative") {
            if(tol_conv > .Machine$double.eps) tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term is >",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the series",
                    ",\n  suggesting non-convergence. Consider using larger m")
        }
    }
    list(p = sum(ansseq), terms = ansseq)
}
