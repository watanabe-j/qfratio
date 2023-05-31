##### pqfr #####
#' Probability distribution of ratio of quadratic forms
#'
#' This functionality is **experimental**.
#'
#' The user is supposed to use the exported functions \code{dqfr()} and
#' \code{pqfr()}, which are (pseudo-)vectorized with respect to \code{quantile}
#' using \code{\link[base]{sapply}()}.  The actual calculations are done by the
#' *internal* functions \code{dqfr_A1I1()}, \code{pqfr_A1B1()},
#' \code{pqfr_imhof()}, and \code{pqfr_davies()} which only
#' accommodate a length-one \code{quantile}.  The internal functions skip
#' most checks on argument structures and do not accommodate \code{Sigma}
#' (where applicable) to reduce execution time.
#'
#' \code{dqfr_A1I1()} and \code{pqfr_A1B1()} evaluate the probability density
#' and (cumulative) distribution function, respectively,
#' as a partial sum of infinite series involving top-order zonal or
#' invariant polynomials.  The density is from Hillier
#' (2001, corollary 2), the distribution function for the central case is
#' from Forchini (2002), and that for the noncentral case is from Forchini
#' (2005) (with apparent errors corrected).  As in other functions of
#' this package, these are evaluated with the recursive algorithm for
#' \eqn{d} in Hillier et al. (2009, 2014).
#'
#' \code{pqfr_imhof()} and \code{pqfr_davies()} evaluate the distribution
#' function by numerical inversion of the characteristic function based on
#' Imhof (1961) or Davies (1973, 1980), calling
#' \code{\link[CompQuadForm]{imhof}()} and \code{\link[CompQuadForm]{davies}()},
#' respectively, from the package \pkg{CompQuadForm}.  This is from the fact
#' that, for nonnegative definite \eqn{\mathbf{B}}{B},
#' \eqn{ P \left\{
#'         \frac{ \mathbf{x^\mathit{T} A x} }{ \mathbf{x^\mathit{T} B x} }
#'         > q \right\} =
#'       P \left\{ \mathbf{x}^\mathit{T} (\mathbf{A} - q \mathbf{B}) \mathbf{x}
#'         > 0 \right\} }{ P\{(x^T A x) / (x^T B x) > q\} =
#'                         P\{x^T (A - q B) x > 0\}}, so that the distribution
#' function of the ratio can be expressed in terms of that of the latter
#' (indefinite) quadratic form.  (The same argument was used in the derivation
#' of Forchini (2002, 2005) as well.)  Optional arguments for
#' \code{\link[CompQuadForm]{imhof}()} and \code{\link[CompQuadForm]{davies}()}
#' can be passed via \code{...}, except for \code{sigma} in the latter,
#' which is not applicable.
#'
#' The density is undefined at the eigenvalues of \eqn{\mathbf{A}}{A}
#' (Hillier 2001) and the distribution function has points of nonanalyticity
#' at the relative eigenvalues of the argument matrices (Forchini 2002).  Around
#' these points, the series expression tends to be *very slow* to converge.  In
#' this case, use a large \code{m}, avoid using \code{method = "forchini"} in
#' \code{pqfr()}, or seek alternative ways of evaluation.
#'
#' @inheritParams qfrm
#'
#' @param quantile
#'   Numeric vector of quantiles \eqn{q}
#' @param A,B
#'   Argument matrices.  Should be square.  \code{B} should be nonnegative
#'   definite.  Will be automatically symmetrized in \code{dqfr()} and
#'   \code{pqfr()}.
#' @param LA
#'   Eigenvalues of \eqn{\mathbf{A}}{|}
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}}{\Sigma} for
#'   \eqn{\mathbf{x}}{x}
#' @param log,lower.tail,log.p
#'   Logical; as in regular probability distribution functions.  But these are
#'   for convenience only, and not meant for accuracy.
#' @param method
#'   Method to specify an internal function in \code{pqfr()}
#'   (see \dQuote{Details}).  Options:
#'   \itemize{
#'     \item{\code{"imhof"}: }{default; uses \code{pqfr_imhof()}, numerical
#'           inversion of Imhof (1961)}
#'     \item{\code{"davies"}: }{uses \code{pqfr_davies()}, numerical inversion
#'           of Davies (1973, 1980)}
#'     \item{\code{"forchini"}: }{uses \code{pqfr_A1B1()}, series expression
#'           of Forchini (2002, 2005)}
#'    }
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
#'           experimental but supposedly robust}
#'    }
#' @param nthreads
#'   Number of threads used in \proglang{OpenMP}-enabled \proglang{C++}
#'   functions.  \code{0} or any negative value is special and means one-half of
#'   the number of processors detected.  See \dQuote{Multithreading} in
#'   \code{\link{qfrm}}.
#' @param ...
#'   Additional arguments passed to internal functions
#'
#' @return
#' \code{dqfr()} gives the density, and \code{pqfr()} gives the
#' distribution function or \eqn{p}-values corresponding to \code{quantile}.
#'
#' \code{dqfr_A1I1()} and \code{pqfr_A1B1()} return a list containing
#' the following elements:
#' \itemize{
#'   \item{\code{$d} or \code{$p}: }{probability density or
#'         lower \eqn{p}-value, respectively (\code{sum(terms)})}
#'   \item{\code{$terms}: }{vector of \eqn{0}th to \eqn{m}th order terms}
#'  }
#' Only \code{$d} or \code{$p} is passed to the external functions.
#'
#' \code{pqfr_imhof()} and \code{pqfr_davies()} simply return the full output
#' (a list) of \code{CompQuadForm::\link[CompQuadForm]{imhof}()} and
#' \code{CompQuadForm::\link[CompQuadForm]{davies}()}, respectively.  This may
#' be useful for debugging purposes.  Note that
#' \code{$Qq} there is the *upper* \eqn{p}-value.  When
#' \code{lower.tail = TRUE}, \code{pqfr()} takes its complement to return
#' the lower \eqn{p}-value.
#'
#' @references
#' Davis, R. B. (1973) Numerical inversion of a characteristic function.
#'   *Biometrika*, **60**, 415--417.
#'   \doi{10.1093/biomet/60.2.415}.
#'
#' Davis, R. B. (1980) Algorithm AS 155: The distribution of a linear
#'   combination of \eqn{\chi^2} random variables.
#'   *Journal of the Royal Statistical Society, Series C---Applied Statistics*,
#'   **29**, 323--333.
#'   \doi{10.2307/2346911}.
#'
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
#' Hillier, G. (2001) The density of quadratic form in a vector uniformly
#'   distributed on the \eqn{n}-sphere.
#'   *Econometric Theory*, **17**, 1--28.
#'   \doi{10.1017/S026646660117101X}.
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
#' Imhof, J. P. (1961) Computing the distribution of quadratic forms in normal
#'   variables.
#'   *Biometrika*, **48**, 419--426.
#'   \doi{10.1093/biomet/48.3-4.419}.
#'
#' @seealso \code{\link{rqfr}}, a Monte Carlo random number generator
#'
#' @name pqfr
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
#' ## f(2.5) for (x^T A x) / (x^T x) where x ~ N(0, I)
#' dqfr(2.5, A)
#'
#' ## P{ (x^T A x) / (x^T x) <= 2.5} where x ~ N(0, I)
#' pqfr(2.5, A)
#'
#' ## P{ (x^T A x) / (x^T B x) <= 1.5} where x ~ N(0, I)
#' pqfr(1.5, A, B)
#'
#' ## P{ (x^T A x) / (x^T B x) <= 1.5} where x ~ N(mu, I)
#' pqfr(1.5, A, B, mu = mu)
#'
#' ## P{ (x^T A x) / (x^T B x) <= 1.5} where x ~ N(mu, Sigma)
#' pqfr(1.5, A, B, mu = mu, Sigma = Sigma)
#'
#' ## These are (pseudo-)vectorized
#' qs <- 0:nv + 0.5
#' dqfr(qs, A)
#' pqfr(qs, A, B, mu = mu, Sigma = Sigma)
#'
NULL

##### pqfr #####
#' Probability distribution of ratio of quadratic forms
#'
#' \code{pqfr()}: \sQuote{Exact} distribution function of the ratio of
#' quadratic forms,
#' \eqn{\frac{ \mathbf{x^\mathit{T} A x} }{ \mathbf{x^\mathit{T} B x} }
#'      }{ (x^T A x) / (x^T B x) }, where
#' \eqn{\mathbf{x} \sim N_n(\bm{\mu}, \mathbf{\Sigma})}{x ~ N_n(\mu, \Sigma)}.
#'
#' @rdname pqfr
#' @order 2
#'
#' @export
#'
pqfr <- function(quantile, A, B, m = 100L, mu = rep.int(0, n), Sigma = diag(n),
                 lower.tail = TRUE, log.p = FALSE,
                 method = c("imhof", "davies", "forchini"),
                 tol_zero = .Machine$double.eps * 100,
                 tol_sing = .Machine$double.eps * 100,
                 ...) {
    method <- match.arg(method)
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
                     "for A, B, mu.\n  See documentation for details")
            }
        }
        return(pqfr(quantile, KtAK, KtBK, m = m, mu = iKmu,
                    lower.tail = lower.tail, log.p = log.p, method = method,
                    tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    LB <- eigen(B, symmetric = TRUE)$values
    ## Check basic requirements for arguments
    stopifnot(
        "A and B must be square matrices" = all(c(dim(A), dim(B)) == n),
        "B must be nonnegative definite" = all(LB >= -tol_sing),
        "quantile must be numeric" = is.numeric(quantile)
    )
    if(method == "forchini") {
        ans <- sapply(quantile,
                      function(q) pqfr_A1B1(q, A, B, m, mu = mu,
                                            tol_zero = tol_zero,
                                            tol_sing = tol_sing, ...)$p)
        if(!lower.tail) ans <- 1 - ans
    } else {
        if(method == "davies") {
            ans <- sapply(quantile,
                          function(q) pqfr_davies(q, A, B, mu = mu,
                                                  tol_zero = tol_zero, ...)$Qq)
        } else {
            ans <- sapply(quantile,
                          function(q) pqfr_imhof(q, A, B, mu = mu,
                                                 tol_zero = tol_zero, ...)$Qq)
        }
        if(lower.tail) ans <- 1 - ans
    }
    if(log.p) ans <- log(ans)
    attributes(ans) <- attributes(quantile)
    return(ans)
}

##### pqfr_A1B1 #####
#' Probability distribution of ratio of quadratic forms
#'
#' \code{pqfr_A1B1()}: internal function used within \code{pqfr()},
#' based on Forchini (2002, 2005).
#'
#' @rdname pqfr
#' @order 4
#'
pqfr_A1B1 <- function(quantile, A, B, m = 100L,
                      mu = rep.int(0, n),
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
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        A <- diag(n)
    } else {
        n <- dim(A)[1L]
        if(missing(B)) B <- diag(n)
    }
    ## Check basic requirements for arguments
    stopifnot(
        "In pqfr_A1B1, quantile must be length-one" = (length(quantile) == 1)
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
            cppres <- p_A1B1_Ec(D1, D2, mu1, mu2, m,
                              thr_margin, nthreads, tol_zero)
        } else if(cpp_method == "long_double") {
            cppres <- p_A1B1_El(D1, D2, mu1, mu2, m,
                              thr_margin, nthreads, tol_zero)
        } else {
            cppres <- p_A1B1_Ed(D1, D2, mu1, mu2, m,
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
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  Result truncated before first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        attr(ansseq, "truncated") <- TRUE
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

##### pqfr_imhof #####
#' Probability distribution of ratio of quadratic forms
#'
#' \code{pqfr_imhof()}: internal function used within \code{pqfr()},
#' based on Imhof (1961).
#'
#' @rdname pqfr
#' @order 5
#'
pqfr_imhof <- function(quantile, A, B, mu = rep.int(0, n),
                       tol_zero = .Machine$double.eps * 100, ...) {
    if(!requireNamespace("CompQuadForm", quietly = TRUE)) {
        stop("Package 'CompQuadForm' is required to use this function")
    }
    ## If A or B is missing, let it be an identity matrix
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        A <- diag(n)
    } else {
        n <- dim(A)[1L]
        if(missing(B)) B <- diag(n)
    }
    ## Check basic requirements for arguments
    stopifnot(
        "In pqfr_imhof, quantile must be length-one" = (length(quantile) == 1)
    )
    eigA_qB <- eigen(A - quantile * B, symmetric = TRUE)
    lambda <- eigA_qB$values
    delta2 <- c(crossprod(eigA_qB$vectors, mu)) ^ 2
    if(all(lambda <=  tol_zero)) return(list(Qq = 0, abserr = 0))
    if(all(lambda >= -tol_zero)) return(list(Qq = 1, abserr = 0))
    CompQuadForm::imhof(0, lambda = lambda, h = rep.int(1, n),
                        delta = delta2, ...)
}

##### pqfr_davies #####
#' Probability distribution of ratio of quadratic forms
#'
#' \code{pqfr_davies()}: internal function used within \code{pqfr()},
#' based on Davies (1973, 1980).
#'
#' @rdname pqfr
#' @order 6
#'
pqfr_davies <- function(quantile, A, B, mu = rep.int(0, n),
                        tol_zero = .Machine$double.eps * 100, ...) {
    if(!requireNamespace("CompQuadForm", quietly = TRUE)) {
        stop("Package 'CompQuadForm' is required to use this function")
    }
    ## If A or B is missing, let it be an identity matrix
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        A <- diag(n)
    } else {
        n <- dim(A)[1L]
        if(missing(B)) B <- diag(n)
    }
    ## Check basic requirements for arguments
    stopifnot(
        "In pqfr_davies, quantile must be length-one" = (length(quantile) == 1)
    )
    eigA_qB <- eigen(A - quantile * B, symmetric = TRUE)
    lambda <- eigA_qB$values
    delta2 <- c(crossprod(eigA_qB$vectors, mu)) ^ 2
    if(all(lambda <=  tol_zero)) return(list(trace = rep.int(0, 7), ifault = 0, Qq = 0))
    if(all(lambda >= -tol_zero)) return(list(trace = rep.int(0, 7), ifault = 0, Qq = 1))
    CompQuadForm::davies(0, lambda = lambda, h = rep.int(1, n),
                         delta = delta2, sigma = 0, ...)
}


##### dqfr #####
#' Probability distribution of ratio of quadratic forms
#'
#' \code{dqfr()}: \sQuote{Exact} density of the ratio of quadratic forms,
#' \eqn{\frac{ \mathbf{x^\mathit{T} A x} }{ \mathbf{x^\mathit{T} x} }
#'      }{ (x^T A x) / (x^T x) }, where
#' \eqn{\mathbf{x} \sim N_n(\mathbf{0}_n, \mathbf{I}_n)}{x ~ N_n(0_n, I_n)}.
#'
#' @rdname pqfr
#' @order 1
#'
#' @export
#'
dqfr <- function(quantile, A, m = 100L, log = FALSE, ...) {
    n <- dim(A)[1L]
    A <- (A + t(A)) / 2
    ## Check basic requirements for arguments
    stopifnot(
        "A must be a square matrix" = all(c(dim(A)) == n),
        "quantile must be numeric" = is.numeric(quantile),
        "dqfr() does not accommodate B, mu, or Sigma" =
            !("B" %in% ...names()) && !("mu" %in% ...names()) &&
            !("Sigma" %in% ...names())
    )
    LA <- eigen(A, symmetric = TRUE)$values
    ans <- sapply(quantile,
                  function(q) dqfr_A1I1(q, LA, m, ...)$d)
    if(log) ans <- log(ans)
    attributes(ans) <- attributes(quantile)
    return(ans)
}

##### dqfr_A1I1 #####
#' Probability distribution of ratio of quadratic forms
#'
#' \code{dqfr_A1I1()}: internal function used within \code{dqfr()},
#' based on Hillier (2001).
#'
#' @rdname pqfr
#' @order 3
#'
dqfr_A1I1 <- function(quantile, LA, m = 100L,
                      check_convergence = c("relative", "strict_relative",
                                            "absolute", "none"),
                      use_cpp = TRUE,
                      tol_conv = .Machine$double.eps ^ (1/4),
                      thr_margin = 100) {
    if(isTRUE(check_convergence)) check_convergence <- "strict_relative"
    if(isFALSE(check_convergence)) check_convergence <- "none"
    check_convergence <- match.arg(check_convergence)
    n <- length(LA)
    ## Check basic requirements for arguments
    stopifnot(
        "In dqfr_A1I1, quantile must be length-one" = (length(quantile) == 1)
    )
    L1 <- max(LA)
    Ls <- min(LA)
    if(quantile >= L1 || quantile <= Ls) {
        return(list(d = 0, terms = rep.int(0, m + 1)))
    }
    f <- (quantile - Ls) / (L1 - quantile)
    n1 <- sum(LA == L1)
    ns <- sum(LA == Ls)
    if(n1 + ns == n) {
        ## When LA has only two distinct elements, the series vanishes and
        ## the distribution of f reduces to a beta prime distribution
        ans <- f ^ (n1 / 2 - 1) / (1 + f) ^ (n / 2) / beta(n1 / 2, (n - n1) / 2)
        ans <- ans * (L1 - Ls) / (L1 - quantile) ^ 2
        return(list(d = ans, terms = c(ans, rep.int(0, m))))
    }
    ind_psi <- which(LA != L1 & LA != Ls)
    psi <- (LA[ind_psi] - Ls) / (L1 - LA[ind_psi])
    if(use_cpp) {
        cppres <- d_A1I1_Ed(quantile, f, L1, Ls, psi, n, n1, ns, m, thr_margin)
        ansseq <- cppres$ansseq
    } else {
        ind_r1 <- psi > f
        pr <- sum(ind_r1) + n1
        if((pr == n1) || (pr == n - ns)) {
            if(pr == n1) {
                D <- psi / f
                nt <- n1
            } else {
                D <- f / psi
                nt <- ns
            }
            dks <- d1_i(D, m, thr_margin)
            lscf <- attr(dks, "logscale")
            ansseq <- hgs_1d(dks, (2 - nt) / 2, (n - nt) / 2, -lscf)
        } else {
            D1 <- f / psi[ind_r1]
            D2 <- psi[!ind_r1] / f
            dk1 <- d1_i(D1, m, thr_margin)
            dk2 <- d1_i(D2, m, thr_margin)
            lscf1 <- attr(dk1, "logscale")
            lscf2 <- attr(dk2, "logscale")
            alpha <- pr / 2 - 1
            beta <- (n - pr) / 2 - 1
            ansmat <- outer(log(dk1), log(dk2), "+")
            ## c_r(j,k) in Hillier (2001, lemma 2) is
            ## (-1)^(j - k) * gamma(alpha + 1) * gamma(beta + 1) /
            ##   (gamma(alpha + 1 + j - k) * gamma(beta + 1 - j + k)), unless
            ## any of the arguments in the denominator are negative integer or zero
            ordmat <- outer(seq.int(0, m), seq.int(0, m), "-") # j - k
            ansmat <- ansmat + lgamma(alpha + 1) + lgamma(beta + 1) -
                      lgamma(alpha + 1 + ordmat) - lgamma(beta + 1 - ordmat)
            ansmat <- ansmat - lscf1
            ansmat <- t(t(ansmat) - lscf2)
            ansmat <- exp(ansmat)
            sgnmat <- suppressWarnings(sign(gamma(alpha + 1 + ordmat) *
                                            gamma(beta + 1 - ordmat)))
            sgnmat[is.nan(sgnmat)] <- 0
            ansmat <- ansmat * sgnmat
            ansseq <- sum_counterdiag(ansmat)
            ansseq <- ansseq * rep_len(c(1, -1), m + 1L)
        }
        ansseq <- ansseq * exp(-lbeta(pr / 2, (n - pr) / 2) +
                               (-sum(log(psi[ind_r1])) + sum(log(1 + psi)) +
                                (pr - 2) * log(f) - n * log(1 + f)) / 2)
        ansseq <- ansseq * (L1 - Ls) / (L1 - quantile) ^ 2
    }
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  Result truncated before first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        attr(ansseq, "truncated") <- TRUE
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
    list(d = sum(ansseq), terms = ansseq)
}
