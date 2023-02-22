##### qfrm #####
#' Moment of ratio of quadratic forms in normal variables
#'
#' \code{qfrm()} is a front-end function to obtain the (compound) moment
#' of a ratio of quadratic forms in normal variables, i.e.,
#' \eqn{ \mathrm{E} \left(
#'   \frac{(\mathbf{x^\mathit{T} A x})^p }{(\mathbf{x^\mathit{T} B x})^q}
#'   \right) },
#' where \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{\Sigma})}.
#' Internally, \code{qfrm()} calls one of the following functions which does
#' the actual calculation, depending on \eqn{\mathbf{A}}, \eqn{\mathbf{B}},
#' and \eqn{p}. Usually the best one is automatically selected.
#'
#' These functions use infinite series expressions based on the joint
#' moment-generating function (with the top-order zonal/invariant polynomials)
#' (see Smith 1989, Hillier et al. 2009, 2014; Bao & Kan 2013), and the results
#' are typically partial (truncated) sums from these infinite series,
#' which necessarily involve truncation errors.
#' (An exception is when \eqn{\mathbf{B} = \mathbf{I}_n} and \eqn{p} is a
#' positive integer, the case handled by \code{qfrm_ApIq_int()}.)
#'
#' The returned value is a list consisting of the truncated sequence
#' up to the order specified by \code{m}, its sum,
#' and error bounds corresponding to these (see "Values").
#' The \code{print} method only displays the terminal partial sum and its
#' error bound (when available).
#' Use \code{plot()} for visual inspection, or the ordinary list
#' element access as required.
#'
#' In most cases, \code{p} and \code{q} should be nonnegative
#' (in addition, \code{p} should be an integer in
#' \code{qfrm_ApIq_int()} and \code{qfrm_ApBq_int()} when used directly),
#' and an error is thrown otherwise. The only exception is
#' \code{qfrm_ApIq_npi()} which accepts negative exponents to accommodate
#' \eqn{\frac{(\mathbf{x^\mathit{T} x})^q }{(\mathbf{x^\mathit{T} A x})^p}}.
#' Even in the latter case, the exponents should have the same sign.
#' (Technically, not all of these conditions are necessary for the mathematical
#' results to hold, but they are enforced for simplicity).
#'
#' When \code{error_bound = TRUE} (default), \code{qfrm_ApBq_int()} evaluates
#' a truncation error bound following Hillier et al. (2009: theorem 6) or
#' Hillier et al. (2014: theorem 7) (for zero and nonzero means, respectively).
#' \code{qfrm_ApIq_npi()} implements similar error bounds.
#' No error bound is known for \code{qfrm_ApBq_npi()} to the
#' author's knowledge.
# #' See \code{vignette("qfratio")} for further technical details.
#'
#' For situations when the error bound is unavailable, a *very rough* check of
#' numerical convergence is also conducted; a warning is thrown if
#' the magnitude of the last term does not look small enough.
#' By default, its relative magnitude to the sum is compared with
#' the tolerance controlled by \code{tol_conv}, whose default is
#' \code{.Machine$double.eps^(1/4)} (= ~\code{1.2e-04})
#' (see \code{check_convergence}).
#'
#' When \code{Sigma} is provided, the quadratic forms are transformed into
#' a canonical form; that is, using the decomposition
#' \eqn{\mathbf{\Sigma} = \mathbf{K} \mathbf{K}^T}, where the number of
#' columns \eqn{m} of \eqn{\mathbf{K}} equals the rank of \eqn{\mathbf{\Sigma}},
#' \eqn{\mathbf{A}_\mathrm{new} = \mathbf{K^\mathit{T} A K}},
#' \eqn{\mathbf{B}_\mathrm{new} = \mathbf{K^\mathit{T} B K}}, and
#' \eqn{\mathbf{x}_\mathrm{new} = \mathbf{K}^{-} \mathbf{x}
#'      \sim N(\mathbf{K}^{-} \bm{\mu}, \mathbf{I}_m)}.
#' \code{qfrm()} handles this by transforming \code{A}, \code{B},
#' and \code{mu} and calling itself recursively with these new arguments.
#' Note that the ``internal'' functions do not accommodate \code{Sigma}
#' (the error for unused arguments will happen).
#' For singular \eqn{\mathbf{\Sigma}}, one of the following conditions should
#' be met for the above transformation to be valid:
#' **1**) \eqn{\bm{\mu}} is in the range of \eqn{\mathbf{\Sigma}};
#' **2**) \eqn{\mathbf{A}} and \eqn{\mathbf{B}} are in the range of
#' \eqn{\mathbf{\Sigma}}; or
#' **3**) \eqn{\mathbf{A} \bm{\mu} = \mathbf{B} \bm{\mu} = \mathbf{0}}.
#' An error is thrown if none is met with a singular \code{Sigma}.
#'
#' The existence of the moment is assessed by the eigenstructures of
#' \eqn{\mathbf{A}} and \eqn{\mathbf{B}}, \eqn{p}, and \eqn{q}, according to
#' Bao & Kan (2013: proposition 1). An error will result if the conditions
#' are not met.
#'
#' Straightforward implementation of the original recursive algorithms can
#' suffer from numerical overflow when the problem is large.
#' Internal functions (\code{\link{d1_i}}, \code{\link{d2_ij}},
#' \code{\link{d3_ijk}}) are desinged to avoid overflow by order-wise scaling.
#' However, when evaluation of multiple series is required
#' (\code{qfrm_ApIq_npi()} with nonzero \code{mu} and \code{qfrm_ApBq_npi()}),
#' the scaling occasionally yields underflow/diminishing of some terms to
#' numerical \code{0}, causing inaccuracy. A warning is
#' thrown in this case. (See also "Scaling" in \code{\link{d1_i}}.)
#' To avoid this problem, the \code{C++} versions of these functions have two
#' workarounds, as controlled by \code{cpp_method}.
#' **1**) The \code{"long_double"} option uses the \code{long double} variable
#' type instead of the regular \code{double}. This is generally slow and
#' most memory-inefficient.
#' **2**) The \code{"coef_wise"} option uses a coefficient-wise scaling
#' algorithm with the \code{double} variable type. This is generally robust
#' against underflow issues. Computational time varies a lot with conditions;
#' generally only modestly slower than the \code{"double"} option, but can be
#' the slowest in some extreme conditions.
#'
#' For the sake of completeness (only), the scaling parameters \eqn{\alpha} and
#' \eqn{\beta} (see, e.g., Bau & Kan 2013: eqs. 10 and 12) can be modified via
#' the arguments \code{alphaA} and \code{alphaB}. These are the factors for
#' the inverses of the largest eigenvalues of \eqn{\mathbf{A}} and
#' \eqn{\mathbf{B}}, respectively, and should be between 0 and 2.
#' The default is 1, which should suffice for most purposes.
#' Values larger than 1 often yield faster convergence, but are *not*
#' recommended as the error bound will not strictly hold
#' (see Hillier et al. 2009, 2014).
#'
#' ## Multithreading:
#' All these functions use \code{C++} versions to speed up computation
#' by default.
#' Furthermore, some of the \code{C++} functions, in particular those
#' using more than one matrix arguments, are parallelized with \code{OpenMP}
#' (when available). Use the argument \code{nthreads} to control the number
#' of \code{OpenMP} threads. By default (\code{nthreads = 0}), one-half of
#' the processors detected with \code{omp_get_num_procs()} are used.
#' This is except when all the argument matrices share the same eigenvectors
#' and hence the calculation only involves element-wise operations of
#' eigenvalues. In that case, the calculation is typically fast without
#' parallelization, so \code{nthreads} is automatically set to \code{1}
#' unless explicitly specified otherwise; the user can still specify
#' a larger value or \code{0} for (typically marginal) speed gains in large
#' problems.
#'
#' @param A,B
#'   Argument matrices. Should be square. Will be automatically symmetrized.
#' @param p,q
#'   Exponents corresponding to \eqn{\mathbf{A}} and \eqn{\mathbf{B}},
#'   respectively. When only one is provided, the other is set to the same
#'   value. Should be length-one numeric (see "Details" for further conditions).
#' @param m
#'   Order of polynomials at which the series expression is truncated.
#'   \eqn{M} in Hillier et al. (2009, 2014).
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}} for \eqn{\mathbf{x}}.
#'   Accommodated only by the front-end \code{qfrm()}. See "Details".
#' @param tol_zero
#'   Tolerance against which numerical zero is determined.  Used to determine,
#'   e.g., whether \code{mu} is a zero vector, \code{A} or \code{B} equals
#'   the identity matrix, etc.
#' @param tol_sing
#'   Tolerance against which matrix singularity and rank are determined.
#'   The eigenvalues smaller than this are considered zero.
#' @param ...
#'   Additional arguments in the front-end \code{qfrm()} will be passed to
#'   the appropriate ``internal'' function.
#' @param alphaA,alphaB
#'   Factors for the scaling constants for \eqn{\mathbf{A}} and
#'   \eqn{\mathbf{B}}, respectively. See "Details".
#' @param use_cpp
#'   Logical to specify whether the calculation is done with \code{C++}
#'   functions via \code{Rcpp}. \code{TRUE} by default.
#' @param cpp_method
#'   Method used in \code{C++} calculations to avoid numerical
#'   overflow/underflow (see "Details"). Options:
#'   \itemize{
#'     \item{\code{"double"}: }{default; fastest but prone to underflow in
#'           some conditions}
#'     \item{\code{"long_double"}: }{same algorithm but using the
#'           \code{long double} variable type; robust but slow and
#'           memory-inefficient}
#'     \item{\code{"coef_wise"}: }{coefficient-wise scaling algorithm;
#'           most robust but variably slow}
#'    }
#' @param error_bound
#'   Logical to specify whether an error bound is returned (if available).
#' @param check_convergence
#'   Specifies how numerical convergence is checked (see "Details"). Options:
#'   \itemize{
#'     \item{\code{"relative"}: }{default; magnitude of the last term of
#'           the series relative to the sum is compared with \code{tol_conv}}
#'     \item{\code{"strict_relative"} or \code{TRUE}: }{same, but stricter than
#'           default by setting \code{tol_conv = .Machine$double.eps}}
#'     \item{\code{"absolute"}: }{absolute magnitude of the last term is
#'           compared with \code{tol_conv}}
#'     \item{\code{"none"} or \code{FALSE}: }{skips convergence check}
#'   }
#' @param tol_conv
#'   Tolerance against which numerical convergence of series is checked.
#'   Used with \code{check_convergence}.
#' @param thr_margin
#'   Optional argument to adjust the threshold for scaling (see "Scaling"
#'   in \code{\link{d1_i}}). Passed to internal functions (\code{\link{d1_i}},
#'   \code{\link{d2_ij}}, \code{\link{d3_ijk}}) or their \code{C++} equivalents.
#' @param nthreads
#'   Number of threads used in OpenMP-enabled \code{C++} functions.
#'   \code{0} or any negative value is special and means one-half of
#'   the number of processors detected. See "Multithreading" in "Details".
#'
#' @return
#' A \code{\link[=new_qfrm]{qfrm}} object consisting of the following:
#' \itemize{
#'   \item{\code{$statistic}: }{evaluation result (\code{sum(terms)})}
#'   \item{\code{$terms}: }{vector of \eqn{0}th to \eqn{m}th order terms}
#'   \item{\code{$error_bound}: }{error bound of \code{statistic}}
#'   \item{\code{$seq_error}: }{vector of error bounds corresponding to
#'                            partial sums (\code{cumsum(terms)})}
#'  }
#'
#' @references
#' Bao, Y. & Kan, R. (2013). On the moments of ratios of quadratic forms in
#'   normal random variables. *Journal of Multivariate Analysis*, **117**,
#'   229--245.
#'   doi:[10.1016/j.jmva.2013.03.002](https://doi.org/10.1016/j.jmva.2013.03.002).
#'
#' Hillier, G., Kan, R, & Wang, X. (2009). Computationally efficient recursions
#'   for top-order invariant polynomials with applications.
#'   *Econometric Theory*, **25**, 211--242.
#'   doi:[10.1017/S0266466608090075](https://doi.org/10.1017/S0266466608090075).
#'
#' Hillier, G., Kan, R, & Wang, X. (2014). Generating functions and
#'   short recursions, with applications to the moments of quadratic forms
#'   in noncentral normal vectors. *Econometric Theory*, **30**, 436--473.
#'   doi:[10.1017/S0266466613000364](https://doi.org/10.1017/S0266466613000364).
#'
#' Smith, M. D. (1989). On the expectation of a ratio of quadratic forms
#'   in normal variables. *Journal of Multivariate Analysis*, **31**, 244--257.
#'   doi:[10.1016/0047-259X(89)90065-1](https://doi.org/10.1016/0047-259X(89)90065-1).
#'
#' Smith, M. D. (1993). Expectations of ratios of quadratic forms in normal
#'   variables: evaluating some top-order invariant polynomials.
#'   *Australian Journal of Statistics*, **35**, 271--282.
#'   doi:[10.1111/j.1467-842X.1993.tb01335.x](https://doi.org/10.1111/j.1467-842X.1993.tb01335.x).
#'
#' @seealso \code{\link{qfmrm}} for multiple ratio
#'
#' @name qfrm
#'
#' @export
#'
#' @examples
#' ## Some symmetric matrices and parameters
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(sqrt(1:nv))
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#'
#' ## Expectation of (x^T A x)^2 / (x^T x)^2 where x ~ N(0, I)
#' ## An exact expression is available
#' (res1 <- qfrm(A, p = 2))
#'
#' # The above internally calls the following:
#' qfrm_ApIq_int(A, p = 2) ## The same
#'
#' # Similar result with different expression
#' # This is a suboptimal option and throws a warning
#' qfrm_ApIq_npi(A, p = 2)
#'
#' ## Expectation of (x^T A x)^1/2 / (x^T x)^1/2 where x ~ N(0, I)
#' ## Note how quickly the series converges in this case
#' (res2 <- qfrm(A, p = 1/2))
#' plot(res2)
#'
#' # The above calls:
#' qfrm_ApIq_npi(A, p = 0.5)
#'
#' # This is not allowed (throws an error):
#' \dontrun{qfrm_ApIq_int(A, p = 0.5)}
#'
#' ## (x^T A x)^2 / (x^T B x)^3 where x ~ N(0, I)
#' (res3 <- qfrm(A, B, 2, 3))
#' plot(res3)
#'
#' ## (x^T A x)^2 / (x^T B x)^2 where x ~ N(mu, I)
#' ## Note the two-sided error bound
#' (res4 <- qfrm(A, B, 2, 2, mu = mu))
#' plot(res4)
#'
#' ## (x^T A x)^2 / (x^T B x)^2 where x ~ N(mu, Sigma)
#' (res5 <- qfrm(A, B, p = 2, q = 2, mu = mu, Sigma = Sigma))
#' plot(res5)
#'
#' # Sigma is not allowed in the "internal" functions:
#' \dontrun{qfrm_ApBq_int(A, B, p = 2, q = 2, Sigma = Sigma)}
#'
#' # In res5 above, the error bound didn't converge
#' # Use larger m to evaluate higher-order terms
#' plot(print(qfrm(A, B, p = 2, q = 2, mu = mu, Sigma = Sigma, m = 300)))
#'
qfrm <- function(A, B, p = 1, q = p, m = 100L,
                 mu = rep.int(0, n), Sigma = diag(n),
                 tol_zero = .Machine$double.eps * 100,
                 tol_sing = .Machine$double.eps * 100, ...) {
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
    if(missing(p) && !missing(q)) p <- q
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
                stop("For singular Sigma, certain condition need to be met ",
                     "for A, B, or mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details.")
            }
        }
        return(qfrm(KtAK, KtBK, p, q, m = m, mu = iKmu,
                    tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    if(iseq(B, In, tol_zero)) {
        if((p %% 1) == 0 && p > 0) {
            return(qfrm_ApIq_int(A = A, p = p, q = q, m = m, mu = mu,
                                 tol_zero = tol_zero, ...))
        } else {
            return(qfrm_ApIq_npi(A = A, p = p, q = q, m = m, mu = mu,
                                 tol_zero = tol_zero, tol_sing = tol_sing, ...))
        }
    } else {
        if(iseq(A, In, tol_zero)) {
            return(qfrm_ApIq_npi(A = B, p = -q, q = -p, m = m, mu = mu,
                                 tol_zero = tol_zero, tol_sing = tol_sing, ...))
        }
    }
    if((p %% 1) == 0) {
        return(qfrm_ApBq_int(A = A, B = B, p = p, q = q, m = m, mu = mu,
                             tol_zero = tol_zero, tol_sing = tol_sing, ...))
    } else {
        return(qfrm_ApBq_npi(A = A, B = B, p = p, q = q, m = m, mu = mu,
                             tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
}
##### qfmrm #####
#' Moment of multiple ratio of quadratic forms in normal variables
#'
#' \code{qfmrm()} is a front-end function to obtain the (compound) moment
#' of a multiple ratio of quadratic forms in normal variables in the following
#' special form:
#' \eqn{ \mathrm{E} \left(
#'   \frac{(\mathbf{x^\mathit{T} A x})^p }
#'        {(\mathbf{x^\mathit{T} B x})^q (\mathbf{x^\mathit{T} D x})^r}
#'   \right) },
#' where \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{\Sigma})}.
#' Like \code{qfrm()}, this function calls one of the following ``internal''
#' functions for actual calculation, as appropriate.
#'
#' The usage of these functions is similar to \code{\link{qfrm}}, to which
#' the user is referred for documentation.
#' It is assumed that \eqn{\mathbf{B} \neq \mathbf{D}}
#' (otherwise, the problem reduces to a simple ratio).
#'
#' When \code{B} is identity or missing, this and its exponent \code{q} will
#' be swapped with \code{D} and \code{r}, respectively, before
#' \code{qfmrm_ApBIqr_***()} is called.
#'
#' The error bound is only available for \code{qfmrm_ApBIqr_int()}.
#' This is similar to that in \code{qfrm_ApBq_int()} (Watanabe, 2022).
# #' See \code{vignette("qfratio")} for technical details.
#'
#' The existence conditions for the moments of this multiple ratio can be
#' reduced to those for a simple ratio, provided that one of the null spaces
#' of \eqn{\mathbf{B}} and \eqn{\mathbf{D}} is a subspace of the other
#' (including the case they are null).
#' The conditions of Bao & Kan (2013: proposition 1) can then be
#' applied by replacing \eqn{q} and \eqn{m} there by \eqn{q + r} and
#' \eqn{\min{( \mathrm{rank}(\mathbf{B}), \mathrm{rank}(\mathbf{D}) )}},
#' respectively (see also Smith 1989: p. 258 for
#' nonsingular \eqn{\mathbf{B}}, \eqn{\mathbf{D}}).
#' An error is thrown if these conditions are not met in this case.
#' Otherwise (i.e., \eqn{\mathbf{B}} and \eqn{\mathbf{D}} are both singular
#' and neither of their null spaces is a subspace of the other), it seems
#' difficult to define general moment existence conditions.  A sufficient
#' condition can be obtained by applying the same proposition with a new
#' denominator matrix whose null space is union of those of \eqn{\mathbf{B}}
#' and \eqn{\mathbf{D}} (Watanabe, 2022).  A warning is thrown if that
#' condition is not met in this case.
#'
#' Note that these functions may take a substantially longer computational time
#' than those for a simple ratio, because thee former involves
#' multiple infinite series along which summation is to be taken.
#' Expect the computational time to scale with \code{m^2} for
#' \code{qfmrm_IpBDqr_gen()} (when \code{mu} is zero),
#' \code{qfmrm_ApBIqr_int()}, and \code{qfmrm_ApBDqr_int()}, and \code{m^3} for
#' the rest.
#'
#' Most of these functions, excepting \code{qfmrm_ApBIqr_int()} with zero
#' \code{mu}, involve evaluation of multiple series, which can suffer
#' from numerical overflow and underflow (see "Scaling" in
#' \code{\link{d1_i}} and "Details" in \code{\link{qfrm}}). To avoid this,
#' \code{cpp_method = "long_double"} or \code{"coef_wise"} options can be used
#' (see "Details" in \code{\link{qfrm}}).
#'
#' @inheritParams qfrm
#'
#' @param A,B,D
#'   Argument matrices. Should be square. Will be automatically symmetrized.
#' @param p,q,r
#'   Exponents for \eqn{\mathbf{A}}, \eqn{\mathbf{B}}, and \eqn{\mathbf{D}},
#'   respectively. By default, \code{q} equals \code{p/2} and
#'   \code{r} equals \code{q}. If unsure, specify all explicitly.
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}} for \eqn{\mathbf{x}}.
#'   Accommodated only by the front-end \code{qfmrm()}. See "Details"
#'   in \code{\link{qfrm}}.
#' @param alphaA,alphaB,alphaD
#'   Factors for the scaling constants for \eqn{\mathbf{A}},
#'   \eqn{\mathbf{B}}, and \eqn{\mathbf{D}}, respectively. See "Details" in
#'   \code{\link{qfrm}}.
#' @param nthreads
#'   Number of threads used in OpenMP-enabled \code{C++} functions.
#'   See "Multithreading" in \code{\link{qfrm}}.
#' @param ...
#'   Additional arguments in the front-end \code{qfmrm()} will be passed to
#'   the appropriate ``internal'' function.
#'
#' @references
#' Bao, Y. & Kan, R. (2013). On the moments of ratios of quadratic forms in
#'   normal random variables. *Journal of Multivariate Analysis*, **117**,
#'   229--245.
#'   doi:[10.1016/j.jmva.2013.03.002](https://doi.org/10.1016/j.jmva.2013.03.002).
#'
#' Smith, M. D. (1989). On the expectation of a ratio of quadratic forms
#'   in normal variables. *Journal of Multivariate Analysis*, **31**, 244--257.
#'   doi:[10.1016/0047-259X(89)90065-1](https://doi.org/10.1016/0047-259X(89)90065-1).
#'
#' Watanabe, J. (2022). Exact expressions and numerical evaluation of average
#'   evolvability measures for characterizing and comparing **G** matrices.
#'   *bioRxiv* preprint, 2022.11.02.514929.
#'   doi:[10.1101/2022.11.02.514929](https://doi.org/10.1101/2022.11.02.514929).
#'
#' @seealso \code{\link{qfrm}} for simple ratio
#'
#' @name qfmrm
#'
#' @export
#'
#' @examples
#' ## Some symmetric matrices and parameters
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(sqrt(1:nv))
#' D <- diag((1:nv)^2 / nv)
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#'
#' ## Expectation of (x^T A x)^2 / (x^T B x) (x^T x) where x ~ N(0, I)
#' (res1 <- qfmrm(A, B, p = 2, q = 1, r = 1))
#' plot(res1)
#'
#' # The above internally calls the following:
#' qfmrm_ApBIqr_int(A, B, p = 2, q = 1, r = 1) ## The same
#'
#' # Similar result with different expression
#' # This is a suboptimal option and throws a warning
#' qfmrm_ApBIqr_npi(A, B, p = 2, q = 1, r = 1)
#'
#' ## Expectation of (x^T A x) / (x^T B x)^(1/2) (x^T D x)^(1/2) where x ~ N(0, I)
#' (res2 <- qfmrm(A, B, D, p = 1, q = 1/2, r = 1/2))
#' plot(res2)
#'
#' # The above internally calls the following:
#' qfmrm_ApBDqr_int(A, B, D, p = 1, q = 1/2, r = 1/2) ## The same
#'
#' ## Average response correlation between A and B
#' (res3 <- qfmrm(crossprod(A, B), crossprod(A), crossprod(B),
#'                p = 1, q = 1/2, r = 1/2))
#' plot(res3)
#'
#' ## Same, but with x ~ N(mu, Sigma)
#' (res4 <- qfmrm(crossprod(A, B), crossprod(A), crossprod(B),
#'                p = 1, q = 1/2, r = 1/2, mu = mu, Sigma = Sigma))
#' plot(res4)
#'
#' ## Average autonomy of D
#' (res5 <- qfmrm(B = D, D = solve(D), p = 2, q = 1, r = 1))
#' plot(res5)
#'
qfmrm <- function(A, B, D, p = 1, q = p / 2, r = q, m = 100L,
                  mu = rep.int(0, n), Sigma = diag(n),
                  tol_zero = .Machine$double.eps * 100,
                  tol_sing = .Machine$double.eps * 100, ...) {
    ## If A, B, or D is missing, let it be an identity matrix
    ## If they are given, symmetrize
    if(missing(A)) {
        if(missing(B)) {
            if(missing(D)) {
                stop("Provide at least one of A, B, and D")
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
    if(missing(p) && (!missing(q) || !missing(r))) p <- q + r
    zeros <- rep.int(0, n)
    ## If any pair of the three arguments are equal,
    ## reduce the problem to a simple ratio
    if(iseq(B, D, tol_zero)) {
        return(qfrm(A, B, p, q + r, m = m, mu = mu, Sigma = Sigma,
                    tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    if(iseq(A, D, tol_zero) && p >= r) {
        return(qfrm(A, B, p - r, q, m = m, mu = mu, Sigma = Sigma,
                    tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    if(iseq(A, B, tol_zero) && p >= q) {
        return(qfrm(A, D, p - q, r, m = m, mu = mu, Sigma = Sigma,
                    tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    ## If Sigma is given, transform A, B, D, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma) && !iseq(Sigma, In, tol_zero)) {
        KiKS <- KiK(Sigma, tol_sing)
        K <- KiKS$K
        iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        KtBK <- t(K) %*% B %*% K
        KtDK <- t(K) %*% D %*% K
        iKmu <- iK %*% mu
        ## If Sigma is singular, check conditions for A, B, D, mu, and Sigma
        if(ncol(K) != n) {
            okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
                    (iseq(A %*% mu, zeros, tol_zero) &&
                     iseq(B %*% mu, zeros, tol_zero) &&
                     iseq(D %*% mu, zeros, tol_zero)) ||
                    (iseq(crossprod(iK, KtAK %*% iK), A) &&
                     iseq(crossprod(iK, KtBK %*% iK), B) &&
                     iseq(crossprod(iK, KtDK %*% iK), D))
            if(!okay) {
                stop("For singular Sigma, certain condition need to be met ",
                     "for A, B, D, or mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details.")
            }
        }
        return(qfmrm(KtAK, KtBK, KtDK, p, q, r, m = m, mu = iKmu,
                     tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    if(iseq(A, In, tol_zero)) {
        return(qfmrm_IpBDqr_gen(B = B, D = D, p = p, q = q, r = r, m = m,
                                mu = mu,
                                tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
    ## If B == In, swap B and D
    if(iseq(B, In, tol_zero)) {
        B <- D
        D <- In
        qtemp <- q
        q <- r
        r <- qtemp
    }
    if(iseq(D, In, tol_zero)) {
        if((p %% 1) == 0 && p > 0) {
            return(qfmrm_ApBIqr_int(A = A, B = B, p = p, q = q, r = r, m = m,
                                    mu = mu,
                                    tol_zero = tol_zero, tol_sing = tol_sing,
                                    ...))
        } else {
            return(qfmrm_ApBIqr_npi(A = A, B = B, p = p, q = q, r = r, m = m,
                                    mu = mu,
                                    tol_zero = tol_zero, tol_sing = tol_sing,
                                    ...))
        }
    }
    if((p %% 1) == 0) {
        return(qfmrm_ApBDqr_int(A = A, B = B, D = D, p = p, q = q, r = r, m = m,
                                mu = mu,
                                tol_zero = tol_zero, tol_sing = tol_sing,
                                ...))
    } else {
        return(qfmrm_ApBDqr_npi(A = A, B = B, D = D, p = p, q = q, r = r, m = m,
                                mu = mu,
                                tol_zero = tol_zero, tol_sing = tol_sing, ...))
    }
}


###############################
## Function for positive integer moment of a quadratic form
###############################
##### qfpm (documentation) #####
#' Moment of (product of) quadratic forms in normal variables
#'
#' Functions to obtain (compound) moments
#' of a product of quadratic forms in normal variables, i.e.,
#' \eqn{ \mathrm{E} \left(
#'   (\mathbf{x^\mathit{T} A x})^p (\mathbf{x^\mathit{T} B x})^q
#'   (\mathbf{x^\mathit{T} D x})^r \right) },
#' where \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{\Sigma})}.
#'
#' These functions implement the super-short recursion algorithms described in
#' Hillier et al. (2014: sec. 3.1--3.2 and 4).
#' At present, only positive integers are accepted as exponents
#' (negative exponents yield ratios, of course). All these yield exact results.
#'
#' @inheritParams qfmrm
#'
#' @param p,q,r
#'   Exponents for \eqn{\mathbf{A}}, \eqn{\mathbf{B}}, and \eqn{\mathbf{D}},
#'   respectively. By default, these are set to the same value.
#'   If unsure, specify all explicitly.
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}} for \eqn{\mathbf{x}}
#' @param cpp_method
#'   Variable type used in \code{C++} calculations.
#'   In these functions this is ignored.
#'
#' @return
#' A \code{\link[=new_qfrm]{qfpm}} object which has the same elements as those
#' returned by the \code{\link{qfrm}} functions.
#' Use \code{$statistic} to access the value of the moment.
#'
#' @seealso
#' \code{\link{qfrm}} and \code{\link{qfmrm}} for moments of ratios
#'
#' @name qfpm
#'
#' @examples
#' ## Some symmetric matrices and parameters
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(sqrt(1:nv))
#' D <- diag((1:nv)^2 / nv)
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#'
#' ## Expectation of (x^T A x)^2 where x ~ N(0, I)
#' qfm_Ap_int(A, 2)
#'
#' ## This is the same but obviously less efficient
#' qfpm_ABpq_int(A, p = 2, q = 0)
#'
#' ## Expectation of (x^T A x) (x^T B x) (x^T D x) where x ~ N(0, I)
#' qfpm_ABDpqr_int(A, B, D, 1, 1, 1)
#'
#' ## Expectation of (x^T A x) (x^T B x) (x^T D x) where x ~ N(mu, Sigma)
#' qfpm_ABDpqr_int(A, B, D, 1, 1, 1, mu = mu, Sigma = Sigma)
#'
#' ## Expectations of (x^T x)^2 where x ~ N(0, I) and x ~ N(mu, I)
#' ## i.e., roundabout way to obtain moments of
#' ## central and noncentral chi-square variables
#' qfm_Ap_int(diag(nv), 2)
#' qfm_Ap_int(diag(nv), 2, mu = mu)
#'
NULL

##### qfm_Ap_int #####
#' Moment of quadratic forms in normal variables
#'
#' \code{qfm_Ap_int()} is for \eqn{q = r = 0} (simple moment)
#'
#' @rdname qfpm
#'
#' @export
#'
qfm_Ap_int <- function(A, p = 1, mu = rep.int(0, n), Sigma = diag(n),
                       use_cpp = TRUE, cpp_method = "double",
                       tol_zero = .Machine$double.eps * 100,
                       tol_sing = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    n <- ncol(A)
    stopifnot(
        "A should be a square matrix" = all(c(dim(A)) == n),
        "p should be a nonnegative integer" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 0
        },
        "mu should be an n-vector" = length(mu) == n
    )
    zeros <- rep.int(0, n)
    ## If Sigma is given, transform A, B, D, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma) && !iseq(Sigma, diag(n), tol_zero)) {
        KiKS <- KiK(Sigma, tol_sing)
        K <- KiKS$K
        iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        iKmu <- iK %*% mu
        ## If Sigma is singular, check conditions for A, B, D, mu, and Sigma
        if(ncol(K) != n) {
            okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
                    (iseq(A %*% mu, zeros, tol_zero)) ||
                    (iseq(crossprod(iK, KtAK %*% iK), A))
            if(!okay) {
                stop("For singular Sigma, certain condition need to be met ",
                     "for A or mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details.")
            }
        }
        return(qfm_Ap_int(KtAK, p, mu = iKmu, use_cpp = use_cpp,
                           cpp_method = cpp_method, tol_zero = tol_zero))
    }
    A <- (A + t(A)) / 2
    central <- iseq(mu, rep.int(0, n), tol = tol_zero)
    if(use_cpp) {
        if(central) {
            cppres <- Ap_int_cmE(A, p)
        } else {
            cppres <- Ap_int_nmE(A, mu, p)
        }
        ans <- cppres$ans
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        LA <- eigA$values
        if(central) {
            dp <- d1_i(LA, m = p)[p + 1]
        } else {
            mu <- c(crossprod(eigA$vectors, c(mu)))
            dp <- dtil1_i_v(LA, mu, m = p)[p + 1]
        }
        ans <- 2 ^ p * factorial(p) * dp
    }
    new_qfpm(ans)
}

###############################
## Function for product of quadratic forms
###############################
##### qfpm_ABpq_int #####
#' Moment of product of quadratic forms in normal variables
#'
#' \code{qfpm_ABpq_int()} is for \eqn{r = 0}
#'
#' @rdname qfpm
#'
#' @export
#'
qfpm_ABpq_int <- function(A, B, p = 1, q = 1,
                          mu = rep.int(0, n), Sigma = diag(n),
                          use_cpp = TRUE, cpp_method = "double",
                          tol_zero = .Machine$double.eps * 100,
                          tol_sing = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    ## If A or B is missing, let it be an identity matrix
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
    ## Check basic requirements for arguments
    stopifnot(
        "A and B should be square matrices" = all(c(dim(A), dim(B)) == n),
        "p and q should be nonnegative integers" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 0 &&
            length(q) == 1 &&
            (q %% 1) == 0 &&
            q >= 0
        },
        "mu should be an n-vector" = length(mu) == n
    )
    ## When p = 0 and use_cpp = FALSE, out of bound error happens in d2_pj_*
    if(p == 0) {
        return(qfm_Ap_int(A = B, p = q, mu = mu, Sigma = Sigma,
                          use_cpp = use_cpp, cpp_method = cpp_method,
                          tol_zero = tol_zero, tol_sing = tol_sing))
    }
    zeros <- rep.int(0, n)
    ## If Sigma is given, transform A, B, D, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma) && !iseq(Sigma, In, tol_zero)) {
        KiKS <- KiK(Sigma, tol_sing)
        K <- KiKS$K
        iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        KtBK <- t(K) %*% B %*% K
        iKmu <- iK %*% mu
        ## If Sigma is singular, check conditions for A, B, D, mu, and Sigma
        if(ncol(K) != n) {
            okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
                    (iseq(A %*% mu, zeros, tol_zero) &&
                     iseq(B %*% mu, zeros, tol_zero)) ||
                    (iseq(crossprod(iK, KtAK %*% iK), A) &&
                     iseq(crossprod(iK, KtBK %*% iK), B))
            if(!okay) {
                stop("For singular Sigma, certain condition need to be met ",
                     "for A, B, or mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details.")
            }
        }
        return(qfpm_ABpq_int(KtAK, KtBK, p, q, mu = iKmu,
                             use_cpp = use_cpp, cpp_method = cpp_method,
                             tol_zero = tol_zero))
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol = tol_zero)
    if(use_vec) LA <- diag(A)
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ABpq_int_cvE(LA, LB, p, q)
            } else {
                cppres <- ABpq_int_cmE(A, LB, p, q)
            }
        } else {
            if(use_vec) {
                cppres <- ABpq_int_nvE(LA, LB, mu, p, q)
            } else {
                cppres <- ABpq_int_nmE(A, LB, mu, p, q)
            }
        }
        ans <- cppres$ans
    } else {
        if(central) {
            if(use_vec) {
                dpq <- d2_pj_v(LA, LB, m = p + q, p)[p + 1, q + 1]
            } else {
                dpq <- d2_pj_m(A, diag(LB), m = p + q, p)[p + 1, q + 1]
            }
        } else {
            if(use_vec) {
                dpq <- dtil2_pq_v(LA, LB, mu, p, q)[p + 1, q + 1]
            } else {
                dpq <- dtil2_pq_m(A, diag(LB), mu, p, q)[p + 1, q + 1]
            }
        }
        ans <- 2 ^ (p + q) * factorial(p) * factorial(q) * dpq
    }
    new_qfpm(ans)
}

##### qfpm_ABDpqr_int #####
#' Moment of product of quadratic forms in normal variables
#'
#' \code{qfpm_ABDpqr_int()} is for the product of all three powers
#'
#' @rdname qfpm
#'
#' @export
#'
qfpm_ABDpqr_int <- function(A, B, D, p = 1, q = 1, r = 1,
                            mu = rep.int(0, n), Sigma = diag(n),
                            use_cpp = TRUE, cpp_method = "double",
                            tol_zero = .Machine$double.eps * 100,
                            tol_sing = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    ## If A, B, or D is missing, let it be an identity matrix
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
    ## Check basic requirements for arguments
    stopifnot(
        "A, B, and D should be square matrices" =
            all(c(dim(A), dim(B), dim(D)) == n),
        "p, q, and r should be nonnegative integers" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 0 &&
            length(q) == 1 &&
            (q %% 1) == 0 &&
            q >= 0 &&
            length(r) == 1 &&
            (r %% 1) == 0 &&
            r >= 0
        },
        "mu should be an n-vector" = length(mu) == n
    )
    ## When (p = 0 or q = r = 0) and use_cpp = FALSE, out of bound error happens
    if(p == 0) {
        return(qfpm_ABpq_int(A = B, B = D, p = q, q = r, mu = mu, Sigma = Sigma,
                             use_cpp = use_cpp, cpp_method = cpp_method,
                             tol_zero = tol_zero, tol_sing = tol_sing))
    }
    if(q == 0 && r == 0) {
        return(qfm_Ap_int(A = A, p = p, mu = mu, Sigma = Sigma,
                          use_cpp = use_cpp, cpp_method = cpp_method,
                          tol_zero = tol_zero, tol_sing = tol_sing))
    }
    zeros <- rep.int(0, n)
    ## If Sigma is given, transform A, B, D, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma) && !iseq(Sigma, In, tol_zero)) {
        KiKS <- KiK(Sigma, tol_sing)
        K <- KiKS$K
        iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        KtBK <- t(K) %*% B %*% K
        KtDK <- t(K) %*% D %*% K
        iKmu <- iK %*% mu
        ## If Sigma is singular, check conditions for A, B, D, mu, and Sigma
        if(ncol(K) != n) {
            okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
                    (iseq(A %*% mu, zeros, tol_zero) &&
                     iseq(B %*% mu, zeros, tol_zero) &&
                     iseq(D %*% mu, zeros, tol_zero)) ||
                    (iseq(crossprod(iK, KtAK %*% iK), A) &&
                     iseq(crossprod(iK, KtBK %*% iK), B) &&
                     iseq(crossprod(iK, KtDK %*% iK), D))
            if(!okay) {
                stop("For singular Sigma, certain condition need to be met ",
                     "for A, B, D, or mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details.")
            }
        }
        return(qfpm_ABDpqr_int(KtAK, KtBK, KtDK, p, q, r, mu = iKmu,
                               use_cpp = use_cpp, cpp_method = cpp_method,
                               tol_zero = tol_zero))
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero, TRUE) && is_diagonal(D, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol = tol_zero)
    if(use_vec) {
        LA <- diag(A)
        LD <- diag(D)
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ABDpqr_int_cvE(LA, LB, LD, p, q, r)
            } else {
                cppres <- ABDpqr_int_cmE(A, LB, D, p, q, r)
            }
        } else {
            if(use_vec) {
                cppres <- ABDpqr_int_nvE(LA, LB, LD, mu, p, q, r)
            } else {
                cppres <- ABDpqr_int_nmE(A, LB, D, mu, p, q, r)
            }
        }
        ans <- cppres$ans
    } else {
        if(central) {
            if(use_vec) {
                dpqr <- d3_pjk_v(LA, LB, LD, q + r, p)[p + 1, q + 1, r + 1]
            } else {
                dpqr <- d3_pjk_m(A, diag(LB), D, q + r, p)[p + 1, q + 1, r + 1]
            }
        } else {
            if(use_vec) {
                dpqr <- dtil3_pqr_v(LA, LB, LD, mu,
                                    p, q, r)[p + 1, q + 1, r + 1]
            } else {
                dpqr <- dtil3_pqr_m(A, diag(LB), D, mu,
                                    p, q, r)[p + 1, q + 1, r + 1]
            }
        }
        ans <- 2 ^ (p + q + r) *
               factorial(p) * factorial(q) * factorial(r) * dpqr
    }
    new_qfpm(ans)
}




###############################
## Functions for simple ratio
###############################

##### qfrm_ApIq_int #####
## The use_cpp part for the noncentral case is suboptimal as it calls
## gsl::hyperg_1F1 from R.  The original C++ function could be called from
## the GSL library, but this requires "LinkingTo: RcppGSL" and separate
## installation of the library itself (hence is not portable).
## To turn this on, do the following:
## - DESCRIPTION: LinkingTo: RcppGSL (and SystemRequirements: GNU GSL?)
## - src/ratio_funs.cpp: Turn on preprocessor directives (around lines 7-9) and
##     relevant lines in ApIq_int_nmE
## - src/Makevars(.win): Turn flags in PKG_CXXFLAGS and PKG_LIBS
##
#' Positive integer moment of ratio of quadratic forms in normal variables
#'
#' \code{qfrm_ApIq_int()}: For \eqn{\mathbf{B} = \mathbf{I}_n} and
#' positive-integral \eqn{p}.
#'
#' ## Dependency note:
#' An exact expression of the moment is available when
#' \eqn{p} is integer and \eqn{\mathbf{B} = \mathbf{I}_n}
#' (handled by \code{qfrm_ApIq_int()}), but this requires evaluation of
#' a confluent hypergeometric function when \eqn{\bm{\mu}} is nonzero
#' (Hillier et al. 2014: theorem 4).
#' This is done via \code{gsl::hyperg_1F1()} if the package \code{gsl} is
#' installed (which this package \code{Suggests}). Otherwise, the function uses
#' the ordinary infinite series expression (Hillier et al. 2009), which is
#' less accurate and slow, and throws a message (once per session).
#' It is recommended to install that package if an accurate estimate
#' is desired for that case.
#'
# #' @importFrom gsl hyperg_1F1
#'
#' @rdname qfrm
#'
#' @export
#'
qfrm_ApIq_int <- function(A, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                          use_cpp = TRUE, cpp_method = "double",
                          tol_zero = .Machine$double.eps * 100,
                          thr_margin = 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    n <- ncol(A)
    stopifnot(
        "A should be a square matrix" = all(c(dim(A)) == n),
        "p should be a positive integer" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 1
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "mu should be an n-vector" = length(mu) == n
    )
    A <- (A + t(A)) / 2
    central <- iseq(mu, rep.int(0, n), tol = tol_zero)
    exact <- TRUE
    if(use_cpp) {
        if(central) {
            cppres <- ApIq_int_cmE(A, p, q)
            ans <- cppres$ans
            ansseq <- ans
        } else {
            cppres <- ApIq_int_nmE(A, mu, p, q)
            ansseq <- cppres$ansseq
            ans <- sum(ansseq)
        }
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        LA <- eigA$values
        if(central) {
            dp <- d1_i(LA, m = p)[p + 1]
            ans <- exp((p - q) * log(2) + lgamma(p + 1) + lgamma(n/2 + p - q)
                       - lgamma(n/2 + p)) * dp
            ansseq <- ans
        } else {
            mu <- c(mu)
            if(requireNamespace("gsl", quietly = TRUE)) {
                ## This is an exact expression (Hillier et al. 2014, (58))
                mu <- c(crossprod(eigA$vectors, c(mu)))
                aps <- a1_pk(LA, mu, m = p)[p + 1, ]
                ls <- 0:p
                ansseq <-
                    exp((p - q) * log(2) + lgamma(p + 1)
                      + lgamma(n/2 + p - q + ls) - ls * log(2)
                      - lgamma(ls + 1) - lgamma(n/2 + p + ls)) *
                    gsl::hyperg_1F1(q, n / 2 + p + ls, -crossprod(mu) / 2) * aps
            } else {
                ## This is a recursive alternative (Hillier et al. 2014, (53))
                ## which is less accurate (by truncation) and slower
                rlang::inform(
                    paste0("  When using qfrm_ApIq_int() with nonzero mu, it ",
                           "is recommended to\n  install the package \"gsl\", ",
                           "with which an exact result is available.\n  See ",
                           "\"Details\" in documentation of this function."),
                    .frequency = "once", .frequency_id = "use_gsl")
                dks <- d2_ij_m(A, tcrossprod(mu), m, p = p,
                               thr_margin = thr_margin)[p + 1, ]
                ansseq <-
                    exp((p - q) * log(2) + lgamma(1 + p) - c(crossprod(mu)) / 2
                      + lgamma(n/2 + p - q + 0:m) - 0:m * log(2)
                      - lgamma(1/2 + 0:m) + lgamma(1/2) - lgamma(n/2 + p + 0:m)
                      + log(dks))
                exact <- FALSE
            }
            ans <- sum(ansseq)
        }
    }
    if(exact) {
        errseq <- ans - cumsum(ansseq)
    } else {
        errseq <- NA_real_
    }
    new_qfrm(statistic = ans, terms = ansseq, seq_error = errseq, exact = exact)
}

##### qfrm_ApIq_npi #####
#' Non-positive-integer moment of ratio of quadratic forms in normal variables
#'
#' \code{qfrm_ApIq_npi()}: For \eqn{\mathbf{B} = \mathbf{I}_n} and
#' non-positive-integral \eqn{p} (fraction or negative).
#'
#' @rdname qfrm
#'
#' @export
#'
qfrm_ApIq_npi <- function(A, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                    error_bound = TRUE,
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 1, alphaA = 1,
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps * 100,
                    thr_margin = 100) {
    if(isTRUE(check_convergence)) check_convergence <- "strict_relative"
    if(isFALSE(check_convergence)) check_convergence <- "none"
    check_convergence <- match.arg(check_convergence)
    if(!missing(cpp_method)) use_cpp <- TRUE
    cpp_method <- match.arg(cpp_method)
    n <- ncol(A)
    stopifnot(
        "A should be a square matrix" = all(c(dim(A)) == n),
        # "mu should be an n-vector" = length(mu) == n,
        "p should be a real number" = length(p) == 1,
        "q should be a real number" = length(q) == 1,
        "p and q should have the same sign" = p * q >= 0
    )
    if((p %% 1) == 0 && p > 0) {
        warning("For integral p, qfrm_ApIq_int() works better")
    }
    A <- (A + t(A)) / 2
    eigA <- eigen(A, symmetric = TRUE)
    LA <- eigA$values
    UA <- eigA$vectors
    if(any(LA < -tol_sing) && ((p %% 1) != 0 || p < 0)) {
        stop("Detected negative eigenvalue(s) of A (< -tol_sing), ",
             "with which\n  non-integer power of quadratic form is not ",
             "well defined.\n  If you know them to be 0, use larger tol_sing ",
             "to suppress this error.")
    }
    ## Check condition for existence of moment (Bao & Kan, 2013, prop. 1)
    cond_exist <- n / 2 + p > q ## condition(1)
    stopifnot("Moment does not exist in this combination of p, q, and rank(B)" =
                  cond_exist)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    # diminished <- FALSE
    bA <- alphaA / max(abs(LA))
    if(use_cpp) {
        if(central) {
            cppres <- ApIq_npi_cvE(LA, bA, p, q, m, error_bound, thr_margin)
        } else {
            if(cpp_method == "coef_wise") {
                cppres <- ApIq_npi_nvEc(LA, UA, bA, mu, p, q, m,
                                        thr_margin, nthreads)
            } else if(cpp_method == "long_double") {
                cppres <- ApIq_npi_nvEl(LA, UA, bA, mu, p, q, m,
                                        thr_margin, nthreads)
            } else {
                cppres <- ApIq_npi_nvE(LA, UA, bA, mu, p, q, m,
                                       thr_margin, nthreads)
            }
            # diminished <- cppres$diminished
        }
        ansseq <- cppres$ansseq
    } else {
        LAh <- rep.int(1, n) - bA * LA
        if(central) {
            dks <- d1_i(LAh, m = m, thr_margin = thr_margin)
            lscf <- attr(dks, "logscale")
            attributes(dks) <- NULL
            ansseq <- hgs_1d(dks, -p, n / 2, ((p - q) * log(2) - p * log(bA)
                             + lgamma(n/2 + p - q) - lgamma(n/2) - lscf))
        } else {
            mu <- c(crossprod(UA, c(mu)))
            # ## This is based on recursion for d as in Hillier et al. (2009)
            # dks <- d2_ij_m(diag(LAh), tcrossprod(mu), m)
            # ansmat <- hgs_dmu_2d(dks, -p, n / 2 + p - q, n / 2,
            #                      ((p - q) * log(2) - c(crossprod(mu)) / 2
            #                      - p * log(bA) + lgamma(n/2 + p - q) - lgamma(n/2)))
            ## This is based on recursion for h as in Hillier et al. (2014)
            dks <- h2_ij_v(LAh, rep.int(0, n), mu, m, thr_margin = thr_margin)
            lscf <- attr(dks, "logscale")
            ansmat <- hgs_2d(dks, -p, q, n / 2, ((p - q) * log(2) - p * log(bA)
                             + lgamma(n/2 + p - q) - lgamma(n/2) - lscf))
            ansseq <- sum_counterdiag(ansmat)
            # diminished <- any(lscf < 0) && any(diag(dks[(m + 1):1, ]) == 0)
        }
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    # if(diminished) {
    #     warning("Some terms in multiple series numerically diminished to 0 ",
    #             "as they were\n  scaled to avoid numerical overflow. ",
    #             "The result will be inaccurate.",
    #             if(cpp_method == "double")
    #                 paste0("\n  Consider using the option cpp_method = ",
    #                        "\"long_double\" or \"coef_wise\"."))
    # }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    singularA <- any(LA < tol_sing)
    alphaout <- alphaA > 1
    if(error_bound && central) {
        if(use_cpp) {
            errseq <- cppres$errseq
            twosided <- FALSE
        } else {
            Lp <- LAh
            dkst <- dks
            twosided <- FALSE
            lcoefe <- (lgamma(-p + 0:m + 1) - lgamma(-p)
                       - lgamma(n/2 + 0:m + 1) + lgamma(n/2 + p - q)
                       + (p - q) * log(2) - p * log(bA))
            errseq <- exp(lcoefe - sum(log(1 - Lp)) / 2) -
                      exp((lcoefe + log(cumsum(dkst[1:(m + 1)] /
                                        exp(lscf[1:(m + 1)] - lscf[m + 1])))) -
                          lscf[m + 1])
            errseq <- errseq * cumprod(sign(-p + 0:m))
        }
        if(singularA) {
            warning("Argument matrix is numerically close to singular.\n  ",
            "If it is singular, this error bound is invalid.")
        }
        if(alphaout) {
            warning("Error bound is unreliable when alphaA > 1\n  ",
            "It is returned purely for heuristic purpose")
        }
    } else {
        if(error_bound) {
            rlang::inform(paste0("Error bound is unavailable for ",
                                 "qfrm_ApIq_npi() when mu is nonzero"),
                          .frequency = "once",
                          .frequency_id = "errorb_ApIq_npi_noncentral")
            errseq <- NA_real_
        } else {
            errseq <- NULL
        }
        twosided <- NULL
    }
    new_qfrm(terms = ansseq, seq_error = errseq, twosided = twosided,
             alphaout = alphaout, singular_arg = singularA)
}


##### qfrm_ApBq_int #####
#' Positive integer moment of ratio of quadratic forms
#'
#' \code{qfrm_ApBq_int()}: For general \eqn{\mathbf{B}} and
#' positive-integral \eqn{p}.
#'
#'
#' @rdname qfrm
#'
#' @export
#'
qfrm_ApBq_int <- function(A, B, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                    error_bound = TRUE,
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE, cpp_method = "double",
                    alphaB = 1,
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps * 100,
                    thr_margin = 100) {
    if(isTRUE(check_convergence)) check_convergence <- "strict_relative"
    if(isFALSE(check_convergence)) check_convergence <- "none"
    check_convergence <- match.arg(check_convergence)
    if(!missing(cpp_method)) use_cpp <- TRUE
    ## If A or B is missing, let it be an identity matrix
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
    ## Check basic requirements for arguments
    stopifnot(
        "A and B should be square matrices" = all(c(dim(A), dim(B)) == n),
        "p should be a positive integer" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 1
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if(iseq(B, In, tol_zero)) {
        warning("For B = I, qfrm_ApIq_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    ## Check condition for existence of moment (Bao & Kan, 2013, prop. 1)
    rB <- sum(LB > tol_sing)
    if(rB == n) {
        cond_exist <- n / 2 + p > q ## condition(1)
    } else {
        A12z <- all(abs(A[1:rB, (rB + 1):n]) < tol_zero)
        A22z <- all(abs(A[(rB + 1):n, (rB + 1):n]) < tol_zero)
        cond_exist <- if(!A22z) {
                    rB / 2 > q              ## condition(2)(iii)
                } else {
                    if(!A12z) {
                        (rB + p) / 2 > q    ## condition(2)(ii)
                    } else {
                        rB / 2 + p > q      ## condiiton(2)(i)
                    }
                }
    }
    stopifnot(
        "B should be nonnegative definite" = all(LB >= -tol_sing),
        "Moment does not exist in this combination of p, q, and rank(B)" =
            cond_exist)
    use_vec <- is_diagonal(A, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    if(use_vec) {
        LA <- diag(A)
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        UA <- eigA$vectors
        LA <- eigA$values
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ApBq_int_cvE(LA, LB, bB, p, q, m,
                                       error_bound, thr_margin)
            } else {
                cppres <- ApBq_int_cmE(A, LA, UA, LB, bB, p, q, m,
                                       error_bound, thr_margin)
            }
        } else {
            if(use_vec) {
                cppres <- ApBq_int_nvE(LA, LB, bB, mu, p, q, m,
                                       error_bound, thr_margin)
            } else {
                cppres <- ApBq_int_nmE(A, LA, UA, LB, bB, mu, p, q, m,
                                       error_bound, thr_margin)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(use_vec) {
            LBh <- rep.int(1, n) - bB * LB
            if(central) {
                dksm <- d2_pj_v(LA, LBh, m, p = p, thr_margin = thr_margin)
            } else {
                dksm <- htil2_pj_v(LA, LBh, mu, m, p = p,
                                   thr_margin = thr_margin)
            }
        } else {
            Bh <- In - bB * diag(LB, nrow = n)
            if(central) {
                dksm <- d2_pj_m(A, Bh, m, p = p, thr_margin = thr_margin)
            } else {
                dksm <- htil2_pj_m(A, Bh, mu, m, p = p, thr_margin = thr_margin)
            }
        }
        dks <- dksm[p + 1, 1:(m + 1)]
        lscf <- attr(dksm, "logscale")[p + 1, 1:(m + 1)]
        ansseq <- hgs_1d(dks, q, n/2 + p,
                         ((p - q) * log(2) + q * log(bB) + lgamma(p + 1)
                          + lgamma(n/2 + p - q) - lgamma(n/2 + p) - lscf))
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    singularB <- any(LB < tol_sing)
    alphaout <- alphaB > 1
    if(error_bound) {
        if(singularB) {
            warning("Argument matrix B is numerically close to singular.\n  ",
                    "When it is singular, this error bound is invalid.\n  ",
                    "If you know it to be, set \"error_bound = FALSE\" ",
                    "to avoid unexpected behaviors.")
        }
        if(alphaout) {
            warning("Error bound is unreliable ",
                    "when alphaB > 1\n  ",
                    "It is returned purely for heuristic purpose")
        }
        if(use_cpp) {
            errseq <- cppres$errseq
            twosided <- cppres$twosided
        } else {
            LAp <- abs(LA)
            if(central) {
                twosided <- any(LA < 0) && ((p %% 2) == 1)
                deldif2 <- 0
                if(twosided) {
                    if(use_vec) {
                        dkstm <- d2_ij_v(LAp, LBh, m, p = p,
                                         thr_margin = thr_margin)
                    } else {
                        Ap <- S_fromUL(UA, LAp)
                        dkstm <- d2_ij_m(Ap, Bh, m, p = p,
                                         thr_margin = thr_margin)
                    }
                } else {
                    Ap <- A
                    dkstm <- dksm
                }
                if(use_vec) {
                    dp <- d1_i(LAp / LB / bB, p, thr_margin = thr_margin)[p + 1]
                } else {
                    Bisqr <- 1 / sqrt(LB)
                    dp <- d1_i(eigen(t(t(Ap * Bisqr) * Bisqr),
                                     symmetric = TRUE)$values / bB, p,
                                     thr_margin = thr_margin)[p + 1]
                }
            } else {
                twosided <- TRUE
                mub <- sqrt(2 / bB) * mu / sqrt(LB)
                deldif2 <- (sum(mub ^ 2) - sum(mu ^ 2)) / 2
                if(use_vec) {
                    dkstm <- hhat2_pj_v(LAp, LBh, mu, m, p = p,
                                        thr_margin = thr_margin)
                    dp <- dtil1_i_v(LAp / LB / bB, mub, p,
                                    thr_margin = thr_margin)[p + 1]
                } else {
                    Ap <- S_fromUL(UA, LAp)
                    dkstm <- hhat2_pj_m(Ap, Bh, mu, m, p = p,
                                        thr_margin = thr_margin)
                    Bisqr <- 1 / sqrt(LB)
                    dp <- dtil1_i_m(t(t(Ap * Bisqr) * Bisqr) / bB,
                                    mub, p, thr_margin = thr_margin)[p + 1]
                }
            }
            dkst <- dkstm[p + 1, 1:(m + 1)]
            lscft <- attr(dkstm, "logscale")[p + 1, 1:(m + 1)]
            lBdet <- sum(log(LB * bB))
            lcoefe <- (lgamma(q + 0:m + 1) - lgamma(q)
                        - lgamma(n/2 + p + 0:m + 1) + lgamma(n/2 + p - q)
                        + (p - q) * log(2) + q * log(bB) + lgamma(p + 1))
            errseq <- exp(lcoefe + (deldif2 + log(dp) - lBdet / 2)) -
                      exp(lcoefe + log(cumsum(dkst / exp(lscft - lscft[m + 1])))
                          - lscft[m + 1])
        }
    } else {
        errseq <- NULL
        twosided <- NULL
    }
    new_qfrm(terms = ansseq, seq_error = errseq, twosided = twosided,
             alphaout = alphaout, singular_arg = singularB)
}

##### qfrm_ApBq_npi #####
#' Non-positive-integer moment of ratio of quadratic forms
#'
#' \code{qfrm_ApBq_npi()}: For general \eqn{\mathbf{B}} and
#' non-integral \eqn{p}.
#'
#' @rdname qfrm
#'
#' @export
#'
qfrm_ApBq_npi <- function(A, B, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 0, alphaA = 1, alphaB = 1,
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
    ## Check basic requirements for arguments
    stopifnot(
        "A and B should be square matrices" = all(c(dim(A), dim(B)) == n),
        "p and q should be nonnegative real numbers" = {
            length(p) == 1 &&
            length(q) == 1 &&
            p >= 0 &&
            q >= 0
        },
        "alphaA should be a scalar with 0 < alphaA < 2" = {
            is.numeric(alphaA) &&
            length(alphaA) == 1 &&
            alphaA > 0 &&
            alphaA < 2
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if((p %% 1) == 0) {
        warning("For integral p, qfrm_ApBq_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    ## Check condition for existence of moment (Bao & Kan, 2013, prop. 1)
    rB <- sum(LB > tol_sing)
    if(rB == n) {
        cond_exist <- n / 2 + p > q ## condition(1)
    } else {
        A12z <- all(abs(A[1:rB, (rB + 1):n]) < tol_zero)
        A22z <- all(abs(A[(rB + 1):n, (rB + 1):n]) < tol_zero)
        cond_exist <- if(!A22z) {
                    rB / 2 > q              ## condition(2)(iii)
                } else {
                    if(!A12z) {
                        (rB + p) / 2 > q    ## condition(2)(ii)
                    } else {
                        rB / 2 + p > q      ## condiiton(2)(i)
                    }
                }
    }
    stopifnot(
        "B should be nonnegative definite" = all(LB >= -tol_sing),
        "Moment does not exist in this combination of p, q, and rank(B)" =
            cond_exist)
    use_vec <- is_diagonal(A, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    if(use_vec && missing(nthreads)) nthreads <- 1
    diminished <- FALSE
    LA <- if(use_vec) diag(A) else eigen(A, symmetric = TRUE)$values
    bA <- alphaA / max(abs(LA))
    if(any(LA < -tol_sing) && (p %% 1) != 0) {
        stop("Detected negative eigenvalue(s) of A (< -tol_sing), ",
             "with which\n  non-integer power of quadratic form is not ",
             "well defined.\n  If you know them to be 0, use larger tol_sing ",
             "to suppress this error.")
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBq_npi_cvEc(LA, LB, bA, bB, p, q, m,
                                            thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBq_npi_cvEl(LA, LB, bA, bB, p, q, m,
                                            thr_margin, nthreads)
                } else {
                    cppres <- ApBq_npi_cvE(LA, LB, bA, bB, p, q, m,
                                           thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBq_npi_cmEc(A, LB, bA, bB, p, q, m,
                                            thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBq_npi_cmEl(A, LB, bA, bB, p, q, m,
                                            thr_margin, nthreads)
                } else {
                    cppres <- ApBq_npi_cmE(A, LB, bA, bB, p, q, m,
                                           thr_margin, nthreads)
                }
            }
        } else {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBq_npi_nvEc(LA, LB, bA, bB, mu, p, q, m,
                                            thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBq_npi_nvEl(LA, LB, bA, bB, mu, p, q, m,
                                            thr_margin, nthreads)
                } else {
                    cppres <- ApBq_npi_nvE(LA, LB, bA, bB, mu, p, q, m,
                                           thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBq_npi_nmEc(A, LB, bA, bB, mu, p, q, m,
                                            thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBq_npi_nmEl(A, LB, bA, bB, mu, p, q, m,
                                            thr_margin, nthreads)
                } else {
                    cppres <- ApBq_npi_nmE(A, LB, bA, bB, mu, p, q, m,
                                           thr_margin, nthreads)
                }
            }
        }
        ansseq <- cppres$ansseq
        diminished <- cppres$diminished
    } else {
        if(use_vec) {
            LAh <- rep.int(1, n) - bA * LA
            LBh <- rep.int(1, n) - bB * LB
            if(central) {
                dksm <- d2_ij_v(LAh, LBh, m, thr_margin = thr_margin)
            } else {
                dksm <- h2_ij_v(LAh, LBh, mu, m, thr_margin = thr_margin)
            }
        } else {
            Ah <- In - bA * A
            Bh <- In - bB * diag(LB, nrow = n)
            if(central) {
                dksm <- d2_ij_m(Ah, Bh, m, thr_margin = thr_margin)
            } else {
                dksm <- h2_ij_m(Ah, Bh, mu, m, thr_margin = thr_margin)
            }
        }
        lscf <- attr(dksm, "logscale")
        ansmat <- hgs_2d(dksm, -p, q, n / 2,
                         ((p - q) * log(2) - p * log(bA) + q * log(bB)
                          + lgamma(n/2 + p - q) - lgamma(n/2) - lscf))
        ansseq <- sum_counterdiag(ansmat)
        diminished <- any(lscf < 0) && any(diag(dksm[(m + 1):1, ]) == 0)
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        # m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    if(diminished) {
        warning("Some terms in multiple series numerically diminished to 0 ",
                "as they were\n  scaled to avoid numerical overflow. ",
                "The result will be inaccurate.",
                if(cpp_method != "coef_wise")
                    paste0("\n  Consider using the option cpp_method = ",
                          if(cpp_method != "long_double") "\"long_double\" or ",
                          "\"coef_wise\"."))
    }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    new_qfrm(terms = ansseq, seq_error = NA_real_)
}



###############################
## Functions for multiple ratio
###############################

##### qfmrm_ApBIqr_int #####
#' Positive integer moment of multiple ratio when D is identity
#'
#' \code{qfmrm_ApBIqr_int()}: For \eqn{\mathbf{D} = \mathbf{I}_n} and
#' positive-integral \eqn{p}
#'
#' @rdname qfmrm
#'
#' @export
#'
qfmrm_ApBIqr_int <- function(A, B, p = 1, q = 1, r = 1, m = 100L,
                    mu = rep.int(0, n),
                    error_bound = TRUE,
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 0, alphaB = 1,
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
    ## Check basic requirements for arguments
    stopifnot(
        "A and B should be square matrices" = all(c(dim(A), dim(B)) == n),
        "p should be a positive integer" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 1
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "r should be a nonnegative real number" = {
            length(r) == 1 &&
            r >= 0
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    ## Check condition for existence of moment (Bao & Kan, 2013, prop. 1)
    rB <- sum(LB > tol_sing)
    if(rB == n) {
        cond_exist <- n / 2 + p > q + r ## condition(1)
    } else {
        A12z <- all(abs(A[1:rB, (rB + 1):n]) < tol_zero)
        A22z <- all(abs(A[(rB + 1):n, (rB + 1):n]) < tol_zero)
        cond_exist <- if(!A22z) {
                    rB / 2 > q + r              ## condition(2)(iii)
                } else {
                    if(!A12z) {
                        (rB + p) / 2 > q + r    ## condition(2)(ii)
                    } else {
                        rB / 2 + p > q + r      ## condiiton(2)(i)
                    }
                }
    }
    stopifnot(
        "B should be nonnegative definite" = all(LB >= -tol_sing),
        "Moment does not exist in this combination of p, q, r, and rank(B)"
        = cond_exist)
    use_vec <- is_diagonal(A, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    # diminished <- FALSE
    if(use_vec) {
        LA <- diag(A)
        LBh <- rep.int(1, n) - bB * LB
        if(missing(nthreads)) nthreads <- 1
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        UA <- eigA$vectors
        LA <- eigA$values
        Bh <- In - bB * diag(LB, nrow = n)
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ApBIqr_int_cvE(LA, LB, bB, p, q, r, m,
                                         error_bound, thr_margin)
            } else {
                cppres <- ApBIqr_int_cmE(A, LA, UA, LB, bB, p, q, r, m,
                                         error_bound, thr_margin)
            }
        } else {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBIqr_int_nvEc(LA, LB, bB, mu, p, q, r, m,
                                              error_bound, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBIqr_int_nvEl(LA, LB, bB, mu, p, q, r, m,
                                              error_bound, thr_margin, nthreads)
                } else {
                    cppres <- ApBIqr_int_nvE(LA, LB, bB, mu, p, q, r, m,
                                             error_bound, thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBIqr_int_nmEc(A, LA, UA, LB, bB, mu, p, q, r, m,
                                              error_bound, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBIqr_int_nmEl(A, LA, UA, LB, bB, mu, p, q, r, m,
                                              error_bound, thr_margin, nthreads)
                } else {
                    cppres <- ApBIqr_int_nmE(A, LA, UA, LB, bB, mu, p, q, r, m,
                                             error_bound, thr_margin, nthreads)
                }
            }
            # diminished <- cppres$diminished
        }
        ansseq <- cppres$ansseq
    } else {
        if(central) {
            if(use_vec) {
                dksm <- d2_pj_v(LA, LBh, m, p = p, thr_margin = thr_margin)
            } else {
                dksm <- d2_pj_m(A, Bh, m, p = p, thr_margin = thr_margin)
            }
            dks <- dksm[p + 1, 1:(m + 1)]
            lscf <- attr(dksm, "logscale")[p + 1, 1:(m + 1)]
            ansseq <- hgs_1d(dks, q, n / 2 + p,
                             ((p - q - r) * log(2) + q * log(bB) + lgamma(p + 1)
                              + lgamma(n/2 + p - q - r) - lgamma(n/2 + p)
                              - lscf))
        } else {
            if(use_vec) {
                dksm <- htil3_pjk_v(LA, LBh, rep.int(0, n), mu, m, p = p, 
                                    thr_margin = thr_margin)
            } else {
                dksm <- htil3_pjk_m(A, Bh, matrix(0, n, n), mu, m, p = p,
                                    thr_margin = thr_margin)
            }
            dks <- dksm[p + 1, , ]
            lscf <- attr(dksm, "logscale")[p + 1, , ]
            ansmat <- hgs_2d(dks, q, r, n / 2 + p,
                             ((p - q - r) * log(2) + q * log(bB) + lgamma(p + 1)
                              + lgamma(n/2 + p - q - r) - lgamma(n/2 + p)
                              - lscf))
            ansseq <- sum_counterdiag(ansmat)
            # diminished <- any(lscf < 0) && any(diag(dks[(m + 1):1, ]) == 0)
        }
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    # if(diminished) {
    #     warning("Some terms in multiple series numerically diminished to 0 ",
    #             "as they were\n  scaled to avoid numerical overflow. ",
    #             "The result will be inaccurate.",
    #             if(cpp_method == "double")
    #                 paste0("\n  Consider using the option cpp_method = ",
    #                        "\"long_double\" or \"coef_wise\"."))
    # }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    singularB <- any(LB < tol_sing)
    alphaout <- alphaB > 1
    if(error_bound) {
        if(singularB) {
            warning("Argument matrix B is numerically close to singular.\n  ",
                    "When it is singular, this error bound is invalid.\n  ",
                    "If you know it to be, set \"error_bound = FALSE\" ",
                    "to avoid unexpected behaviors.")
        }
        if(alphaout) {
            warning("Error bound is unreliable ",
                    "when alphaB > 1\n  ",
                    "It is returned purely for heuristic purpose")
        }
        if(use_cpp) {
            errseq <- cppres$errseq
            twosided <- cppres$twosided
        } else {
            LAp <- abs(LA)
            if(central) {
                twosided <- any(LA < 0) && (p %% 2) == 1
                deldif2 <- 0
                s <- q
                if(twosided) {
                    if(use_vec) {
                        dkstm <- d2_pj_v(LAp, LBh, m, p = p,
                                         thr_margin = thr_margin)
                    } else {
                        Ap <- S_fromUL(UA, LAp)
                        dkstm <- d2_pj_m(Ap, Bh, m, p = p,
                                         thr_margin = thr_margin)
                    }
                } else {
                    Ap <- A
                    dkstm <- dksm
                }
                dkst <- dkstm[p + 1, 1:(m + 1)]
                if(use_vec) {
                    dp <- d1_i(LAp / LB / bB, p, thr_margin = thr_margin)[p + 1]
                } else {
                    Bisqr <- 1 / sqrt(LB)
                    dp <- d1_i(eigen(t(t(Ap * Bisqr) * Bisqr),
                                     symmetric = TRUE)$values / bB, p,
                                     thr_margin = thr_margin)[p + 1]
                }
                lscft <- attr(dkstm, "logscale")[p + 1, ]
            } else {
                twosided <- TRUE
                mub <- sqrt(3 / bB) * mu / sqrt(LB)
                deldif2 <- (sum(mub ^ 2) - sum(mu ^ 2)) / 2
                s <- max(q, r)
                if(use_vec) {
                    dkstm <- hhat3_pjk_v(LAp, LBh, rep.int(0, n), mu, m, p = p,
                                         thr_margin = thr_margin)
                    dp <- dtil1_i_v(LAp / LB / bB, mub, p,
                                    thr_margin = thr_margin)[p + 1]
                } else {
                    Ap <- S_fromUL(UA, LAp)
                    dkstm <- hhat3_pjk_m(Ap, Bh, matrix(0, n, n), mu, m, p = p,
                                         thr_margin = thr_margin)
                    Bisqr <- 1 / sqrt(LB)
                    dp <- dtil1_i_m(t(t(Ap * Bisqr) * Bisqr) / bB,
                                    mub, p, thr_margin = thr_margin)[p + 1]
                }
                dkst <- sum_counterdiag(dkstm[p + 1, , ])
                lscft <- attr(dkstm, "logscale")[p + 1, , 1]
            }
            lBdet <- sum(log(LB * bB))
            lcoefe <- (lgamma(s + 0:m + 1) - lgamma(s)
                       - lgamma(n/2 + p + 0:m + 1) + lgamma(n/2 + p - q - r)
                       + (p - q - r) * log(2) + q * log(bB) + lgamma(p + 1))
            errseq <- exp(lcoefe + (deldif2 + log(dp) - lBdet / 2)) -
                      exp(lcoefe + log(cumsum(dkst / exp(lscft - lscft[m + 1])))
                          - lscft[m + 1])
        }
    } else {
        errseq <- NULL
        twosided <- NULL
    }
    new_qfrm(terms = ansseq, seq_error = errseq,
             twosided = twosided, alphaout = alphaout, singular_arg = singularB)
}

##### qfmrm_ApBIqr_npi #####
#' Non-positive-integer moment of multiple ratio when D is identity
#'
#' \code{qfmrm_ApBIqr_npi()}: For \eqn{\mathbf{D} = \mathbf{I}_n} and
#' non-integral \eqn{p}
#'
#'
#' @rdname qfmrm
#'
#' @export
#'
qfmrm_ApBIqr_npi <- function(A, B, p = 1, q = 1, r = 1, m = 100L,
                    mu = rep.int(0, n),
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 0, alphaA = 1, alphaB = 1,
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
    ## Check basic requirements for arguments
    stopifnot(
        "A and B should be square matrices" = all(c(dim(A), dim(B)) == n),
        "p should be a nonnegative real number" = {
            length(p) == 1 &&
            p >= 0
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "r should be a nonnegative real number" = {
            length(r) == 1 &&
            r >= 0
        },
        "alphaA should be a scalar with 0 < alphaA < 2" = {
            is.numeric(alphaA) &&
            length(alphaA) == 1 &&
            alphaA > 0 &&
            alphaA < 2
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if((p %% 1) == 0) {
        warning("For integral p, qfmrm_ApBIqr_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    ## Check condition for existence of moment (Bao & Kan, 2013, prop. 1)
    rB <- sum(LB > tol_sing)
    if(rB == n) {
        cond_exist <- n / 2 + p > q + r ## condition(1)
    } else {
        A12z <- all(abs(A[1:rB, (rB + 1):n]) < tol_zero)
        A22z <- all(abs(A[(rB + 1):n, (rB + 1):n]) < tol_zero)
        cond_exist <- if(!A22z) {
                    rB / 2 > q + r              ## condition(2)(iii)
                } else {
                    if(!A12z) {
                        (rB + p) / 2 > q + r    ## condition(2)(ii)
                    } else {
                        rB / 2 + p > q + r      ## condiiton(2)(i)
                    }
                }
    }
    stopifnot(
        "B should be nonnegative definite" = all(LB >= -tol_sing),
        "Moment does not exist in this combination of p, q, r, and rank(B)"
            = cond_exist)
    use_vec <- is_diagonal(A, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    diminished <- FALSE
    if(use_vec) {
        LA <- diag(A)
        bA <- alphaA / max(abs(LA))
        LAh <- rep.int(1, n) - bA * LA
        LBh <- rep.int(1, n) - bB * LB
        if(missing(nthreads)) nthreads <- 1
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        LA <- eigA$values
        bA <- alphaA / max(abs(LA))
        Ah <- In - bA * A
        Bh <- In - bB * diag(LB, nrow = n)
    }
    if(any(LA < -tol_sing) && (p %% 1) != 0) {
        stop("Detected negative eigenvalue(s) of A (< -tol_sing), ",
             "with which\n  non-integer power of quadratic form is not ",
             "well defined.\n  If you know them to be 0, use larger tol_sing ",
             "to suppress this error.")
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBIqr_npi_cvEc(LA, LB, bA, bB, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBIqr_npi_cvEl(LA, LB, bA, bB, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- ApBIqr_npi_cvE(LA, LB, bA, bB, p, q, r, m,
                                             thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBIqr_npi_cmEc(A, LB, bA, bB, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBIqr_npi_cmEl(A, LB, bA, bB, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- ApBIqr_npi_cmE(A, LB, bA, bB, p, q, r, m,
                                             thr_margin, nthreads)
                }
            }
        } else {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBIqr_npi_nvEc(LA, LB, bA, bB, mu, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBIqr_npi_nvEl(LA, LB, bA, bB, mu, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- ApBIqr_npi_nvE(LA, LB, bA, bB, mu, p, q, r, m,
                                             thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBIqr_npi_nmEc(A, LB, bA, bB, mu, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBIqr_npi_nmEl(A, LB, bA, bB, mu, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- ApBIqr_npi_nmE(A, LB, bA, bB, mu, p, q, r, m,
                                             thr_margin, nthreads)
                }
            }
        }
        ansseq <- cppres$ansseq
        diminished <- cppres$diminished
    } else {
        if(central) {
            if(use_vec) {
                dksm <- d2_ij_v(LAh, LBh, m, thr_margin = thr_margin)
            } else {
                dksm <- d2_ij_m(Ah, Bh, m, thr_margin = thr_margin)
            }
            lscf <- attr(dksm, "logscale")
            ansmat <- hgs_2d(dksm, -p, q, n / 2, ((p - q - r) * log(2)
                        - p * log(bA) + q * log(bB) + lgamma(n/2 + p - q - r)
                        - lgamma(n/2) - lscf))
            ansseq <- sum_counterdiag(ansmat)
            diminished <- any(lscf < 0) && any(diag(dksm[(m + 1):1, ]) == 0)
        } else {
            if(use_vec) {
                dksm <- h3_ijk_v(LAh, LBh, rep.int(0, n), mu, m,
                                 thr_margin = thr_margin)
            } else {
                dksm <- h3_ijk_m(Ah, Bh, matrix(0, n, n), mu, m,
                                 thr_margin = thr_margin)
            }
            lscf <- attr(dksm, "logscale")
            ansarr <- hgs_3d(dksm, -p, q, r, n / 2, ((p - q - r) * log(2)
                        - p * log(bA) + q * log(bB) + lgamma(n/2 + p - q - r)
                        - lgamma(n/2) - lscf))
            ansseq <- sum_counterdiag3D(ansarr)
            if(any(lscf < 0 )) {
                for(k in 1:(m + 1)) {
                    diminished <-
                        any(diag(dksm[(m + 2 - k):1, 1:(m + 2 - k), k]) == 0)
                    if(diminished) break
                }
            }
        }
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        # m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    if(diminished) {
        warning("Some terms in multiple series numerically diminished to 0 ",
                "as they were\n  scaled to avoid numerical overflow. ",
                "The result will be inaccurate.",
                if(cpp_method != "coef_wise")
                    paste0("\n  Consider using the option cpp_method = ",
                          if(cpp_method != "long_double") "\"long_double\" or ",
                          "\"coef_wise\"."))
    }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    new_qfrm(terms = ansseq, seq_error = NA_real_)
}

##### qfmrm_IpBDqr_gen #####
#' Moment of multiple ratio when A is identity
#'
#' \code{qfmrm_IpBDqr_gen()}: For \eqn{\mathbf{A} = \mathbf{I}_n}
#'
#' @rdname qfmrm
#'
#' @export
#'
qfmrm_IpBDqr_gen <- function(B, D, p = 1, q = 1, r = 1, mu = rep.int(0, n),
                    m = 100L,
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 0, alphaB = 1, alphaD = 1,
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
    if(missing(B)) {
        if(missing(D)) stop("Provide at least one of B and D")
        n <- dim(D)[1L]
        In <- diag(n)
        B <- In
    } else {
        n <- dim(B)[1L]
        In <- diag(n)
        B <- (B + t(B)) / 2
    }
    if(missing(D)) {
        D <- In
    } else {
        D <- (D + t(D)) / 2
    }
    ## Check basic requirements for arguments
    stopifnot(
        "B and D should be square matrices" = all(c(dim(B), dim(D)) == n),
        "p should be a nonnegative real number" = {
            length(p) == 1 &&
            p >= 0
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "r should be a nonnegative real number" = {
            length(r) == 1 &&
            r >= 0
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "alphaD should be a scalar with 0 < alphaD < 2" = {
            is.numeric(alphaD) &&
            length(alphaD) == 1 &&
            alphaD > 0 &&
            alphaD < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate D and mu with eigenvectors of B
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(D, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    diminished <- FALSE
    if(use_vec) {
        LD <- diag(D)
        bD <- alphaD / max(LD)
        LBh <- rep.int(1, n) - bB * LB
        LDh <- rep.int(1, n) - bD * LD
        if(missing(nthreads)) nthreads <- 1
    } else {
        eigD <- eigen(D, symmetric = TRUE)
        LD <- eigD$values
        bD <- alphaD / max(LD)
        Bh <- In - bB * diag(LB, nrow = n)
        Dh <- In - bD * D
    }
    stopifnot("B should be nonnegative definite" = all(LB >= -tol_sing),
              "D should be nonnegative definite" = all(LD >= -tol_sing))
    ## Check condition for existence of moment
    nzB <- (LB > tol_sing)
    nzD <- (LD > tol_sing)
    if(use_vec) {
        ## Common range of B and D
        nzBD <- nzB * nzD
    } else {
        if(all(nzB) && all(nzD)) {
            nzBD <- rep.int(TRUE, n)
        } else {
            zerocols <- cbind(In[, !nzB], eigD$vectors[, !nzD])
            projmat <- zerocols %*% MASS::ginv(crossprod(zerocols)) %*%
                       t(zerocols)
            eigBD <- eigen(In - projmat, symmetric = TRUE)
            nzBD <- eigBD$values > tol_sing
        }
    }
    rBD <- sum(nzBD)
    ## When A == I, the condition simplifies as A12 == 0 and A22 != 0
    if(rBD == n) {
        cond_exist <- n / 2 + p > q + r ## condition(1)
        necess_cond <- TRUE
    } else {
        cond_exist <- rBD / 2 > q + r              ## condition(2)(iii)
        necess_cond <- (rBD == sum(nzB)) || (rBD == sum(nzD))
    }
    if(!cond_exist) {
        if(necess_cond) {
            stop("Moment does not exist in this combination of p, q, r, and",
                 "\n  eigenstructures of A, B, and D")
        } else {
            warning("Moment may not exist in this combination of p, q, r, and",
                    "\n  eigenstructures of A, B, and D")
        }
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- IpBDqr_gen_cvEc(LB, LD, bB, bD, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- IpBDqr_gen_cvEl(LB, LD, bB, bD, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- IpBDqr_gen_cvE(LB, LD, bB, bD, p, q, r, m,
                                             thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- IpBDqr_gen_cmEc(LB, D, bB, bD, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- IpBDqr_gen_cmEl(LB, D, bB, bD, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- IpBDqr_gen_cmE(LB, D, bB, bD, p, q, r, m,
                                             thr_margin, nthreads)
                }
            }
        } else {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- IpBDqr_gen_nvEc(LB, LD, bB, bD, mu, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- IpBDqr_gen_nvEl(LB, LD, bB, bD, mu, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- IpBDqr_gen_nvE(LB, LD, bB, bD, mu, p, q, r, m,
                                             thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- IpBDqr_gen_nmEc(LB, D, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- IpBDqr_gen_nmEl(LB, D, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- IpBDqr_gen_nmE(LB, D, bB, bD, mu,
                                             p, q, r, m, thr_margin, nthreads)
                }
            }
        }
        ansseq <- cppres$ansseq
        diminished <- cppres$diminished
    } else {
        if(central) {
            if(use_vec) {
                dksm <- d2_ij_v(LBh, LDh, m, thr_margin = thr_margin)
            } else {
                dksm <- d2_ij_m(Bh, Dh, m, thr_margin = thr_margin)
            }
            lscf <- attr(dksm, "logscale")
            ansmat <- hgs_2d(dksm, q, r, n / 2,
                             ((p - q - r) * log(2) + q * log(bB) + r * log(bD)
                              + lgamma(n/2 + p - q - r) - lgamma(n/2)
                              - lscf))
            ansseq <- sum_counterdiag(ansmat)
            diminished <- any(lscf < 0) && any(diag(dksm[(m + 1):1, ]) == 0)
        } else {
            if(use_vec) {
                dksm <- h3_ijk_v(rep.int(0, n), LBh, LDh, mu, m,
                                 thr_margin = thr_margin)
            } else {
                dksm <- h3_ijk_m(matrix(0, n, n), Bh, Dh, mu, m,
                                 thr_margin = thr_margin)
            }
            lscf <- attr(dksm, "logscale")
            ansarr <- hgs_3d(dksm, -p, q, r, n / 2,
                             ((p - q - r) * log(2) + q * log(bB) + r * log(bD)
                              + lgamma(n/2 + p - q - r) - lgamma(n/2)
                              - lscf))
            ansseq <- sum_counterdiag3D(ansarr)
            if(any(lscf < 0)) {
                for(k in 1:(m + 1)) {
                    diminished <-
                        any(diag(dksm[(m + 2 - k):1, 1:(m + 2 - k), k]) == 0)
                    if(diminished) break
                }
            }
        }
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        # m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    if(diminished) {
        warning("Some terms in multiple series numerically diminished to 0 ",
                "as they were\n  scaled to avoid numerical overflow. ",
                "The result will be inaccurate.",
                if(cpp_method != "coef_wise")
                    paste0("\n  Consider using the option cpp_method = ",
                          if(cpp_method != "long_double") "\"long_double\" or ",
                          "\"coef_wise\"."))
    }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    new_qfrm(terms = ansseq, seq_error = NA_real_)
}

##### qfmrm_ApBDqr_int #####
#' Positive integer moment of multiple ratio
#'
#' \code{qfmrm_ApBDqr_int()}: For general \eqn{\mathbf{A}}, \eqn{\mathbf{B}},
#' and \eqn{\mathbf{D}}, and positive-integral \eqn{p}
#'
#' @rdname qfmrm
#'
#' @export
#'
qfmrm_ApBDqr_int <- function(A, B, D, p = 1, q = 1, r = 1, m = 100L,
                    mu = rep.int(0, n),
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 0, alphaB = 1, alphaD = 1,
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
    ## Check basic requirements for arguments
    stopifnot(
        "A, B and D should be square matrices" =
            all(c(dim(A), dim(B), dim(D)) == n),
        "p should be a positive integer" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 1
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "r should be a nonnegative real number" = {
            length(r) == 1 &&
            r >= 0
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "alphaD should be a scalar with 0 < alphaD < 2" = {
            is.numeric(alphaD) &&
            length(alphaD) == 1 &&
            alphaD > 0 &&
            alphaD < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate A, D, and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero, TRUE) && is_diagonal(D, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    if(use_vec) {
        LA <- diag(A)
        LD <- diag(D)
        if(missing(nthreads)) nthreads <- 1
    } else {
        eigD <- eigen(D, symmetric = TRUE)
        LA <- eigen(A, symmetric = TRUE)$values
        LD <- eigD$values
    }
    stopifnot("B should be nonnegative definite" = all(LB >= -tol_sing),
              "D should be nonnegative definite" = all(LD >= -tol_sing))
    ## Check condition for existence of moment
    nzB <- (LB > tol_sing)
    nzD <- (LD > tol_sing)
    if(use_vec) {
        ## common nonzero space of B and D, and A "rotated" with its basis
        nzBD <- nzB * nzD
        Ar <- A
    } else {
        if(all(nzB) && all(nzD)) {
            nzBD <- rep.int(TRUE, n)
            Ar <- A
        } else {
            zerocols <- cbind(In[, !nzB], eigD$vectors[, !nzD])
            projmat <- zerocols %*% MASS::ginv(crossprod(zerocols)) %*%
                       t(zerocols)
            eigBD <- eigen(In - projmat, symmetric = TRUE)
            nzBD <- eigBD$values > tol_sing
            Ar <- with(eigBD, crossprod(crossprod(A, vectors), vectors))
        }
    }
    rBD <- sum(nzBD)
    if(rBD == n) {
        cond_exist <- n / 2 + p > q + r ## condition(1)
        necess_cond <- TRUE
    } else {
        A12z <- all(abs(Ar[nzBD, !nzBD]) < tol_zero)
        A22z <- all(abs(Ar[!nzBD, !nzBD]) < tol_zero)
        cond_exist <- if(!A22z) {
                    rBD / 2 > q + r              ## condition(2)(iii)
                } else {
                    if(!A12z) {
                        (rBD + p) / 2 > q + r    ## condition(2)(ii)
                    } else {
                        rBD / 2 + p > q + r      ## condiiton(2)(i)
                    }
                }
        necess_cond <- (rBD == sum(nzB)) || (rBD == sum(nzD))
    }
    if(!cond_exist) {
        if(necess_cond) {
            stop("Moment does not exist in this combination of p, q, r, and",
                 "\n  eigenstructures of A, B, and D")
        } else {
            warning("Moment may not exist in this combination of p, q, r, and",
                    "\n  eigenstructures of A, B, and D")
        }
    }
    bD <- alphaD / max(LD)
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_int_cvEc(LA, LB, LD, bB, bD, p, q, r, m,
                                              thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_int_cvEl(LA, LB, LD, bB, bD, p, q, r, m,
                                              thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_int_cvE(LA, LB, LD, bB, bD, p, q, r, m,
                                             thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_int_cmEc(A, LB, D, bB, bD,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_int_cmEl(A, LB, D, bB, bD,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_int_cmE(A, LB, D, bB, bD,
                                             p, q, r, m, thr_margin, nthreads)
                }
            }
        } else {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_int_nvEc(LA, LB, LD, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_int_nvEl(LA, LB, LD, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_int_nvE(LA, LB, LD, bB, bD, mu,
                                             p, q, r, m, thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_int_nmEc(A, LB, D, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_int_nmEl(A, LB, D, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_int_nmE(A, LB, D, bB, bD, mu,
                                             p, q, r, m, thr_margin, nthreads)
                }
            }
        }
        ansseq <- cppres$ansseq
        diminished <- cppres$diminished
    } else {
        if(use_vec) {
            LBh <- rep.int(1, n) - bB * LB
            LDh <- rep.int(1, n) - bD * LD
            if(central) {
                dksm <- d3_pjk_v(LA, LBh, LDh, m, p = p,
                                 thr_margin = thr_margin)
            } else {
                dksm <- htil3_pjk_v(LA, LBh, LDh, mu, m, p = p,
                                    thr_margin = thr_margin)
            }
        } else {
            Bh <- In - bB * diag(LB, nrow = n)
            Dh <- In - bD * D
            if(central) {
                dksm <- d3_pjk_m(A, Bh, Dh, m, p = p,
                                 thr_margin = thr_margin)
            } else {
                dksm <- htil3_pjk_m(A, Bh, Dh, mu, m, p = p,
                                    thr_margin = thr_margin)
            }
        }
        dks <- dksm[p + 1, , ]
        lscf <- attr(dksm, "logscale")[p + 1, , ]
        ansmat <- hgs_2d(dks, q, r, n / 2 + p,
                         ((p - q - r) * log(2) + q * log(bB) + r * log(bD)
                          + lgamma(p + 1) + lgamma(n/2 + p - q - r)
                          - lgamma(n/2 + p) - lscf))
        ansseq <- sum_counterdiag(ansmat)
        diminished <- any(lscf < 0) && any(diag(dks[(m + 1):1, ]) == 0)
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        # m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    if(diminished) {
        warning("Some terms in multiple series numerically diminished to 0 ",
                "as they were\n  scaled to avoid numerical overflow. ",
                "The result will be inaccurate.",
                if(cpp_method != "coef_wise")
                    paste0("\n  Consider using the option cpp_method = ",
                          if(cpp_method != "long_double") "\"long_double\" or ",
                          "\"coef_wise\"."))
    }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    new_qfrm(terms = ansseq, seq_error = NA_real_)
}

##### qfmrm_ApBDqr_npi #####
#' Positive integer moment of multiple ratio
#'
#' \code{qfmrm_ApBDqr_npi()}: For general \eqn{\mathbf{A}}, \eqn{\mathbf{B}},
#' and \eqn{\mathbf{D}}, and non-integral \eqn{p}
#'
#' @rdname qfmrm
#'
#' @export
#'
qfmrm_ApBDqr_npi <- function(A, B, D, p = 1, q = 1, r = 1,
                    m = 100L, mu = rep.int(0, n),
                    check_convergence = c("relative", "strict_relative",
                                          "absolute", "none"),
                    use_cpp = TRUE,
                    cpp_method = c("double", "long_double", "coef_wise"),
                    nthreads = 0, alphaA = 1, alphaB = 1, alphaD = 1,
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
    ## Check basic requirements for arguments
    stopifnot(
        "A, B and D should be square matrices" =
            all(c(dim(A), dim(B), dim(D)) == n),
        "p should be a nonnegative real number" = {
            length(p) == 1 &&
            p >= 0
        },
        "q should be a nonnegative real number" = {
            length(q) == 1 &&
            q >= 0
        },
        "r should be a nonnegative real number" = {
            length(r) == 1 &&
            r >= 0
        },
        "alphaA should be a scalar with 0 < alphaA < 2" = {
            is.numeric(alphaA) &&
            length(alphaA) == 1 &&
            alphaA > 0 &&
            alphaA < 2
        },
        "alphaB should be a scalar with 0 < alphaB < 2" = {
            is.numeric(alphaB) &&
            length(alphaB) == 1 &&
            alphaB > 0 &&
            alphaB < 2
        },
        "alphaD should be a scalar with 0 < alphaD < 2" = {
            is.numeric(alphaD) &&
            length(alphaD) == 1 &&
            alphaD > 0 &&
            alphaD < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if((p %% 1) == 0) {
        warning("For integral p, qfmrm_ApBDqr_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    bB <- alphaB / max(LB)
    ## Rotate A, D, and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero, TRUE) && is_diagonal(D, tol_zero, TRUE)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    diminished <- FALSE
    if(use_vec) {
        LA <- diag(A)
        LD <- diag(D)
        if(missing(nthreads)) nthreads <- 1
    } else {
        eigD <- eigen(D, symmetric = TRUE)
        LA <- eigen(A, symmetric = TRUE)$values
        LD <- eigD$values
    }
    stopifnot("B should be nonnegative definite" = all(LB >= -tol_sing),
              "D should be nonnegative definite" = all(LD >= -tol_sing))
    ## Check condition for existence of moment
    nzB <- (LB > tol_sing)
    nzD <- (LD > tol_sing)
    if(use_vec) {
        ## common nonzero space of B and D, and A "rotated" with its basis
        nzBD <- nzB * nzD
        Ar <- A
    } else {
        if(all(nzB) && all(nzD)) {
            nzBD <- rep.int(TRUE, n)
            Ar <- A
        } else {
            zerocols <- cbind(In[, !nzB], eigD$vectors[, !nzD])
            projmat <- zerocols %*% MASS::ginv(crossprod(zerocols)) %*%
                       t(zerocols)
            eigBD <- eigen(In - projmat, symmetric = TRUE)
            nzBD <- eigBD$values > tol_sing
            Ar <- with(eigBD, crossprod(crossprod(A, vectors), vectors))
        }
    }
    rBD <- sum(nzBD)
    if(rBD == n) {
        cond_exist <- n / 2 + p > q + r ## condition(1)
        necess_cond <- TRUE
    } else {
        A12z <- all(abs(Ar[nzBD, !nzBD]) < tol_zero)
        A22z <- all(abs(Ar[!nzBD, !nzBD]) < tol_zero)
        cond_exist <- if(!A22z) {
                    rBD / 2 > q + r              ## condition(2)(iii)
                } else {
                    if(!A12z) {
                        (rBD + p) / 2 > q + r    ## condition(2)(ii)
                    } else {
                        rBD / 2 + p > q + r      ## condiiton(2)(i)
                    }
                }
        necess_cond <- (rBD == sum(nzB)) || (rBD == sum(nzD))
    }
    if(!cond_exist) {
        if(necess_cond) {
            stop("Moment does not exist in this combination of p, q, r, and",
                 "\n  eigenstructures of A, B, and D")
        } else {
            warning("Moment may not exist in this combination of p, q, r, and",
                    "\n  eigenstructures of A, B, and D")
        }
    }
    if(any(LA < -tol_sing) && (p %% 1) != 0) {
        stop("Detected negative eigenvalue(s) of A (< -tol_sing), ",
             "with which\n  non-integer power of quadratic form is not ",
             "well defined.\n  If you know them to be 0, use larger tol_sing ",
             "to suppress this error.")
    }
    bA <- alphaA / max(abs(LA))
    bD <- alphaD / max(LD)
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_npi_cvEc(LA, LB, LD, bA, bB, bD,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_npi_cvEl(LA, LB, LD, bA, bB, bD,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_npi_cvE(LA, LB, LD, bA, bB, bD,
                                             p, q, r, m, thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_npi_cmEc(A, LB, D, bA, bB, bD,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_npi_cmEl(A, LB, D, bA, bB, bD,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_npi_cmE(A, LB, D, bA, bB, bD,
                                             p, q, r, m, thr_margin, nthreads)
                }
            }
        } else {
            if(use_vec) {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_npi_nvEc(LA, LB, LD, bA, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_npi_nvEl(LA, LB, LD, bA, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_npi_nvE(LA, LB, LD, bA, bB, bD, mu,
                                             p, q, r, m, thr_margin, nthreads)
                }
            } else {
                if(cpp_method == "coef_wise") {
                    cppres <- ApBDqr_npi_nmEc(A, LB, D, bA, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else if(cpp_method == "long_double") {
                    cppres <- ApBDqr_npi_nmEl(A, LB, D, bA, bB, bD, mu,
                                              p, q, r, m, thr_margin, nthreads)
                } else {
                    cppres <- ApBDqr_npi_nmE(A, LB, D, bA, bB, bD, mu,
                                             p, q, r, m, thr_margin, nthreads)
                }
            }
        }
        ansseq <- cppres$ansseq
        diminished <- cppres$diminished
    } else {
        if(use_vec) {
            LAh <- rep.int(1, n) - bA * LA
            LBh <- rep.int(1, n) - bB * LB
            LDh <- rep.int(1, n) - bD * LD
            if(central) {
                dksm <- d3_ijk_v(LAh, LBh, LDh, m, thr_margin = thr_margin)
            } else {
                dksm <- h3_ijk_v(LAh, LBh, LDh, mu, m, thr_margin = thr_margin)
            }
        } else {
            Ah <- In - bA * A
            Bh <- In - bB * diag(LB, nrow = n)
            Dh <- In - bD * D
            if(central) {
                dksm <- d3_ijk_m(Ah, Bh, Dh, m, thr_margin = thr_margin)
            } else {
                dksm <- h3_ijk_m(Ah, Bh, Dh, mu, m, thr_margin = thr_margin)
            }
        }
        lscf <- attr(dksm, "logscale")
        ansarr <- hgs_3d(dksm, -p, q, r, n / 2, ((p - q - r) * log(2)
                         - p * log(bA) + q * log(bB) + r * log(bD)
                         + lgamma(n/2 + p - q - r) - lgamma(n/2) - lscf))
        ansseq <- sum_counterdiag3D(ansarr)
        if(any(lscf < 0)) {
            for(k in 1:(m + 1)) {
                diminished <-
                    any(diag(dksm[(m + 2 - k):1, 1:(m + 2 - k), k]) == 0)
                if(diminished) break
            }
        }
    }
    ## If there's any NaN, truncate series before summing up
    nans_ansseq <- is.nan(ansseq) | is.infinite(ansseq)
    if(any(nans_ansseq)) {
        warning("NaNs detected at k = ", which(nans_ansseq)[1L],
                if(sum(nans_ansseq) > 1) " ..." else NULL,
                "\n  The result is truncated before the first NaN")
        ansseq <- ansseq[-(which(nans_ansseq)[1L]:length(ansseq))]
        # m <- length(ansseq) - 1L
        attr(ansseq, "truncated") <- TRUE
    }
    if(diminished) {
        warning("Some terms in multiple series numerically diminished to 0 ",
                "as they were\n  scaled to avoid numerical overflow. ",
                "The result will be inaccurate.",
                if(cpp_method != "coef_wise")
                    paste0("\n  Consider using the option cpp_method = ",
                          if(cpp_method != "long_double") "\"long_double\" or ",
                          "\"coef_wise\"."))
    }
    if(check_convergence != "none") {
        if(missing(tol_conv) && check_convergence == "strict_relative") {
            tol_conv <- .Machine$double.eps
        }
        non_convergence <- if(check_convergence == "absolute") {
            abs(ansseq[length(ansseq)]) > tol_conv
        } else {
            abs(ansseq[length(ansseq)] / sum(ansseq)) > tol_conv
        }
        if(non_convergence) {
            warning("Last term of series is more than ",
                    sprintf("%.1e", tol_conv),
                    if(check_convergence != "absolute")
                        " times as large as the sum",
                    ",\n  suggesting non-convergence. Consider using larger m.")
        }
    }
    new_qfrm(terms = ansseq, seq_error = NA_real_)
}
