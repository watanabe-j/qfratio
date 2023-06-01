## Documentation of the package
#' qfratio: Moments of Ratios of Quadratic Forms Using Recursion
#'
#' This package is for evaluating moments of ratios (and products) of quadratic
#' forms in normal variables, specifically using recursive algorithms developed
#' by Bao et al. (2013) and Hillier et al. (2014) (see also Smith, 1989, 1993;
#' Hillier et al., 2009).  It was originally developed as a supplement to
#' Watanabe (2023) for evaluating average evolvability measures in evolutionary
#' quantitative genetics, but can be used for a broader class of moments.
#'
#' The primary front-end functions of this package are
#' \code{\link{qfrm}()} and \code{\link{qfmrm}()} for evaluating moments of
#' ratios of quadratic forms.  These pass arguments to one of the several
#' \dQuote{internal} (though exported) functions which do actual calculations,
#' depending on the argument matrices and exponents.  In addition, there are
#' a few functions to calculate moments of products
#' of quadratic forms (integer exponents only; \code{\link{qfpm}}).
#'
#' There are many internal functions for calculating coefficients in
#' power-series expansion of generating functions for these moments
#' (\code{\link{d1_i}}, \code{\link{d2_ij}}, \code{\link{d3_ijk}},
#' \code{\link{dtil2_pq}}) using \dQuote{super-short} recursions
#' (Bao and Kan, 2013; Hillier et al. 2014).  Some of these coefficients are
#' related to the top-order zonal and invariant
#' polynomials of matrix arguments.
#'
#' See package vignette (\code{vignette("qfratio")}) for more details.
#'
#' The DESCRIPTION file:
#' \packageDESCRIPTION{qfratio}
#' \packageIndices{qfratio}
#'
#' @section Author/Maintainer:
#' Junya Watanabe <jw2098@cam.ac.uk>
#'
#' @references
#' Bao, Y. and Kan, R. (2013) On the moments of ratios of quadratic forms in
#'   normal random variables. *Journal of Multivariate Analysis*, **117**,
#'   229--245.
#'   \doi{10.1016/j.jmva.2013.03.002}.
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
#' Smith, M. D. (1989) On the expectation of a ratio of quadratic forms
#'   in normal variables. *Journal of Multivariate Analysis*, **31**, 244--257.
#'   \doi{10.1016/0047-259X(89)90065-1}.
#'
#' Smith, M. D. (1993) Expectations of ratios of quadratic forms in normal
#'   variables: evaluating some top-order invariant polynomials.
#'   *Australian Journal of Statistics*, **35**, 271--282.
#'   \doi{10.1111/j.1467-842X.1993.tb01335.x}.
#'
#' Watanabe, J. (2023) Exact expressions and numerical evaluation of average
#'   evolvability measures for characterizing and comparing **G** matrices.
#'   *Journal of Mathematical Biology*, **86**, 95.
#'   \doi{10.1007/s00285-023-01930-8}.
#'
#' @seealso
#'   \code{\link{qfrm}}: Moment of simple ratio of quadratic forms
#'
#'   \code{\link{qfmrm}}: Moment of multiple ratio of quadratic forms
#'
#'   \code{\link{qfpm}}: Moment of product of quadratic forms
#'
#'   \code{\link{rqfr}}: Monte Carlo sampling of ratio/product of
#'                       quadratic forms
#'
#' @examples
#' ## Symmetric matrices
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(sqrt(1:nv))
#' D <- diag((1:nv)^2 / nv)
#' mu <- nv:1 / nv
#' Sigma <- matrix(0.5, nv, nv)
#' diag(Sigma) <- 1
#'
#' ## Expectation of (x^T A x)^2 / (x^T x)^2 where x ~ N(0, I)
#' qfrm(A, p = 2)
#' ## And a Monte Carlo mean of the same
#' mean(rqfr(1000, A = A, p = 2))
#'
#' ## Expectation of (x^T A x)^1/2 / (x^T x)^1/2 where x ~ N(0, I)
#' (res1 <- qfrm(A, p = 1/2))
#' plot(res1)
#' ## A Monte Carlo mean
#' mean(rqfr(1000, A = A, p = 1/2))
#'
#' ## (x^T A x)^2 / (x^T B x)^3 where x ~ N(mu, Sigma)
#' (res2 <- qfrm(A, B, p = 2, q = 3, mu = mu, Sigma = Sigma))
#' plot(res2)
#' ## A Monte Carlo mean
#' mean(rqfr(1000, A = A, B = B, p = 2, q = 3, mu = mu, Sigma = Sigma))
#'
#' ## Expectation of (x^T A x)^2 / (x^T B x) (x^T x) where x ~ N(0, I)
#' (res3 <- qfmrm(A, B, p = 2, q = 1, r = 1))
#' plot(res3)
#' ## A Monte Carlo mean
#' mean(rqfmr(1000, A = A, B = B, p = 2, q = 1, r = 1))
#'
#' ## Expectation of (x^T A x)^2 where x ~ N(0, I)
#' qfm_Ap_int(A, 2)
#' ## A Monte Carlo mean
#' mean(rqfp(1000, A = A, p = 2, q = 0, r = 0))
#'
#' ## Expectation of (x^T A x) (x^T B x) (x^T D x) where x ~ N(mu, I)
#' qfpm_ABDpqr_int(A, B, D, 1, 1, 1, mu = mu)
#' ## A Monte Carlo mean
#' mean(rqfp(1000, A = A, B = B, D = D, p = 1, q = 1, r = 1, mu = mu))
#'
#' @docType package
#' @name qfratio-package
#'
## Setting Rcpp-related import
#' @useDynLib qfratio, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
NULL

## Documentation of C++ functions
#' Internal C++ functions
#'
#' These are internal \proglang{C++} functions called from corresponding \R
#' functions when \code{use_cpp = TRUE}.  Direct access by the user is **not**
#' assumed.  All parameters are assumed to be appropriately structured.
#'
#' \code{ApIq_int_nmE()} calls the \proglang{C} function
#' \code{gsl_sf_hyperg_1F1()} from \proglang{GSL} via \pkg{RcppGSL}.
#'
#' @param A,B,D
#'   Argument matrices passed as \code{Eigen::Matrix}.
#'   Symmetry is assumed.
#' @param LA,LB
#'   Eigenvalues of the argument matrices passed as \code{Eigen::Array}
#' @param bA,bB,bD
#'   Scaling coefficients for \eqn{\mathbf{A}}{A}, \eqn{\mathbf{B}}{B},
#'   and \eqn{\mathbf{D}}{D}.  Passed as \code{double} or \code{long double}.
#' @param mu
#'   Mean vector \eqn{\bm{\mu}}{\mu} for \eqn{\mathbf{x}}{x}
#'   passed as \code{Eigen::Array}
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}}{\Sigma} for \eqn{\mathbf{x}}{x}.
#'   Passed as \code{Eigen::Matrix}.
#' @param p_,q_,r_
#'   Exponents for \eqn{\mathbf{A}}{A}, \eqn{\mathbf{B}}{B}, and
#'   \eqn{\mathbf{D}}{D}.
#'   Passed as \code{double} or \code{long double}.
#' @param m
#'   Integer to specify the order of polynomials at which the series
#'   expression is truncated.  Passed as \code{Eigen::Index}
#'   (aka \code{std::ptrdiff_t} or \code{long long int})
#' @param error_bound
#'   \code{bool} to specify whether the error bound is returned
#' @param nthreads
#'   \code{int} to specify the number of threads in \proglang{OpenMP}-enabled
#'   functions.  See \dQuote{Multithreading} in \code{\link{qfrm}}.
#' @param thr_margin
#'   Optional argument to adjust the threshold for scaling.  See
#'   \dQuote{Scaling} in \code{\link{d1_i}}.
#' @param tol_zero
#'   Tolerance against which numerical zero is determined
#' @param nit
#'   \code{int} to specify the number of iteration or sample size
#' @param D1,D2
#'   Positive and (positivized) negative eigenvalues of
#'   \eqn{\mathbf{A} - q \mathbf{B}}{A - qB} passed as \code{Eigen::Array}
#' @param mu1,mu2
#'   Parts of mean vector \eqn{\bm{\mu}}{\mu} corresponding to
#'   \code{D1} and \code{D2}, respectively
#' @param quantile,f
#'   Scalar of quantile \eqn{q} and corresponding \eqn{f},
#'   passed as \code{double}
#' @param psi
#'   Vector of transformed eigenvalues \eqn{\psi_i},
#'   passed as \code{Eigen::Array}
#' @param L1,Ls
#'   Largest and smallest eigenvalues of \eqn{\mathbf{A}}{A},
#'   \eqn{\lambda_1} and \eqn{\lambda_s}, passed as \code{double}
#' @param n_,n1_,ns_
#'   Number of eigenvalues \eqn{n} and degrees of multiplicity of
#'   largest and smallest eigenvalues, \eqn{n_1} and \eqn{n_s},
#'   passed as \code{double}
#'
#' @return
#'   All return a list via \code{Rcpp::List} of the following (as appropriate):
#'   \itemize{
#'      \item{\code{$ans}: }{Exact moment, from \code{double} or
#'                           \code{long double}}
#'      \item{\code{$ansseq}: }{Series for the moment, from
#'                              \code{Eigen::Array}}
#'      \item{\code{$errseq}: }{Series of errors, from \code{Eigen::Array}}
#'      \item{\code{$twosided}: }{Logical, from \code{bool}}
#'      \item{\code{$dimnished}: }{Logical, from \code{bool}}
#'   }
#'
#' @name qfrm_cpp
#'
NULL
