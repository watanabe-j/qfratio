## Documentation of the package
#' qfratio: Moments of Ratios of Quadratic Forms Using Recursion
#'
#' This package provides functions to evaluate moments of ratios (and products)
#' of quadratic forms of normal variables, specifically using recursive
#' algorithms developed by Hillier, Bao, Kan and colleagues.
#'
#' The primary front-end functions of this package are
#' \code{\link{qfrm}()} and \code{\link{qfmrm}()} for evaluating moments of
#' ratios of quadratic forms.
#' These pass arguments to one of the several ``internal'' (though exported)
#' functions which do actual calculations, depending on the argument matrices
#' and exponents.
#' In addition, there are a few functions to calculate moments of products
#' of quadratic forms (integer exponents only; \code{\link{qfpm}}).
#'
#' There are many internal functions for calculating coefficients in
#' power-series expansion of generating functions for these moments
#' (\code{\link{d1_i}}, \code{\link{d2_ij}}, \code{\link{d3_ijk}},
#' \code{\link{dtil2_pq}}) using ``super-short'' recursions
#' (Bao & Kan, 2013; Hillier et al. 2014).
#' Some of these coefficients are related to the top-order zonal and invariant
#' polynomials of matrix arguments.
#'
#' The DESCRIPTION file:
#' \packageDESCRIPTION{qfratio}
#' \packageIndices{qfratio}
#'
#' @section Author/Maintainer:
#' Junya Watanabe <jw2098@cam.ac.uk>
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
#' Watanabe, J. (2022). Exact expressions and numerical evaluation of average
#'   evolvability measures for characterizing and comparing **G** matrices.
#'   *bioRxiv* preprint, 2022.11.02.514929.
#'   doi:[10.1101/2022.11.02.514929](https://doi.org/10.1101/2022.11.02.514929).
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
# #' @aliases qfratio
#'
## Setting Rcpp-related import
#' @useDynLib qfratio, .registration = TRUE
# #' @import RcppEigen
#' @importFrom Rcpp evalCpp
#'
# #' @importFrom Rcpp sourceCpp
# #' @exportPattern "^[[:alpha:]]+"
NULL

## Documentation of C++ functions
#' Internal C++ functions
#'
#' These are internal \code{C++} functions called from
#' corresponding \code{R} functions when \code{use_cpp = TRUE}.
#' Direct access by the user is **not** assumed.
#' All parameters are assumed to be appropriately structured.
#'
#' At present, \code{ApIq_int_nmE()} calls the \code{R} function
#' \code{gsl::hyperg_1F1()}, so will not be much faster than
#' the \code{R} equivalent.
#' Ideally, the \code{C++} library \code{gsl} (or the like) should be used with
#' \code{RcppGSL}, but this is not done (commented out in the source code)
#' to increase portability.
#'
#' \code{rqfpE} uses \code{Rcpp::rnorm()},
#' which may not be particularly efficient.
#'
#' @param A,B,D
#'   Argument matrices passed as \code{Eigen::MatrixXd}.
#'   Symmetry is assumed.
#' @param LA,LB,LD
#'   Eigenvalues of the argument matrices passed as \code{Eigen::ArrayXd}
#' @param UA
#'   Matrix whose columns are eigenvectors of \eqn{\mathbf{A}} corresponding to
#'   \code{LA}.  Passed as \code{Eigen::MatrixXd}.
#' @param bA,bB,bD
#'   Scaling coefficients for \eqn{\mathbf{A}}, \eqn{\mathbf{B}},
#'   and \eqn{\mathbf{D}}.  Passed as \code{double}.
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}
#'   passed as \code{Eigen::ArrayXd}
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}} for \eqn{\mathbf{x}}.
#'   Passed as \code{Eigen::MatrixXd}.
#' @param p,q,r
#'   Exponents for \eqn{\mathbf{A}}, \eqn{\mathbf{B}}, and \eqn{\mathbf{D}}.
#'   Passed as \code{double}.
#' @param m
#'   \code{int} to specify the order of polynomials at which the series
#'   expression is truncated
#' @param error_bound
#'   \code{bool} to specify whether the error bound is returned
#' @param nthreads
#'   \code{int} to specify the number of threads in \code{OpenMP}-enabled
#'   functions.  See "Details" in \code{\link{qfmrm}}.
#' @param thr_margin
#'   Optional argument to adjust the threshold for scaling. See "Scaling"
#'   in \code{\link{d1_i}}.
#' @param nit
#'   \code{int} to specify the number of iteration or sample size
#'
#' @return
#'   All return a list via \code{Rcpp::List} of the following (as appropriate):
#'   \itemize{
#'      \item{\code{$ans}: }{Exact moment, from \code{double}}
#'      \item{\code{$ansseq}: }{Series for the moment, from
#'                              \code{Eigen::ArrayXd}}
#'      \item{\code{$errseq}: }{Series of errors, from \code{Eigen::ArrayXd}}
#'      \item{\code{$twosided}: }{Logical, from \code{bool}}
#'   }
#'
#' @name qfrm_cpp
#'
NULL
