## Documentation of the package
#' qfratio: Moments of Ratios of Quadratic Forms Using Recursion
#'
#' This package provides functions to evaluate moments of ratios
#' of quadratic forms of normal variables, specifically using recursive
#' algorithms developed by Chikuse, Hillier, Bao, Kan and colleagues.
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
#' @seealso
#'   \code{\link{qfratio}}: Moment of simple ratio of quadratic forms
#'
#' @examples
#' ## Symmetric matrices
#' nv <- 4
#' A <- diag(nv:1)
#' B <- diag(sqrt(1:nv))
#'
#' ## Expectation of (x^T A x)^2 / (x^T x)^2 where x ~ N(0, I)
#' qfrm(A, p = 2)
#'
#' ## Expectation of (x^T A x)^1/2 / (x^T x)^1/2 where x ~ N(0, I)
#' qfrm(A, p = 1/2)
#'
#' ## (x^T A x)^2 / (x^T B x)^3 where x ~ N(0, I)
#' qfrm(A, B, p = 2, q = 3)
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
