##### d1_i (documentation) #####
#' Coefficients in polynomial expansion of generating function---single matrix
#'
#' These are internal functions to calculate the coefficients
#' in polynomial expansion of generating functions for quadratic forms
#' in multivariate normal variables.
#'
#' \code{d1_i()} calculates \eqn{d_k(\mathbf{A})}, and \code{dtil1_i_v()} and
#' \code{dtil1_i_m()} calculate \eqn{\tilde{d}_k(\mathbf{A})} in
#' Hillier et al. (2009, 2014) and Bao & Kan (2013).
#' The former is related to the top-order zonal polynomial
#' \eqn{C_{[k]}(\mathbf{A})} in the following way:
#' \eqn{ d_k(\mathbf{A}) = \frac{1}{k!} \left( \frac{1}{2} \right)_k C_{[k]}(\mathbf{A}) },
#' where \eqn{(x)_k = x (x + 1) \dots (x + k - 1)}.
#'
#' These functions calculate the coefficients based on the super-short
#' recursion algorithm described in Hillier et al. (2014: 3.2, eqs. 28--30).
#'
#' ## Scaling:
#' The coefficients described herein (and in \code{\link{d2_ij}} and
#' \code{\link{d3_ijk}}) can become very large for higher-order terms,
#' so there is a practical risk of numerical overflow when applied to
#' large (\eqn{n > 100}, say) matrices.
#' To avoid numerical overflow, these functions automatically scale
#' coefficients (and temporary objects used to calculate them) by a large number
#' (\code{1e10} at present) when any value in the temporary objects
#' exceeds a predefined threshold (\code{.Machine$double.xmax / 100 / n}).
#' This scaling happens order-wise; i.e., it influences all the coefficients
#' of the same order in multidimensional coefficients (in \code{\link{d2_ij}}
#' and \code{\link{d3_ijk}}) and the coefficients of the subsequent orders.
#'
#' These scaling factors are recorded in the attribute \code{"logscale"} of the
#' return value, which is a vector/matrix/array whose size is identical to the
#' return value, so that \code{value / exp(attr(value, "logscale"))} equals
#' the original quantities to be obtained (if there were no overflow).
#'
#' The \code{qfrm} and \code{qfmrm} functions handle return values of these
#' functions by first multiplying them with hypergeometric coefficients
#' (which are typically \eqn{\ll 1}) and then scaling the products back
#' to the original scale using the recorded scaling factors.
#' (To be precise, this typically happens within \code{\link{hgs}} functions.)
#' The \code{C++} functions handle the problem similarly (but by using
#' separate objects rather than attributes).
#'
#' However, this procedure does not always mitigate the problem in
#' multiple series; when there are very large and very small
#' coefficients in the same order, smaller ones can diminish/underflow to
#' the numerical \code{0} after repeated scaling.
#' (The \code{qfrm} and \code{qfmrm} functions try to detect and warn against
#' this by examining whether any of the highest-order terms is \code{0}.)
#' As this seems technically difficult to avoid without implementing
#' cumbersome and inefficient coefficient-wise scaling, the only workaround
#' implemented in this package is to use the \code{long double} variable type
#' in \code{C++} versions (see \code{\link{qfrm}} and \code{\link{qfmrm}}).
#'
#' @return
#' Vector of length \code{m + 1}, corresponding to
#' the 0th, 1st, ..., and mth order terms.
#' Hence, the \code{[m + 1]}-th element should be extracted
#' when the coefficient for the mth order term is required.
#'
#' Has the attribute \code{"logscale"} as described in "Scaling" above.
#'
#' @param L
#'   Vector of eigenvalues of the argument matrix
#' @param m
#'   Integer-alike to specify the order of polynomials
#' @param A
#'   Argument matrix.  Assumed to be symmetric in these functions.
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}
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
#' @seealso
#' \code{\link{qfpm}}, \code{\link{qfrm}}, and \code{\link{qfmrm}} are
#' major front-end functions that utilize these functions
#'
#' \code{\link{dtil2_pq}} for \eqn{\tilde{d}}
#' used for moments of a product of quadratic forms
#'
#' \code{\link{d2_ij}} and \code{\link{d3_ijk}} for \eqn{d}, \eqn{h},
#' \eqn{\tilde{h}}, and \eqn{\hat{h}} used for moments of ratios
#' of quadratic forms
#'
#' @name d1_i
#' @aliases dtil1_i
#'
NULL

##### dtil2_pq (documentation) #####
#' Coefficients in polynomial expansion of generating function---for products
#'
#' These are internal functions to calculate the coefficients
#' in polynomial expansion of joint generating functions for two or three
#' quadratic forms in potentially noncentral multivariate normal variables,
#' \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{I}_n)}.
#' They are primarily used to calculate moments of a product of two or
#' three quadratic forms.
#'
#' \code{dtil2_pq_m()} and \code{dtil2_pq_v()} calculate
#' \eqn{\tilde{d}_{p,q}(\mathbf{A}_1, \mathbf{A}_2)} in
#' Hillier et al. (2014).
#' \code{dtil2_1q_m()} and \code{dtil2_1q_v()} are fast versions
#' for the commonly used case where \eqn{p = 1}.
#' Similarly, \code{dtil3_pqr_m()} and \code{dtil3_pqr_v()} are for
#' \eqn{\tilde{d}_{p,q,r}(\mathbf{A}_1, \mathbf{A}_2, \mathbf{A}_3)}.
#'
#' Those ending with \code{_m} take matrices as arguments, whereas
#' those with \code{_v} take eigenvalues.
#' The latter can be used only when the argument matrices share the same
#' eigenvectors, to which the eigenvalues correspond in the orders given,
#' but is substantially faster.
#'
#' These functions calculate the coefficients based on the super-short
#' recursion algorithm described in Hillier et al. (2014: sec. 4.2).
#'
#' @return
#' A \code{(p + 1) * (q + 1)} matrix for the 2D functions,
#' or a \code{(p + 1) * (q + 1) * (r + 1)} array for the 3D functions.
#'
#' The 1st, 2nd, and 3rd dimensions correspond to increasing orders for
#' \eqn{\mathbf{A}_1}, \eqn{\mathbf{A}_2}, and \eqn{\mathbf{A}_3},
#' respectively. And the 1st row/column of each dimension corresponds
#' to the 0th order (hence \code{[p + 1, q + 1]} for the \eqn{(p,q)}-th moment).
#'
#' @param A1,A2,A3
#'   Argument matrices. Assumed to be symmetric and of the same order.
#' @param L1,L2,L3
#'   Eigenvalues of the argument matrices
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}
#' @param p,q,r
#'   Integer-alikes to specify the order along the three argument matrices
#'
#' @references
#' Hillier, G., Kan, R, & Wang, X. (2014). Generating functions and
#'   short recursions, with applications to the moments of quadratic forms
#'   in noncentral normal vectors. *Econometric Theory*, **30**, 436--473.
#'   doi:[10.1017/S0266466613000364](https://doi.org/10.1017/S0266466613000364).
#'
#' @seealso
#' \code{\link{qfpm}} is a front-end functions that utilizes these functions
#'
#' \code{\link{d1_i}} for a single-matrix equivalent of \eqn{\tilde{d}}
#'
#' @name dtil2_pq
#' @aliases dtil3_pqr
#'
NULL

##### d2_ij (documentation) #####
#' Coefficients in polynomial expansion of generating function---for
#' ratios with two matrices
#'
#' These are internal functions to calculate the coefficients
#' in polynomial expansion of joint generating functions for two
#' quadratic forms in potentially noncentral multivariate normal variables,
#' \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{I}_n)}.
#' They are primarily used in calculations around moments of a ratio
#' involving two or three quadratic forms.
#'
#' \code{d2_**_*()} functions calculate
#' \eqn{d_{i,j}(\mathbf{A}_1, \mathbf{A}_2)} in
#' Hillier et al. (2009, 2014) and Bao & Kan (2013).
#' These are also related to the top-order invariant polynomials
#' \eqn{C_{[k_1],[k_2]}(\mathbf{A}_1, \mathbf{A}_2)} in the following way:
#' \eqn{ d_{i,j}(\mathbf{A}_1, \mathbf{A}_2) =
#'      \frac{1}{k_1! k_2!} \left( \frac{1}{2} \right)_{k_1 + k_2}
#'      C_{[k_1],[k_2]}(\mathbf{A}_1, \mathbf{A}_2) },
#' where \eqn{(x)_k = x (x + 1) \dots (x + k - 1)}
#' (Chikuse 1987; Hillier et al. 2009).
#'
#' \code{h2_ij_*()} and \code{htil2_pj_*()} functions calculate
#' \eqn{h_{i,j}(\mathbf{A}_1, \mathbf{A}_2)} and
#' \eqn{\tilde{h}_{i,j}(\mathbf{A}_1, \mathbf{A}_2)}, respectively,
#' in Bao & Kan (2013). Note that the latter is denoted by the symbol
#' \eqn{h_{i,j}} in Hillier et al. (2014).
#' \code{hhat2_pj_*()} functions are for
#' \eqn{\hat{h}_{i,j}(\mathbf{A}_1, \mathbf{A}_2)}
#' in Hillier et al. (2014), used to calculate an error bound for
#' truncated sum for moments of a ratio of quadratic forms.
#' The mean vector \eqn{\bm{\mu}} is a parameter in all these.
#'
#' There are two different situations in which these coefficients are used
#' in calculation of moments of ratios of quadratic forms:
#' **1**) within an infinite series for one of the subscripts, with the
#' other subscript fixed (when the exponent \eqn{p} of the numerator
#' is integer); **2**) within a double infinite series for both subscripts
#' (when \eqn{p} is non-integer) (see Bao & Kan 2013).
#' In this package, the situation **1** is handled by
#' the \code{*_pj_*} (and \code{*_1j_*}) functions, and **2** is by
#' the \code{*_ij_*} functions.
#'
#' In particular, the \code{*_pj_*} functions always return a
#' \code{(p + 1) * (m + 1)} matrix where all elements are filled with
#' the relevant coefficients (e.g., \eqn{d_{i,j}}, \eqn{\tilde{h}_{i,j}}),
#' from which, typically, the \code{[p + 1, ]}-th row is used for
#' subsequent calculations.
#' (Those with \code{*_1q_*} are simply fast versions
#' for the commonly used case where \eqn{p = 1}.)
#' On the other hand, the \code{*_ij_*} functions by default return a
#' \code{(m + 1) * (m + 1)} matrix whose upper-left triangular part
#' (including the diagonals) is filled with the coefficients
#' (\eqn{d_{i,j}} or \eqn{h_{i,j}}), the rest being 0, and all the coefficients
#' are used in subsequent calculations.
#'
#' (At present, the \code{*_ij_*} functions also have the functionality to
#' fill all coefficients of a potentially non-square output matrix,
#' but this is less efficient than \code{*_pj_*} functions so may
#' be omitted in the future development.)
#'
#' Those ending with \code{_m} take matrices as arguments, whereas
#' those with \code{_v} take eigenvalues.
#' The latter can be used only when the argument matrices share the same
#' eigenvectors, to which the eigenvalues correspond in the orders given,
#' but is substantially faster.
#'
#' This package also involves \code{C++} equivalents for most of these functions
#' (which are suffixed by \code{E} for \code{Eigen}),
#' but these are exclusively for internal use and not exposed to the user.
#'
#' These functions calculate the coefficients based on the super-short
#' recursion algorithm described in Hillier et al. (2014: sec. 5.4) abd
#' Bao & Kan (2014: sec. 5).
#' The algorithm for \eqn{\hat{h}_{i,j}} was said to be ``very similar'' to
#' that of \eqn{\tilde{h}_{i,j}} by Hillier et al. (2014), but differs
#' in the signs of some terms.
#'
#' @return
#' \code{(p + 1) * (m + 1)} matrix for the \code{*_pj_*} functions.
#'
#' \code{(m + 1) * (m + 1)} matrix for the \code{*_ij_*} functions.
#'
#' The rows and columns correspond to increasing orders for
#' \eqn{\mathbf{A}_1} and \eqn{\mathbf{A}_2}, respectively.
#' And the 1st row/column of each dimension corresponds
#' to the 0th order (hence \code{[p + 1, q + 1]} for the \eqn{(p,q)}-th order).
#'
#' Has the attribute \code{"logscale"} as described in the "Scaling" section
#' in \code{\link{d1_i}}. This is a matrix of the same size as the return itself.
#'
#' @param A1,A2
#'   Argument matrices. Assumed to be symmetric and of the same order.
#' @param L1,L2
#'   Eigenvalues of the argument matrices
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}
#' @param m
#'   Integer-alike to specify the desired order along \code{A2}/\code{L2}
#' @param p,q
#'   Integer-alikes to specify the desired orders along
#'   \code{A1}/\code{L1} and \code{A2}/\code{L2}, respectively.
#' @param fill_all
#'   Logical to specify whether all the output matrix should be filled.
#'   See "Details".
#'
#' @references
#' Bao, Y. & Kan, R. (2013). On the moments of ratios of quadratic forms in
#'   normal random variables. *Journal of Multivariate Analysis*, **117**,
#'   229--245.
#'   doi:[10.1016/j.jmva.2013.03.002](https://doi.org/10.1016/j.jmva.2013.03.002).
#'
#' Chikuse, Y. (1987). Methods for constructing top order invariant polynomials.
#'   *Econometric Theory*, **3**, 195--207.
#'   doi:[10.1017/S026646660001029X](https://doi.org/10.1017/S026646660001029X).
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
#' @seealso
#' \code{\link{qfrm}} and \code{\link{qfmrm}} are
#' major front-end functions that utilize these functions
#'
#' \code{\link{dtil2_pq}} for \eqn{\tilde{d}}
#' used for moments of a product of quadratic forms
#'
#' \code{\link{d3_ijk}} for equivalents for three matrices
#'
#' @name d2_ij
#' @aliases d2_pj h2_ij htil2_pj hhat2_pj d2_1j htil2_1j hhat2_1j
#'
NULL

##### d3_ijk (documentation) #####
#' Coefficients in polynomial expansion of generating function---for
#' ratios with three matrices
#'
#' These are internal functions to calculate the coefficients
#' in polynomial expansion of joint generating functions for three
#' quadratic forms in potentially noncentral multivariate normal variables,
#' \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{I}_n)}.
#' They are primarily used in calculations around moments of a ratio
#' involving three quadratic forms.
#'
#' All these functions have equivalents for two matrices (\code{\link{d2_ij}}),
#' to which the user is referred for documentations.
#' The primary difference of these functions from the latter is
#' the addition of arguments for the third matrix \code{A3}/\code{L3}.
#'
#' \code{d3_*jk_*()} functions calculate
#' \eqn{d_{i,j,k}(\mathbf{A}_1, \mathbf{A}_2, \mathbf{A}_3)} in
#' Hillier et al. (2009, 2014) and Bao & Kan (2013).
#' These are also related to the top-order invariant polynomials as described
#' in \code{\link{d2_ij}}.
#'
#' \code{h3_ijk_*()}, \code{htil3_pjk_*()}, and \code{hhat3_pjk_*()} functions
#'  calculate \eqn{h_{i,j,k}(\mathbf{A}_1, \mathbf{A}_2, \mathbf{A}_3)},
#' \eqn{\tilde{h}_{i;j,k}(\mathbf{A}_1; \mathbf{A}_2, \mathbf{A}_3)}, and
#' \eqn{\hat{h}_{i;j,k}(\mathbf{A}_1; \mathbf{A}_2, \mathbf{A}_3)},
#' respectively, as described in Watanabe (forthcoming). These are equivalent
#' to similar coefficients described in Bao & Kan (2013) and
#' Hillier et al. (2014).
#'
#' The difference between the \code{*_pjk_*} and \code{*_ijk_*} functions
#' is as described for \code{*_pj_*} and \code{*_ij_*}
#' (see "Details" in \code{\link{d2_ij}}).
#' The only difference is that these functions return a 3D array.
#' In the \code{*_pjk_*} functions, all the slices along the first dimension
#' (i.e., \code{[i, , ]}) are an upper-left triangular matrix
#' like what the \code{*_ij_*} functions return in the 2D case;
#' in other words, the return has the coefficients for the terms that satisfy
#' \eqn{j + k \leq m} for all \eqn{i = 0, 1, \dots, p}.
#' Typically, the \code{[p + 1, , ]}-th slice is used for subsequent
#' calculations.
#' In the return of the \code{*_ijk_*} functions, only the triangular prism
#' close to the \code{[1, 1, 1]} is filled with coefficients, which
#' correspond to the terms satisfying \eqn{i + j + k \leq m}.
#'
#' (At present, the \code{*_ijk_*} functions have the functionality to
#' fill all coefficients of a potentially non-cubic output array,
#' but this is less efficient than \code{*_pjk_*} functions so may
#' be omitted in the future development.)
#'
#' This package also involves \code{C++} equivalents for most of these functions
#' (which are suffixed by \code{E} for \code{Eigen}),
#' but these are exclusively for internal use and not exposed to the user.
#'
#' These functions calculate the coefficients based on the super-short
#' recursion algorithm described in Hillier et al. (2014: sec. 5.4) abd
#' Bao & Kan (2014: sec. 5).
#'
#' @return
#'
#' \code{(p + 1) * (m + 1) * (m + 1)} array for the \code{*_pjk_*} functions
#'
#' \code{(m + 1) * (m + 1) * (m + 1)} array for the \code{*_ijk_*} functions
#' (by default; see "Details").
#'
#' The 1st, 2nd, and 3rd dimensions correspond to increasing orders for
#' \eqn{\mathbf{A}_1}, \eqn{\mathbf{A}_2}, and \eqn{\mathbf{A}_3}, respectively.
#' And the 1st row/column of each dimension corresponds
#' to the 0th order (hence \code{[p + 1, q + 1, r + 1]} for
#' the \eqn{(p,q,r)}-th order).
#'
#' Has the attribute \code{"logscale"} as described in the "Scaling" section
#' in \code{\link{d1_i}}. This is an array of the same size as the return itself.
#'
#' @inheritParams d2_ij
#' @param A1,A2,A3
#'   Argument matrices. Assumed to be symmetric and of the same order.
#' @param L1,L2,L3
#'   Eigenvalues of the argument matrices
#' @param m
#'   Integer-alike to specify the desired order along \code{A2}/\code{L2}
#'   and \code{A3}/\code{L3}
#' @param p,q,r
#'   Integer-alikes to specify the desired orders along
#'   \code{A1}/\code{L1}, \code{A2}/\code{L2}, and \code{A3}/\code{L3},
#'   respectively.
#' @param fill_across
#'   Logical vector of length 3, to specify whether each dimension of
#'   the output matrix should be filled.
#'
#' @references
#' Bao, Y. & Kan, R. (2013). On the moments of ratios of quadratic forms in
#'   normal random variables. *Journal of Multivariate Analysis*, **117**,
#'   229--245.
#'   doi:[10.1016/j.jmva.2013.03.002](https://doi.org/10.1016/j.jmva.2013.03.002).
#'
#' Hillier, G., Kan, R, & Wang, X. (2014). Generating functions and
#'   short recursions, with applications to the moments of quadratic forms
#'   in noncentral normal vectors. *Econometric Theory*, **30**, 436--473.
#'   doi:[10.1017/S0266466613000364](https://doi.org/10.1017/S0266466613000364).
#'
#' @seealso
#' \code{\link{qfmrm}} is a
#' major front-end function that utilizes these functions
#'
#' \code{\link{dtil2_pq}} for \eqn{\tilde{d}}
#' used for moments of a product of quadratic forms
#'
#' \code{\link{d2_ij}} for equivalents for two matrices
#'
#' @name d3_ijk
#' @aliases d3_pjk h3_ijk htil3_pjk hhat3_pjk
#'
NULL

##### d1_i #####
#' Coefficients in polynomial expansion of generating function
#'
#' \code{d1_i()} is for standard multivariate normal variables,
#' \eqn{\mathbf{x} \sim N(\mathbf{0}, \mathbf{I}_n)}.
#'
#' @rdname d1_i
#'
d1_i <- function(L, m = 100L) {
    n <- length(L)
    dks <- rep.int(c(1, 0), c(1L, m))
    lscf <- rep.int(0, length(dks))
    thr <- .Machine$double.xmax / 100 / n
    uk <- rep.int(0, n)
    for(k in seq_len(m)) {
        uk <- L * (dks[k] + uk)
        dks[k + 1L] <- sum(uk) / (2 * k)
        if(max(uk) > thr) {
            dks[k + 1L] <- dks[k + 1L] / 1e10
            uk <- uk / 1e10
            lscf[(k + 1L):length(lscf)] <- lscf[(k + 1L):length(lscf)] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### dtil1_i_v #####
#' Coefficients in polynomial expansion of generating function
#'
#' \code{dtil1_i_v()} is for noncentral multivariate normal variables,
#' \eqn{\mathbf{x} \sim N(\bm{\mu}, \mathbf{I}_n)}.
#'
#' @rdname d1_i
#'
dtil1_i_v <- function(L, mu = rep.int(0, n), m = 100L) {
    n <- length(L)
    D <- mu ^ 2
    dks <- rep.int(c(1, 0), c(1L, m))
    lscf <- rep.int(0, length(dks))
    thr <- .Machine$double.xmax / 100 / n
    uk <- rep.int(0, n)
    vk <- rep.int(0, n)
    for(k in seq_len(m)) {
        uk <- L * (dks[k] + uk)
        vk <- D * uk + L * vk
        dks[k + 1L] <- sum(uk + vk) / (2 * k)
        if(max(uk) > thr || max(vk) > thr) {
            dks[k + 1L] <- dks[k + 1L] / 1e10
            uk <- uk / 1e10
            vk <- vk / 1e10
            lscf[(k + 1L):length(lscf)] <- lscf[(k + 1L):length(lscf)] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### dtil1_i_m #####
#' Coefficients in polynomial expansion of generating function
#'
#' \code{dtil1_i_m()} is a wrapper for \code{dtil1_i_v()}
#' and takes the argument matrix rather than its eigenvalues.
#'
#' @rdname d1_i
#'
dtil1_i_m <- function(A, mu = rep.int(0, n), m = 100L) {
    eigA <- eigen(A, symmetric = TRUE)
    L <- eigA$values
    D <- c(crossprod(eigA$vectors, mu)) ^ 2
    n <- length(L)
    dks <- rep.int(c(1, 0), c(1L, m))
    lscf <- rep.int(0, length(dks))
    thr <- .Machine$double.xmax / 100 / n
    uk <- rep.int(0, n)
    vk <- rep.int(0, n)
    for(k in seq_len(m)) {
        uk <- L * (dks[k] + uk)
        vk <- D * uk + L * vk
        dks[k + 1L] <- sum(uk + vk) / (2 * k)
        if(max(uk) > thr || max(vk) > thr) {
            dks[k + 1L] <- dks[k + 1L] / 1e10
            uk <- uk / 1e10
            vk <- vk / 1e10
            lscf[(k + 1L):length(lscf)] <- lscf[(k + 1L):length(lscf)] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

# dtil2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, p = m, q = m,
#                          fill_all = !missing(p) || !missing(q)) {
#     il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     zerovec <- matrix(0, n, 1)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     # attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeromat)
#     hs <- list(zerovec)
#     order_mat <- outer(0:p, 0:q, "+")
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         Gs[which(order_mat == k - 2L)] <- 0
#         hs[which(order_mat == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeromat
#             if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
#             if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
#             Gs[[il2(i1, i2)]] <- tG
#             ht <- tG %*% mu
#             if(i1 >= 1L) ht <- ht + A1 %*% hs[[il2(i1 - 1L, i2)]]
#             if(i2 >= 1L) ht <- ht + A2 %*% hs[[il2(i1, i2 - 1L)]]
#             hs[[il2(i1, i2)]] <- ht
#             dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, ht))) / (2 * k)
#         }
#         # if(max(unlist(Gs)) > thr || max(unlist(hs)) > thr) {
#         #     dks <- dks / 1e10
#         #     Gs <- lapply(Gs, function(x) x / 1e10)
#         #     hs <- lapply(hs, function(x) x / 1e10)
#         #     attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         # }
#     }
#     return(dks)
# }

# dtil2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, p = m, q = m,
#                          fill_all = !missing(p) || !missing(q)) {
#     il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     # attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeros)
#     hs <- list(zeros)
#     order_mat <- outer(0:p, 0:q, "+")
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         Gs[which(order_mat == k - 2L)] <- 0
#         hs[which(order_mat == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeros
#             if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
#             if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
#             Gs[[il2(i1, i2)]] <- tG
#             ht <- tG * mu
#             if(i1 >= 1L) ht <- ht + L1 * hs[[il2(i1 - 1L, i2)]]
#             if(i2 >= 1L) ht <- ht + L2 * hs[[il2(i1, i2 - 1L)]]
#             hs[[il2(i1, i2)]] <- ht
#             dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * ht)) / (2 * k)
#         }
#         # if(max(unlist(Gs)) > thr || max(ht) > thr) {
#         #     dks <- dks / 1e10
#         #     Gs <- lapply(Gs, function(x) x / 1e10)
#         #     hs <- lapply(hs, function(x) x / 1e10)
#         #     attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         # }
#     }
#     return(dks)
# }


# dtil3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
#                  fill_across = c(!missing(p), !missing(q), !missing(r))) { # , verbose = m > 200L) {
#     il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     zerovec <- matrix(0, n, 1)
#     dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
#     dks[1L] <- 1
#     # attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeromat)
#     hs <- list(zerovec)
#     order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
#     kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
#     for(k in 1L:kmax) {
#         Gs[which(order_array == k - 2L)] <- 0
#         hs[which(order_array == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
#             for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
#                 i3 <- k - i1 - i2
#                 if(fill_across[3L] && i1 + i2 + r > kmax) next
#                 if(fill_across[2L] && i1 + i3 + q > kmax) next
#                 # if(fill_across[1] && i2 + i3 + p > kmax) next
#                 tG <- zeromat
#                 if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
#                 if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
#                 if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
#                 Gs[[il3(i1, i2, i3)]] <- tG
#                 th <- tG %*% mu
#                 if(i1 >= 1L) th <- th + A1 %*% hs[[il3(i1 - 1L, i2, i3)]]
#                 if(i2 >= 1L) th <- th + A2 %*% hs[[il3(i1, i2 - 1L, i3)]]
#                 if(i3 >= 1L) th <- th + A3 %*% hs[[il3(i1, i2, i3 - 1L)]]
#                 hs[[il3(i1, i2, i3)]] <- th
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, th))) / (2 * k)
#             }
#         }
#         # if(max(unlist(Gs)) > thr || max(unlist(hs)) > thr) {
#         #     dks <- dks / 1e10
#         #     Gs <- lapply(Gs, function(x) x / 1e10)
#         #     hs <- lapply(hs, function(x) x / 1e10)
#         #     attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         # }
#         # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
#     }
#     return(dks)
# }

# dtil3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
#                  fill_across = c(!missing(p), !missing(q), !missing(r))) { # , verbose = m > 200L) {
#     il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
#     dks[1L] <- 1
#     # attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeros)
#     hs <- list(zeros)
#     order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
#     kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
#     for(k in 1L:kmax) {
#         Gs[which(order_array == k - 2L)] <- 0
#         hs[which(order_array == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
#             for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
#                 i3 <- k - i1 - i2
#                 if(fill_across[3L] && i1 + i2 + r > kmax) next
#                 if(fill_across[2L] && i1 + i3 + q > kmax) next
#                 # if(fill_across[1] && i2 + i3 + p > kmax) next
#                 tG <- zeros
#                 if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
#                 if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
#                 if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
#                 Gs[[il3(i1, i2, i3)]] <- tG
#                 th <- tG * mu
#                 if(i1 >= 1L) th <- th + L1 * hs[[il3(i1 - 1L, i2, i3)]]
#                 if(i2 >= 1L) th <- th + L2 * hs[[il3(i1, i2 - 1L, i3)]]
#                 if(i3 >= 1L) th <- th + L3 * hs[[il3(i1, i2, i3 - 1L)]]
#                 hs[[il3(i1, i2, i3)]] <- th
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * th)) / (2 * k)
#             }
#         }
#         # if(max(unlist(Gs)) > thr || max(unlist(hs)) > thr) {
#         #     dks <- dks / 1e10
#         #     Gs <- lapply(Gs, function(x) x / 1e10)
#         #     hs <- lapply(hs, function(x) x / 1e10)
#         #     attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         # }
#         # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
#     }
#     return(dks)
# }

##### dtil2_pq_m #####
#' Coefficients in polynomial expansion of generating function---for
#' product of two matrices
#'
#' @rdname dtil2_pq
#'
dtil2_pq_m <- function(A1, A2, mu = rep.int(0, n), p = 1L, q = 1L) {
    if(p == 1L) return(dtil2_1q_m(A1, A2, mu, q))
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, p + 1L, q + 1L)
    dks[1L, 1L] <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(p + 1L)] <- list(zeromat)
    g_k_i[1L:(p + 1L)] <- list(zerovec)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu + A1 %*% g_k_i[[i]]
        dks[i + 1L, 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * i)
    }
    for(k in 1L:q) {
        G_k_i[[1L]] <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        g_k_i[[1L]] <- G_k_i[[1L]] %*% mu + A2 %*% g_k_i[[1L]]
        dks[1L, k + 1L] <- (tr(G_k_i[[1L]]) + c(crossprod(mu, g_k_i[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            G_k_i[[i + 1L]] <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu +
                               A1 %*% g_k_i[[i]] + A2 %*% g_k_i[[i + 1L]]
            dks[i + 1L, k + 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * (k + i))
        }
    }
    return(dks)
}

##### dtil2_1q_m #####
#' Coefficients in polynomial expansion of generating function---for
#' product of two matrices
#'
#' @rdname dtil2_pq
#'
dtil2_1q_m <- function(A1, A2, mu = rep.int(0, n), q = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, 2L, q + 1L)
    dks[1L, 1L] <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1
    g_k_0 <- matrix(0, n, 1)
    g_k_1 <- G_k_1 %*% mu
    dks[2L, 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / 2
    for(k in 1L:q) {
        G_k_0 <- A2 %*% (dks[1L, k] * In + G_k_0)
        g_k_0 <- G_k_0 %*% mu + A2 %*% g_k_0
        dks[1L, k + 1L] <- (tr(G_k_0) + c(crossprod(mu, g_k_0))) / (2 * k)
        G_k_1 <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        g_k_1 <- G_k_1 %*% mu + A1 %*% g_k_0 + A2 %*% g_k_1
        dks[2L, k + 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / (2 * (k + 1))
    }
    return(dks)
}

##### dtil2_pq_v #####
#' Coefficients in polynomial expansion of generating function---for
#' product of two matrices
#'
#' @rdname dtil2_pq
#'
dtil2_pq_v <- function(L1, L2, mu = rep.int(0, n), p = 1L, q = 1L) {
    if(p == 1L) return(dtil2_1q_v(L1, L2, mu, q))
    n <- length(L1)
    dks <- matrix(0, p + 1L, q + 1L)
    dks[1L, 1L] <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(p + 1L)] <- list(zeros)
    g_k_i[1L:(p + 1L)] <- list(zeros)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu + L1 * g_k_i[[i]]
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:q) {
        G_k_i[[1L]] <- L2 * (dks[1L, k] + G_k_i[[1L]])
        g_k_i[[1L]] <- G_k_i[[1L]] * mu + L2 * g_k_i[[1L]]
        dks[1L, k + 1L] <- sum(G_k_i[[1L]] + mu * g_k_i[[1L]]) / (2 * k)
        for(i in 1L:p) {
            G_k_i[[i + 1L]] <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu +
                               L1 * g_k_i[[i]] + L2 * g_k_i[[i + 1L]]
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * (k + i))
        }
    }
    return(dks)
}

##### dtil2_1q_v #####
#' Coefficients in polynomial expansion of generating function---for
#' product of two matrices
#'
#' @rdname dtil2_pq
#'
dtil2_1q_v <- function(L1, L2, mu = rep.int(0, n), q = 1L) {
    n <- length(L1)
    dks <- matrix(0, 2L, q + 1L)
    dks[1L, 1L] <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1
    g_k_0 <- rep.int(0, n)
    g_k_1 <- G_k_1 * mu
    dks[2L, 1L] <- sum(G_k_1 + mu * g_k_1) / 2
    for(k in 1L:q) {
        G_k_0 <- L2 * (dks[1L, k] + G_k_0)
        g_k_0 <- G_k_0 * mu + L2 * g_k_0
        dks[1L, k + 1L] <- sum(G_k_0 + mu * g_k_0) / (2 * k)
        G_k_1 <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        g_k_1 <- G_k_1 * mu + L1 * g_k_0 + L2 * g_k_1
        dks[2L, k + 1L] <- sum(G_k_1 + mu * g_k_1) / (2 * (k + 1))
    }
    return(dks)
}


##### dtil3_pqr_m #####
#' Coefficients in polynomial expansion of generating function---for
#' product of three matrices
#'
#' @rdname dtil2_pq
#'
dtil3_pqr_m <- function(A1, A2, A3, mu = rep.int(0, n), p = 1L, q = 1L, r = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    m <- q + r
    dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
    dks[1L] <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    Gc <- list()
    gc <- list()
    Gc[1L:(p + 1L)] <- list(zeromat)
    gc[1L:(p + 1L)] <- list(zerovec)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu + A1 %*% gc[[i]]
        dks[i + 1L, 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        if(k <= q) {
            Gc[[1L]] <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
            gc[[1L]] <- Gc[[1L]] %*% mu + A2 %*% go[[1L]][[1L]]
            dks[1L, k + 1L, 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                                A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu +
                                   A1 %*% gc[[i]] + A2 %*% go[[1L]][[i + 1L]]
                dks[i + 1L, k + 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
            }
            Gn <- list(Gc)
            gn <- list(gc)
        } else {
            Gn <- list(0)
            gn <- list(0)
        }
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                if(k - j <= q && j <= r) {
                    Gc[[1L]] <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                    gc[[1L]] <- Gc[[1L]] %*% mu + A2 %*% go[[j + 1L]][[1L]] + A3 %*% go[[j]][[1L]]
                    dks[1L, k - j + 1L, j + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
                    for(i in 1L:p) {
                        Gc[[i + 1L]] <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                                        A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                                        A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu +
                                        A1 %*% gc[[i]] + A2 %*% go[[j + 1L]][[i + 1L]] + A3 %*% go[[j]][[i + 1L]]
                        dks[i + 1L, k - j + 1L, j + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
                    }
                    Gn <- c(Gn, list(Gc))
                    gn <- c(gn, list(gc))
                } else {
                    Gn <- c(Gn, list(0))
                    gn <- c(gn, list(0))
                }
            }
        }
        if(k <= r) {
            Gc[[1L]] <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
            gc[[1L]] <- Gc[[1L]] %*% mu + A3 %*% go[[k]][[1L]]
            dks[1L, 1L, k + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                                A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu +
                                A1 %*% gc[[i]] + A3 %*% go[[k]][[i + 1L]]
                dks[i + 1L, 1L, k + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
            }
            Gn <- c(Gn, list(Gc))
            gn <- c(gn, list(gc))
        } else {
            Gn <- c(Gn, list(0))
            gn <- c(gn, list(0))
        }
    }
    return(dks)
}

##### dtil3_pqr_v #####
#' Coefficients in polynomial expansion of generating function---for
#' product of three matrices
#'
#' @rdname dtil2_pq
#'
dtil3_pqr_v <- function(L1, L2, L3, mu = rep.int(0, n), p = 1L, q = 1L, r = 1L) {
    n <- length(L1)
    m <- q + r
    dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
    dks[1L] <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    Gc <- list()
    gc <- list()
    Gc[1L:(p + 1L)] <- list(zeros)
    gc[1L:(p + 1L)] <- list(zeros)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] * mu + L1 * gc[[i]]
        dks[i + 1L, 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        if(k <= q) {
            Gc[[1L]] <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
            gc[[1L]] <- Gc[[1L]] * mu + L2 * go[[1L]][[1L]]
            dks[1L, k + 1L, 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                                L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] * mu +
                                   L1 * gc[[i]] + L2 * go[[1L]][[i + 1L]]
                dks[i + 1L, k + 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
            }
            Gn <- list(Gc)
            gn <- list(gc)
        } else {
            Gn <- list(0)
            gn <- list(0)
        }
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                if(k - j <= q && j <= r) {
                    Gc[[1L]] <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                    gc[[1L]] <- Gc[[1L]] * mu + L2 * go[[j + 1L]][[1L]] + L3 * go[[j]][[1L]]
                    dks[1L, k - j + 1L, j + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
                    for(i in 1L:p) {
                        Gc[[i + 1L]] <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                                        L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                                        L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                        gc[[i + 1L]] <- Gc[[i + 1L]] * mu +
                                        L1 * gc[[i]] + L2 * go[[j + 1L]][[i + 1L]] + L3 * go[[j]][[i + 1L]]
                        dks[i + 1L, k - j + 1L, j + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
                    }
                    Gn <- c(Gn, list(Gc))
                    gn <- c(gn, list(gc))
                } else {
                    Gn <- c(Gn, list(0))
                    gn <- c(gn, list(0))
                }
            }
        }
        if(k <= r) {
            Gc[[1L]] <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
            gc[[1L]] <- Gc[[1L]] * mu + L3 * go[[k]][[1L]]
            dks[1L, 1L, k + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                                L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] * mu +
                                L1 * gc[[i]] + L3 * go[[k]][[i + 1L]]
                dks[i + 1L, 1L, k + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
            }
            Gn <- c(Gn, list(Gc))
            gn <- c(gn, list(gc))
        } else {
            Gn <- c(Gn, list(0))
            gn <- c(gn, list(0))
        }
    }
    return(dks)
}


##### arl #####
#' Recursion for a_{r,l}
#'
#' \code{arl()} is an internal function to calculate \eqn{a_{r,l}} as defined
#' in Hillier et al. (2014; eq. 24), which is used in the calculation of
#' the moment of such a ratio of quadratic forms in normal variables
#' where the denominator matrix is identity.
#'
#' This function implements the super-short recursion described in
#' Hillier et al. (2014  eqs. 31--32).
#' Note that \eqn{w_{r,i}} there should be understood as \eqn{w_{r,l,i}} with
#' the index \eqn{l} fixed for each \eqn{a_{r,l}}.
#'
#' The \code{matrix} method just calculates \code{L} and \code{D} from
#' \code{A} and \code{mu} and passes them to the \code{default} method.
#'
#' @param L
#'   Eigenvalues of the argument matrix; vector of \eqn{\lambda_i}
#' @param A
#'   Argument matrix.  Assumed to be symmetric.
#' @param mu
#'   Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}
#' @param D
#'   Squared norm of the mean vector projected on the eigenvalues of
#'   the argument matrix: vectors of \eqn{\delta_i}
#' @param m
#'   Scalar to specify the desired order
#' @param ...
#'   Additional arguments passed to internal methods
#'
#' @seealso
#' \code{\link{qfrm}}; this function is used in \code{qfrm_ApIq_int()}
#' (for noncentral cases only)
#'
#' @name arl
#' @order 1
#'
arl <- function(L, ...) {
    UseMethod("arl", L)
}

#' @rdname arl
#' @order 3
#'
#' @exportS3Method
#'
arl.matrix <- function(A, mu, m = 10L) {
    if(any(dim(A) == 1L)) return(arl.default(A, D, m))
    eigA <- eigen(A, symmetric = TRUE)
    L <- eigA$values
    D <- c(crossprod(eigA$vectors, mu)) ^ 2
    return(arl.default(L, D, m))
}

#' @rdname arl
#' @order 2
#'
#' @exportS3Method
#'
arl.default <-  function(L, D, m = 10L) {
    n <- length(L)
    arls <- matrix(0, m + 1L, m + 1L)
    arls[, 1L] <- d1_i(L, m)
    wrls <- matrix(0, n, m)
    for(k in 1:m) {
        for(l in 1:k) {
            wrls[, l] <- L * (arls[k, l] + wrls[, l])
            arls[k + 1L, l + 1L] <- sum(D * wrls[, l])
        }
    }
    return(arls)
}


# d2_ij_mb <- function(A1, A2, m = 100L, p = m, q = m, fill_all = !missing(p) || !missing(q)) {
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     Gs <- array(0, dim = c(n, n, dim(dks)))
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeromat
#             if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[, , i1, i2 + 1L])
#             if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[, , i1 + 1L, i2])
#             Gs[, , i1 + 1L, i2 + 1L] <- tG
#             dks[i1 + 1L, i2 + 1L] <- sum(tG) / (2 * k)
#         }
#     }
#     return(dks)
# }

# d2_ij_vb <- function(L1, L2, m = 100L, p = m, q = m, fill_all = !missing(p) || !missing(q)) {
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- array(0, dim = c(n, dim(dks)))
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeros
#             if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[, i1, i2 + 1L])
#             if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[, i1 + 1L, i2])
#             Gs[, i1 + 1L, i2 + 1L] <- tG
#             dks[i1 + 1L, i2 + 1L] <- sum(tG) / (2 * k)
#         }
#         if(max(unlist(Gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- Gs / 1e10
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#     }
#     return(dks)
# }

##### d2_ij_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
d2_ij_m <- function(A1, A2, m = 100L, p = m, q = m, fill_all = !missing(p) || !missing(q)) {
    il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    order_mat <- outer(0:p, 0:q, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            dks[i1 + 1L, i2 + 1L] <- tr(tG) / (2 * k)
        }
        if(max(unlist(Gs)) > thr) {
            ind_dks <- which(order_mat == k)
            ind_lscf <- which(order_mat >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d2_ij_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
d2_ij_v <- function(L1, L2, m = 100L, p = m, q = m, fill_all = !missing(p) || !missing(q)) {
    il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    order_mat <- outer(0:p, 0:q, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            dks[i1 + 1L, i2 + 1L] <- sum(tG) / (2 * k)
        }
        if(max(unlist(Gs)) > thr) {
            ind_dks <- which(order_mat == k)
            ind_lscf <- which(order_mat >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d2_pj_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
d2_pj_m <- function(A1, A2, m = 100L, p = 1L) {
    if(p == 1L) return(d2_1j_m(A1, A2, m))
    n <- ncol(A1)
    p1 <- p + 1L
    m1 <- m + 1L
    In <- diag(n)
    dks <- matrix(0, p1, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    G_k_i <- list()
    G_k_i[1L:p1] <- list(zeromat)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        dks[i + 1L, 1L] <- tr(G_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        G_k_i[[1L]] <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        dks[1L, k + 1L] <- tr(G_k_i[[1L]]) / (2 * k)
        for(i in 1L:p) {
            G_k_i[[i + 1L]] <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            dks[i + 1L, k + 1L] <- tr(G_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(G_k_i[[i + 1L]]) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d2_1j_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
d2_1j_m <- function(A1, A2, m = 100L) {
    n <- ncol(A1)
    m1 <- m + 1L
    In <- diag(n)
    dks <- matrix(0, 2L, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1
    dks[2L, 1L] <- tr(G_k_1) / 2
    for(k in 1L:m) {
        G_k_0 <- A2 %*% (dks[1L, k] * In + G_k_0)
        dks[1L, k + 1L] <- tr(G_k_0) / (2 * k)
        G_k_1 <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        dks[2L, k + 1L] <- tr(G_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d2_pj_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
d2_pj_v <- function(L1, L2, m = 100L, p = 1L) {
    if(p == 1L) return(d2_1j_v(L1, L2, m))
    n <- length(L1)
    p1 <- p + 1L
    m1 <- m + 1L
    dks <- matrix(0, p1, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    G_k_i[1L:p1] <- list(zeros)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        G_k_i[[1L]] <- L2 * (dks[1L, k] + G_k_i[[1L]])
        dks[1L, k + 1L] <- sum(G_k_i[[1L]]) / (2 * k)
        for(i in 1L:p) {
            G_k_i[[i + 1L]] <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(G_k_i[[i + 1L]]) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d2_1j_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
d2_1j_v <- function(L1, L2, m = 100L) {
    n <- length(L1)
    m1 <- m + 1L
    dks <- matrix(0, 2L, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1
    dks[2L, 1L] <- sum(G_k_1) / 2
    for(k in 1L:m) {
        G_k_0 <- L2 * (dks[1L, k] + G_k_0)
        dks[1L, k + 1L] <- sum(G_k_0) / (2 * k)
        G_k_1 <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        dks[2L, k + 1L] <- sum(G_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}


# d2_i1_m <- function(A1, A2, m = 101L) {
#     # tr <- function(X) sum(diag(X))
#     n <- ncol(A1)
#     In <- diag(n)
#     dks <- matrix(0, m + 1L, 2L)
#     dks[1L, 1L] <- 1
#     dks[1L, 2L] <- tr(A2) / 2
#     dks[2L, 1L] <- tr(A1) / 2
#     G_k2_1 <- A2
#     G_k1_0 <- A1
#     for(k in 2L:m) {
#         G_k1_1 <- A1 %*% (dks[k - 1L, 2L] * In + G_k2_1) +
#                   A2 %*% (dks[k, 1L] * In + G_k1_0)
#         dks[k, 2L] <- tr(G_k1_1) / (2 * k)
#         G_k_0 <- A1 %*% (dks[k, 1L] * In + G_k1_0)
#         dks[k + 1L, 1L] <- tr(G_k_0) / (2 * k)
#         G_k2_1 <- G_k1_1
#         G_k1_0 <- G_k_0
#     }
#     return(dks)
# }
#
# d2_i1_v <- function(L1, L2, m = 101L) {
#     n <- length(L1)
#     dks <- matrix(0, m + 1L, 2L)
#     dks[1L, 1L] <- 1
#     dks[1L, 2L] <- sum(L2) / 2
#     dks[2L, 1L] <- sum(L1) / 2
#     G_k2_1 <- L2
#     G_k1_0 <- L1
#     for(k in 2L:m) {
#         G_k1_1 <- L1 * (dks[k - 1L, 2L] + G_k2_1) +
#                   L2 * (dks[k, 1L] + G_k1_0)
#         dks[k, 2L] <- sum(G_k1_1) / (2 * k)
#         G_k_0 <- L1 * (dks[k, 1L] + G_k1_0)
#         dks[k + 1L, 1L] <- sum(G_k_0) / (2 * k)
#         G_k2_1 <- G_k1_1
#         G_k1_0 <- G_k_0
#     }
#     return(dks)
# }

# d2_ip_v <- function(L1, L2, m = 100L, p = 1L) {
#     n <- length(L1)
#     dks <- matrix(0, m + 1L, p + 1L)
#     dks[1L, 1L] <- 1
#     zeros <- rep.int(0, n)
#     G_k_i <- matrix(0, n, p + 1L)
#     for(i in 1L:p) {
#         G_k_i[, i + 1L] <- L2 * (dks[1L, i] + G_k_i[, i])
#         dks[1L, i + 1L] <- sum(G_k_i[, i + 1L]) / (2 * i)
#     }
#     for(k in 1L:m) {
#         G_k_i[, 1L] <- L1 * (dks[k, 1L] + G_k_i[, 1L])
#         dks[k + 1L, 1L] <- sum(G_k_i[, 1L]) / (2 * k)
#         for(i in 1L:p) {
#             G_k_i[, i + 1L] <- L1 * (dks[k, i + 1L] + G_k_i[, i + 1L]) +
#                                L2 * (dks[k + 1L, i] + G_k_i[, i])
#             dks[k + 1L, i + 1L] <- sum(G_k_i[, i + 1L]) / (2 * (k + i))
#         }
#     }
#     return(dks)
# }
#
# d2_ip_v <- function(L1, L2, m = 100L, p = 1L) {
#     n <- length(L1)
#     dks <- matrix(0, m + 1L, p + 1L)
#     dks[1L, 1L] <- 1
#     zeros <- rep.int(0, n)
#     G_k_i <- list()
#     G_k_i[1L:(p + 1L)] <- list(zeros)
#     for(i in 1L:p) {
#         G_k_i[[i + 1L]] <- L2 * (dks[1L, i] + G_k_i[[i]])
#         dks[1L, i + 1L] <- sum(G_k_i[[i + 1L]]) / (2 * i)
#     }
#     for(k in 1L:m) {
#         G_k_i[[1L]] <- L1 * (dks[k, 1L] + G_k_i[[1L]])
#         dks[k + 1L, 1L] <- sum(G_k_i[[1L]]) / (2 * k)
#         for(i in 1L:p) {
#             G_k_i[[i + 1L]] <- L1 * (dks[k, i + 1L] + G_k_i[[i + 1L]]) +
#                                L2 * (dks[k + 1L, i] + G_k_i[[i]])
#             dks[k + 1L, i + 1L] <- sum(G_k_i[[i + 1L]]) / (2 * (k + i))
#         }
#     }
#     return(dks)
# }


##### d3_ijk_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
d3_ijk_m <- function(A1, A2, A3, m = 100L, p = m, q = m, r = m,
                 fill_across = c(!missing(p), !missing(q), !missing(r))) {
    il3 <- function(i1, i2, i3) i1 + i2 * p1 + i3 * p1 * q1 + 1L
    n <- ncol(A1)
    p1 <- p + 1L
    q1 <- q + 1L
    r1 <- r + 1L
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    dks <- array(0, dim = c(p1, q1, r1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
            for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + r > kmax) next
                if(fill_across[2L] && i1 + i3 + q > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- tr(tG) / (2 * k)
            }
        }
        if(max(unlist(Gs)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d3_ijk_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
d3_ijk_v <- function(L1, L2, L3, m = 100L, p = m, q = m, r = m,
                 fill_across = c(!missing(p), !missing(q), !missing(r))) {
    il3 <- function(i1, i2, i3) i1 + i2 * p1 + i3 * p1 * q1 + 1L
    n <- length(L1)
    p1 <- p + 1L
    q1 <- q + 1L
    r1 <- r + 1L
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(p1, q1, r1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
            for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + r > kmax) next
                if(fill_across[2L] && i1 + i3 + q > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- sum(tG) / (2 * k)
            }
        }
        if(max(unlist(Gs)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d3_pjk_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
d3_pjk_m <- function(A1, A2, A3, m = 100L, p = 1L) {
    n <- ncol(A1)
    p1 <- p + 1L
    m1 <- m + 1L
    In <- diag(n)
    dks <- array(0, dim = c(p1, m1, m1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    order_array <- outer(outer(rep.int(0, p1), 0:m, "+"), 0:m, "+")
    zeromat <- matrix(0, n, n)
    Gc <- list()
    Gc[1L:p1] <- list(zeromat)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        dks[i + 1L, 1L, 1L] <- tr(Gc[[i + 1L]]) / (2 * i)
    }
    Gn <- list(Gc)
    for(k in 1L:m) {
        Go <- Gn
        Gc[[1L]] <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
        dks[1L, k + 1L, 1L] <- tr(Gc[[1L]]) / (2 * k)
        for(i in 1L:p) {
            Gc[[i + 1L]] <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                            A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
            dks[i + 1L, k + 1L, 1L] <- tr(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- list(Gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                Gc[[1L]] <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                dks[1L, k - j + 1L, j + 1L] <- tr(Gc[[1L]]) / (2 * k)
                for(i in 1L:p) {
                    Gc[[i + 1L]] <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                                    A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                                    A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                    dks[i + 1L, k - j + 1L, j + 1L] <- tr(Gc[[i + 1L]]) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
            }
        }
        Gc[[1L]] <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
        dks[1L, 1L, k + 1L] <- tr(Gc[[1L]]) / (2 * k)
        for(i in 1L:p) {
            Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                            A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
            dks[i + 1L, 1L, k + 1L] <- tr(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        if(max(unlist(Gn)) > thr) {
            ind_dks <- which(order_array == k + 1L)
            ind_lscf <- which(order_array >= k + 1L)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### d3_pjk_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
d3_pjk_v <- function(L1, L2, L3, m = 100L, p = 1L) {
    n <- length(L1)
    p1 <- p + 1L
    m1 <- m + 1L
    dks <- array(0, dim = c(p1, m1, m1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    order_array <- outer(outer(rep.int(0, p1), 0:m, "+"), 0:m, "+")
    zeros <- rep.int(0, n)
    Gc <- list()
    Gc[1L:p1] <- list(zeros)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        dks[i + 1L, 1L, 1L] <- sum(Gc[[i + 1L]]) / (2 * i)
    }
    Gn <- list(Gc)
    for(k in 1L:m) {
        Go <- Gn
        Gc[[1L]] <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
        dks[1L, k + 1L, 1L] <- sum(Gc[[1L]]) / (2 * k)
        for(i in 1L:p) {
            Gc[[i + 1L]] <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                            L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
            dks[i + 1L, k + 1L, 1L] <- sum(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- list(Gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                Gc[[1L]] <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                dks[1L, k - j + 1L, j + 1L] <- sum(Gc[[1L]]) / (2 * k)
                for(i in 1L:p) {
                    Gc[[i + 1L]] <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                                    L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                                    L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                    dks[i + 1L, k - j + 1L, j + 1L] <- sum(Gc[[i + 1L]]) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
            }
        }
        Gc[[1L]] <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
        dks[1L, 1L, k + 1L] <- sum(Gc[[1L]]) / (2 * k)
        for(i in 1L:p) {
            Gc[[i + 1L]] <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                            L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
            dks[i + 1L, 1L, k + 1L] <- sum(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        if(max(unlist(Gn)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}


# d3_ij1_m <- function(A1, A2, A3, m = 101L) {
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     dks <- array(0, dim = c(m + 1L, m + 1L, 2L))
#     dks[1L] <- 1
#     Gs <- array(0, dim = c(n, n, dim(dks)))
#     for(k in 1L:m) {
#         for(i3 in 0L:1L) {
#             for(i1 in (k - i3):0L) {
#                 i2 <- k - i1 - i3
#                 tG <- zeromat
#                 if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[, , i1, i2 + 1L, i3 + 1L])
#                 if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[, , i1 + 1L, i2, i3 + 1L])
#                 if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[, , i1 + 1L, i2 + 1L, i3])
#                 Gs[, , i1 + 1L, i2 + 1L, i3 + 1L] <- tG
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- tr(tG) / (2 * k)
#             }
#         }
#     }
#     return(dks)
# }
#
# d3_ij1_v <- function(A1, A2, A3, m = 101L) {
#     n <- ncol(A1)
#     zeros <- rep.int(0, n)
#     dks <- array(0, dim = c(m + 1L, m + 1L, 2L))
#     dks[1L] <- 1
#     Gs <- array(0, dim = c(n, n, dim(dks)))
#     for(k in 1L:m) {
#         for(i3 in 0L:1L) {
#             for(i1 in (k - i3):0L) {
#                 i2 <- k - i1 - i3
#                 tG <- zeros
#                 if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[, , i1, i2 + 1L, i3 + 1L])
#                 if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[, , i1 + 1L, i2, i3 + 1L])
#                 if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[, , i1 + 1L, i2 + 1L, i3])
#                 Gs[, , i1 + 1L, i2 + 1L, i3 + 1L] <- tG
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- sum(tG) / (2 * k)
#             }
#         }
#     }
#     return(dks)
# }

# h1_i_m <- function(A1, mu = rep.int(0, n), m = 100L) {
#     n <- ncol(A1)
#     In <- diag(n)
#     dks <- dks <- rep.int(c(1, 0), c(1L, m))
#     Go <- matrix(0, n, n)
#     go <- matrix(0, n, 1)
#     for(k in 1L:m) {
#         tG <- A1 %*% (dks[k] * In + Go)
#         tg <- (tG - Go) %*% mu - dks[k] * mu + A1 %*% go
#         dks[k + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
#         Go <- tG
#         go <- tg
#     }
#     return(dks)
# }

##### h2_1j_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
h2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, p = m, q = m,
                      fill_all = !missing(p) || !missing(q)) {
    il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_mat <- outer(0:p, 0:q, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG %*% mu
            if(i1 >= 1L) tg <- tg - Gs[[il2(i1 - 1L, i2)]] %*% mu - dks[i1, i2 + 1L] * mu + A1 %*% gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg - Gs[[il2(i1, i2 - 1L)]] %*% mu - dks[i1 + 1L, i2] * mu + A2 %*% gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
        }
        if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
            ind_dks <- which(order_mat == k)
            ind_lscf <- which(order_mat >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### h2_1j_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
h2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, p = m, q = m,
                      fill_all = !missing(p) || !missing(q)) {
    il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_mat <- outer(0:p, 0:q, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG * mu
            if(i1 >= 1L) tg <- tg - (Gs[[il2(i1 - 1L, i2)]] + dks[i1, i2 + 1L]) * mu + L1 * gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg - (Gs[[il2(i1, i2 - 1L)]] + dks[i1 + 1L, i2]) * mu + L2 * gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
        }
        if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
            ind_dks <- which(order_mat == k)
            ind_lscf <- which(order_mat >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

# htil2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, p = m, q = m,
#                          fill_all = !missing(p) || !missing(q)) {
#     il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     zerovec <- matrix(0, n, 1)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeromat)
#     gs <- list(zerovec)
#     order_mat <- outer(0:p, 0:q, "+")
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         Gs[which(order_mat == k - 2L)] <- 0
#         gs[which(order_mat == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeromat
#             if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
#             if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
#             Gs[[il2(i1, i2)]] <- tG
#             tg <- tG %*% mu
#             if(i1 >= 1L) tg <- tg + A1 %*% gs[[il2(i1 - 1L, i2)]]
#             if(i2 >= 1L) tg <- tg - Gs[[il2(i1, i2 - 1L)]] %*% mu - dks[i1 + 1L, i2] * mu + A2 %*% gs[[il2(i1, i2 - 1L)]]
#             gs[[il2(i1, i2)]] <- tg
#             dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#     }
#     return(dks)
# }

# htil2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, p = m, q = m,
#                          fill_all = !missing(p) || !missing(q)) {
#     il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeros)
#     gs <- list(zeros)
#     order_mat <- outer(0:p, 0:q, "+")
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         Gs[which(order_mat == k - 2L)] <- 0
#         gs[which(order_mat == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeros
#             if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
#             if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
#             Gs[[il2(i1, i2)]] <- tG
#             tg <- tG * mu
#             if(i1 >= 1L) tg <- tg + L1 * gs[[il2(i1 - 1L, i2)]]
#             if(i2 >= 1L) tg <- tg - (Gs[[il2(i1, i2 - 1L)]] + dks[i1 + 1L, i2]) * mu + L2 * gs[[il2(i1, i2 - 1L)]]
#             gs[[il2(i1, i2)]] <- tg
#             dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#     }
#     return(dks)
# }

##### htil2_pj_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
htil2_pj_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, p = 1L) {
    if(p == 1L) return(htil2_1j_m(A1, A2, mu, m))
    n <- ncol(A1)
    p1 <- p + 1L
    m1 <- m + 1L
    In <- diag(n)
    dks <- matrix(0, p1, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:p1] <- list(zeromat)
    g_k_i[1L:p1] <- list(zerovec)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu + A1 %*% g_k_i[[i]]
        dks[i + 1L, 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG - G_k_i[[1L]] - (dks[1L, k] * In)) %*% mu + A2 %*% g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- (tr(G_k_i[[1L]]) + c(crossprod(mu, g_k_i[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            tG <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG - G_k_i[[i + 1L]]
                                - (dks[i + 1L, k] * In)) %*% mu +
                               A1 %*% g_k_i[[i]] + A2 %*% g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### htil2_1j_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
htil2_1j_m <- function(A1, A2, mu = rep.int(0, n), m = 100L) {
    n <- ncol(A1)
    m1 <- m + 1L
    In <- diag(n)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1
    g_k_0 <- matrix(0, n, 1)
    g_k_1 <- G_k_1 %*% mu
    dks[2L, 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / 2
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_0)
        g_k_0 <- (tG - G_k_0 - (dks[1L, k] * In)) %*% mu + A2 %*% g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- (tr(G_k_0) + c(crossprod(mu, g_k_0))) / (2 * k)
        tG <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        g_k_1 <- (tG - G_k_1 - (dks[2L, k] * In)) %*% mu +
                 A1 %*% g_k_0 + A2 %*% g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### htil2_pj_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
htil2_pj_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, p = 1L) {
    if(p == 1L) return(htil2_1j_v(L1, L2, mu, m))
    n <- length(L1)
    p1 <- p + 1L
    m1 <- m + 1L
    dks <- matrix(0, p1, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:p1] <- list(zeros)
    g_k_i[1L:p1] <- list(zeros)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu + L1 * g_k_i[[i]]
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG - G_k_i[[1L]] - dks[1L, k]) * mu + L2 * g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- sum(G_k_i[[1L]] + mu * g_k_i[[1L]]) / (2 * k)
        for(i in 1L:p) {
            tG <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG - G_k_i[[i + 1L]] - dks[i + 1L, k]) * mu +
                               L1 * g_k_i[[i]] + L2 * g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### htil2_1j_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
htil2_1j_v <- function(L1, L2, mu = rep.int(0, n), m = 100L) {
    n <- length(L1)
    m1 <- m + 1L
    dks <- matrix(0, 2L, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1
    g_k_0 <- rep.int(0, n)
    g_k_1 <- G_k_1 * mu
    dks[2L, 1L] <- sum(G_k_1 + mu * g_k_1) / 2
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_0)
        g_k_0 <- (tG - G_k_0 - dks[1L, k]) * mu + L2 * g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- sum(G_k_0 + mu * g_k_0) / (2 * k)
        tG <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        g_k_1 <- (tG - G_k_1 - dks[2L, k]) * mu + L1 * g_k_0 + L2 * g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- sum(G_k_1 + mu * g_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}


##### h3_ijk_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
h3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
                 fill_across = c(!missing(p), !missing(q), !missing(r))) {
    il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
            for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + r > kmax) next
                if(fill_across[2L] && i1 + i3 + q > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG %*% mu
                if(i1 >= 1L) tg <- tg - Gs[[il3(i1 - 1L, i2, i3)]] %*% mu - dks[i1, i2 + 1L, i3 + 1L] * mu + A1 %*% gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg - Gs[[il3(i1, i2 - 1L, i3)]] %*% mu - dks[i1 + 1L, i2, i3 + 1L] * mu + A2 %*% gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg - Gs[[il3(i1, i2, i3 - 1L)]] %*% mu - dks[i1 + 1L, i2 + 1L, i3] * mu + A3 %*% gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
            }
        }
        if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### h3_ijk_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
h3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
                 fill_across = c(!missing(p), !missing(q), !missing(r))) {
    il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
            for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + r > kmax) next
                if(fill_across[2L] && i1 + i3 + q > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG * mu
                if(i1 >= 1L) tg <- tg - (Gs[[il3(i1 - 1L, i2, i3)]] + dks[i1, i2 + 1L, i3 + 1L]) * mu + L1 * gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg - (Gs[[il3(i1, i2 - 1L, i3)]] + dks[i1 + 1L, i2, i3 + 1L]) * mu + L2 * gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg - (Gs[[il3(i1, i2, i3 - 1L)]] + dks[i1 + 1L, i2 + 1L, i3]) * mu + L3 * gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
            }
        }
        if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

# htil3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
#                  fill_across = c(!missing(p), !missing(q), !missing(r))) { # , verbose = m > 200L) {
#     il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     zerovec <- matrix(0, n, 1)
#     dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
#     dks[1L] <- 1
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeromat)
#     gs <- list(zerovec)
#     order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
#     kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
#     for(k in 1L:kmax) {
#         Gs[which(order_array == k - 2L)] <- 0
#         gs[which(order_array == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
#             for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
#                 i3 <- k - i1 - i2
#                 if(fill_across[3L] && i1 + i2 + r > kmax) next
#                 if(fill_across[2L] && i1 + i3 + q > kmax) next
#                 # if(fill_across[1] && i2 + i3 + p > kmax) next
#                 tG <- zeromat
#                 if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
#                 if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
#                 if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
#                 Gs[[il3(i1, i2, i3)]] <- tG
#                 tg <- tG %*% mu
#                 if(i1 >= 1L) tg <- tg + A1 %*% gs[[il3(i1 - 1L, i2, i3)]]
#                 if(i2 >= 1L) tg <- tg - Gs[[il3(i1, i2 - 1L, i3)]] %*% mu - dks[i1 + 1L, i2, i3 + 1L] * mu + A2 %*% gs[[il3(i1, i2 - 1L, i3)]]
#                 if(i3 >= 1L) tg <- tg - Gs[[il3(i1, i2, i3 - 1L)]] %*% mu - dks[i1 + 1L, i2 + 1L, i3] * mu + A3 %*% gs[[il3(i1, i2, i3 - 1L)]]
#                 gs[[il3(i1, i2, i3)]] <- tg
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
#             }
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#         # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
#     }
#     return(dks)
# }

# htil3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
#                  fill_across = c(!missing(p), !missing(q), !missing(r))) { # , verbose = m > 200L) {
#     il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
#     dks[1L] <- 1
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeros)
#     gs <- list(zeros)
#     order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
#     kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
#     for(k in 1L:kmax) {
#         Gs[which(order_array == k - 2L)] <- 0
#         gs[which(order_array == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
#             for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
#                 i3 <- k - i1 - i2
#                 if(fill_across[3L] && i1 + i2 + r > kmax) next
#                 if(fill_across[2L] && i1 + i3 + q > kmax) next
#                 # if(fill_across[1] && i2 + i3 + p > kmax) next
#                 tG <- zeros
#                 if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
#                 if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
#                 if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
#                 Gs[[il3(i1, i2, i3)]] <- tG
#                 tg <- tG * mu
#                 if(i1 >= 1L) tg <- tg + L1 * gs[[il3(i1 - 1L, i2, i3)]]
#                 if(i2 >= 1L) tg <- tg - (Gs[[il3(i1, i2 - 1L, i3)]] + dks[i1 + 1L, i2, i3 + 1L]) * mu + L2 * gs[[il3(i1, i2 - 1L, i3)]]
#                 if(i3 >= 1L) tg <- tg - (Gs[[il3(i1, i2, i3 - 1L)]] + dks[i1 + 1L, i2 + 1L, i3]) * mu + L3 * gs[[il3(i1, i2, i3 - 1L)]]
#                 gs[[il3(i1, i2, i3)]] <- tg
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
#             }
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#         # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
#     }
#     return(dks)
# }

##### htil3_pjk_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
htil3_pjk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, p = 1L) {
    n <- ncol(A1)
    p1 <- p + 1L
    m1 <- m + 1L
    In <- diag(n)
    dks <- array(0, dim = c(p1, m1, m1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    order_array <- outer(outer(rep.int(0, p1), 0:m, "+"), 0:m, "+")
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    Gc <- list()
    gc <- list()
    Gc[1L:p1] <- list(zeromat)
    gc[1L:p1] <- list(zerovec)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu + A1 %*% gc[[i]]
        dks[i + 1L, 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
        gc[[1L]] <- (tG - Go[[1L]][[1L]] - (dks[1L, k, 1L] * In)) %*% mu + A2 %*% go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            tG <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                  A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG - Go[[1L]][[i + 1L]]
                                - (dks[i + 1L, k, 1L] * In)) %*% mu +
                               A1 %*% gc[[i]] + A2 %*% go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                gc[[1L]] <- (tG - Go[[j + 1L]][[1L]] - Go[[j]][[1L]] - ((dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j]) * In)) %*% mu + A2 %*% go[[j + 1L]][[1L]] + A3 %*% go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
                for(i in 1L:p) {
                    tG <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                          A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                          A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG - Go[[j + 1L]][[i + 1L]] - Go[[j]][[i + 1L]]
                                     - ((dks[i + 1L, k - j, j + 1L] + dks[i + 1L, k - j + 1L, j]) * In)) %*% mu +
                                    A1 %*% gc[[i]] + A2 %*% go[[j + 1L]][[i + 1L]] + A3 %*% go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
        gc[[1L]] <- (tG - Go[[k]][[1L]] - (dks[1L, 1L, k] * In)) %*% mu + A3 %*% go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            tG <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                  A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG - Go[[k]][[i + 1L]]
                             - (dks[i + 1L, 1L, k] * In)) %*% mu +
                            A1 %*% gc[[i]] + A3 %*% go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gn)) > thr || max(unlist(gn)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### htil3_pjk_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
htil3_pjk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, p = 1L) {
    n <- length(L1)
    p1 <- p + 1L
    m1 <- m + 1L
    dks <- array(0, dim = c(p1, m1, m1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    order_array <- outer(outer(rep.int(0, p1), 0:m, "+"), 0:m, "+")
    zeros <- rep.int(0, n)
    Gc <- list()
    gc <- list()
    Gc[1L:p1] <- list(zeros)
    gc[1L:p1] <- list(zeros)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] * mu + L1 * gc[[i]]
        dks[i + 1L, 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
        gc[[1L]] <- (tG - (Go[[1L]][[1L]] + dks[1L, k, 1L])) * mu + L2 * go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:p) {
            tG <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                  L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG - (Go[[1L]][[i + 1L]]
                                + dks[i + 1L, k, 1L])) * mu +
                               L1 * gc[[i]] + L2 * go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                gc[[1L]] <- (tG - (Go[[j + 1L]][[1L]] + Go[[j]][[1L]] + dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j])) * mu + L2 * go[[j + 1L]][[1L]] + L3 * go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
                for(i in 1L:p) {
                    tG <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                          L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                          L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG - (Go[[j + 1L]][[i + 1L]] + Go[[j]][[i + 1L]]
                                        + dks[i + 1L, k - j, j + 1L]
                                        + dks[i + 1L, k - j + 1L, j])) * mu +
                                       L1 * gc[[i]] + L2 * go[[j + 1L]][[i + 1L]] + L3 * go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
        gc[[1L]] <- (tG - (Go[[k]][[1L]] + dks[1L, 1L, k])) * mu + L3 * go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:p) {
            tG <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                  L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG - (Go[[k]][[i + 1L]]
                             + dks[i + 1L, 1L, k])) * mu +
                            L1 * gc[[i]] + L3 * go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gn)) > thr || max(unlist(gn)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}



# #' Recursion for hhat
# #'
# #' \code{hhat2_ij_m()} and \code{hhat2_ij_m()} are recursive algorithms
# #' for truncation error in noncentral case.
# #'
# #' These provide recursive algorithms for \eqn{\hat{h}_{i,j}} in
# #' Hillier et al. (2014, theorem 7).  This recursion is said to be very similar
# #' to those for \eqn{h_{i,j}} in the note therein, but differs in the signs
# #' of some terms.
# #'
# hhat2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, p = m, q = m, fill_all = !missing(p) || !missing(q)) {
#     il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     zerovec <- matrix(0, n, 1)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeromat)
#     gs <- list(zerovec)
#     order_mat <- outer(0:p, 0:q, "+")
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         Gs[which(order_mat == k - 2L)] <- 0
#         gs[which(order_mat == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeromat
#             if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
#             if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
#             Gs[[il2(i1, i2)]] <- tG
#             tg <- tG %*% mu
#             if(i1 >= 1L) tg <- tg + A1 %*% gs[[il2(i1 - 1L, i2)]]
#             if(i2 >= 1L) tg <- tg + Gs[[il2(i1, i2 - 1L)]] %*% mu + dks[i1 + 1L, i2] * mu + A2 %*% gs[[il2(i1, i2 - 1L)]]
#             gs[[il2(i1, i2)]] <- tg
#             dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#     }
#     return(dks)
# }

# hhat2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, p = m, q = m, fill_all = !missing(p) || !missing(q)) {
#     il2 <- function(i1, i2) i1 + i2 * (p + 1L) + 1L
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (p + 1) * (q + 1) - 1L)), p + 1, q + 1)
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeros)
#     gs <- list(zeros)
#     order_mat <- outer(0:p, 0:q, "+")
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         Gs[which(order_mat == k - 2L)] <- 0
#         gs[which(order_mat == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q, 0L)) {
#             i2 <- k - i1
#             tG <- zeros
#             if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
#             if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
#             Gs[[il2(i1, i2)]] <- tG
#             tg <- tG * mu
#             if(i1 >= 1L) tg <- tg + L1 * gs[[il2(i1 - 1L, i2)]]
#             if(i2 >= 1L) tg <- tg + (Gs[[il2(i1, i2 - 1L)]] + dks[i1 + 1L, i2]) * mu + L2 * gs[[il2(i1, i2 - 1L)]]
#             gs[[sum(c(i1, i2) * c(1L, p + 1L)) + 1L]] <- tg
#             dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#     }
#     return(dks)
# }

##### hhat2_pj_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
hhat2_pj_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, p = 1L) {
    if(p == 1L) return(hhat2_1j_m(A1, A2, mu, m))
    n <- ncol(A1)
    p1 <- p + 1L
    m1 <- m + 1L
    In <- diag(n)
    dks <- matrix(0, p1, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:p1] <- list(zeromat)
    g_k_i[1L:p1] <- list(zerovec)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu + A1 %*% g_k_i[[i]]
        dks[i + 1L, 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG + G_k_i[[1L]] + (dks[1L, k] * In)) %*% mu + A2 %*% g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- (tr(G_k_i[[1L]]) + c(crossprod(mu, g_k_i[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            tG <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG + G_k_i[[i + 1L]]
                                + (dks[i + 1L, k] * In)) %*% mu +
                               A1 %*% g_k_i[[i]] + A2 %*% g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### hhat2_1j_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
hhat2_1j_m <- function(A1, A2, mu = rep.int(0, n), m = 100L) {
    n <- ncol(A1)
    m1 <- m + 1L
    In <- diag(n)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1
    g_k_0 <- matrix(0, n, 1)
    g_k_1 <- G_k_1 %*% mu
    dks[2L, 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / 2
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_0)
        g_k_0 <- (tG + G_k_0 + (dks[1L, k] * In)) %*% mu + A2 %*% g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- (tr(G_k_0) + c(crossprod(mu, g_k_0))) / (2 * k)
        tG <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        g_k_1 <- (tG + G_k_1 + (dks[2L, k] * In)) %*% mu +
                 A1 %*% g_k_0 + A2 %*% g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### hhat2_pj_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
hhat2_pj_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, p = 1L) {
    if(p == 1L) return(hhat2_1j_v(L1, L2, mu, m))
    n <- length(L1)
    p1 <- p + 1L
    m1 <- m + 1L
    dks <- matrix(0, p1, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:p1] <- list(zeros)
    g_k_i[1L:p1] <- list(zeros)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu + L1 * g_k_i[[i]]
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG + G_k_i[[1L]] + dks[1L, k]) * mu + L2 * g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- sum(G_k_i[[1L]] + mu * g_k_i[[1L]]) / (2 * k)
        for(i in 1L:p) {
            tG <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG + G_k_i[[i + 1L]] + dks[i + 1L, k]) * mu +
                               L1 * g_k_i[[i]] + L2 * g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### hhat2_1j_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with two matrices
#'
#' @rdname d2_ij
#'
hhat2_1j_v <- function(L1, L2, mu = rep.int(0, n), m = 100L) {
    n <- length(L1)
    m1 <- m + 1L
    dks <- matrix(0, 2L, m1)
    dks[1L, 1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1
    g_k_0 <- rep.int(0, n)
    g_k_1 <- G_k_1 * mu
    dks[2L, 1L] <- sum(G_k_1 + mu * g_k_1) / 2
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_0)
        g_k_0 <- (tG + G_k_0 + dks[1L, k]) * mu + L2 * g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- sum(G_k_0 + mu * g_k_0) / (2 * k)
        tG <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        g_k_1 <- (tG + G_k_1 + dks[2L, k]) * mu + L1 * g_k_0 + L2 * g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- sum(G_k_1 + mu * g_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks[, k + 1L] <- dks[, k + 1L] / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            lscf[, (k + 1L):m1] <- lscf[, (k + 1L):m1] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}



# hhat3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
#                  fill_across = c(!missing(p), !missing(q), !missing(r))) { # , verbose = m > 200L) {
#     il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     zerovec <- matrix(0, n, 1)
#     dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
#     dks[1L] <- 1
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeromat)
#     gs <- list(zerovec)
#     order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
#     kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
#     for(k in 1L:kmax) {
#         Gs[which(order_array == k - 2L)] <- 0
#         gs[which(order_array == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
#             for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
#                 i3 <- k - i1 - i2
#                 if(fill_across[3L] && i1 + i2 + r > kmax) next
#                 if(fill_across[2L] && i1 + i3 + q > kmax) next
#                 # if(fill_across[1] && i2 + i3 + p > kmax) next
#                 tG <- zeromat
#                 if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
#                 if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
#                 if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
#                 Gs[[il3(i1, i2, i3)]] <- tG
#                 tg <- tG %*% mu
#                 if(i1 >= 1L) tg <- tg + A1 %*% gs[[il3(i1 - 1L, i2, i3)]]
#                 if(i2 >= 1L) tg <- tg + Gs[[il3(i1, i2 - 1L, i3)]] %*% mu + dks[i1 + 1L, i2, i3 + 1L] * mu + A2 %*% gs[[il3(i1, i2 - 1L, i3)]]
#                 if(i3 >= 1L) tg <- tg + Gs[[il3(i1, i2, i3 - 1L)]] %*% mu + dks[i1 + 1L, i2 + 1L, i3] * mu + A3 %*% gs[[il3(i1, i2, i3 - 1L)]]
#                 gs[[il3(i1, i2, i3)]] <- tg
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
#             }
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#         # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
#     }
#     return(dks)
# }

# hhat3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, p = m, q = m, r = m,
#                  fill_across = c(!missing(p), !missing(q), !missing(r))) { # , verbose = m > 200L) {
#     il3 <- function(i1, i2, i3) i1 + i2 * (p + 1L) + i3 * (p + 1L) * (q + 1L) + 1L
#     n <- length(L1)
#     zeros <- rep.int(0, n)
#     dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
#     dks[1L] <- 1
#     attr(dks, "logscale") <- 0
#     thr <- .Machine$double.xmax / 100 / n
#     Gs <- list(zeros)
#     gs <- list(zeros)
#     order_array <- outer(outer(0:p, 0:q, "+"), 0:r, "+")
#     kmax <- min(any(!fill_across) * m + sum(fill_across * c(p, q, r)), sum(dim(dks) - 1L))
#     for(k in 1L:kmax) {
#         Gs[which(order_array == k - 2L)] <- 0
#         gs[which(order_array == k - 2L)] <- 0
#         for(i1 in min(k, p):max(k - q - r, 0L, fill_across[1L] * (k + p - kmax))) {
#             for(i2 in min(k - i1, q):max(k - i1 - r, 0L)) {
#                 i3 <- k - i1 - i2
#                 if(fill_across[3L] && i1 + i2 + r > kmax) next
#                 if(fill_across[2L] && i1 + i3 + q > kmax) next
#                 # if(fill_across[1] && i2 + i3 + p > kmax) next
#                 tG <- zeros
#                 if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
#                 if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
#                 if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
#                 Gs[[il3(i1, i2, i3)]] <- tG
#                 tg <- tG * mu
#                 if(i1 >= 1L) tg <- tg + L1 * gs[[il3(i1 - 1L, i2, i3)]]
#                 if(i2 >= 1L) tg <- tg + (Gs[[il3(i1, i2 - 1L, i3)]] + dks[i1 + 1L, i2, i3 + 1L]) * mu + L2 * gs[[il3(i1, i2 - 1L, i3)]]
#                 if(i3 >= 1L) tg <- tg + (Gs[[il3(i1, i2, i3 - 1L)]] + dks[i1 + 1L, i2 + 1L, i3]) * mu + L3 * gs[[il3(i1, i2, i3 - 1L)]]
#                 gs[[il3(i1, i2, i3)]] <- tg
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
#             }
#         }
#         if(max(unlist(Gs)) > thr || max(unlist(gs)) > thr) {
#             dks <- dks / 1e10
#             Gs <- lapply(Gs, function(x) x / 1e10)
#             gs <- lapply(gs, function(x) x / 1e10)
#             attr(dks, "logscale") <- attr(dks, "logscale") - log(1e10)
#         }
#         # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
#     }
#     return(dks)
# }

##### hhat3_pjk_m #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
hhat3_pjk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, p = 1L) {
    n <- ncol(A1)
    p1 <- p + 1L
    m1 <- m + 1L
    In <- diag(n)
    dks <- array(0, dim = c(p1, m1, m1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    order_array <- outer(outer(rep.int(0, p1), 0:m, "+"), 0:m, "+")
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    Gc <- list()
    gc <- list()
    Gc[1L:p1] <- list(zeromat)
    gc[1L:p1] <- list(zerovec)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu + A1 %*% gc[[i]]
        dks[i + 1L, 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
        gc[[1L]] <- (tG + Go[[1L]][[1L]] + (dks[1L, k, 1L] * In)) %*% mu + A2 %*% go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            tG <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                  A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG + Go[[1L]][[i + 1L]]
                                + (dks[i + 1L, k, 1L] * In)) %*% mu +
                               A1 %*% gc[[i]] + A2 %*% go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                gc[[1L]] <- (tG + Go[[j + 1L]][[1L]] + Go[[j]][[1L]] + ((dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j]) * In)) %*% mu + A2 %*% go[[j + 1L]][[1L]] + A3 %*% go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
                for(i in 1L:p) {
                    tG <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                          A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                          A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG + Go[[j + 1L]][[i + 1L]] + Go[[j]][[i + 1L]]
                                     + ((dks[i + 1L, k - j, j + 1L] + dks[i + 1L, k - j + 1L, j]) * In)) %*% mu +
                                    A1 %*% gc[[i]] + A2 %*% go[[j + 1L]][[i + 1L]] + A3 %*% go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
        gc[[1L]] <- (tG + Go[[k]][[1L]] + (dks[1L, 1L, k] * In)) %*% mu + A3 %*% go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            tG <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                  A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG + Go[[k]][[i + 1L]]
                             + (dks[i + 1L, 1L, k] * In)) %*% mu +
                            A1 %*% gc[[i]] + A3 %*% go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gn)) > thr || max(unlist(gn)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}

##### hhat3_pjk_v #####
#' Coefficients in polynomial expansion of generating function---for
#' ratio with three matrices
#'
#' @rdname d3_ijk
#'
hhat3_pjk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, p = 1L) {
    n <- length(L1)
    p1 <- p + 1L
    m1 <- m + 1L
    dks <- array(0, dim = c(p1, m1, m1))
    dks[1L] <- 1
    lscf <- array(0, dim(dks))
    thr <- .Machine$double.xmax / 100 / n
    order_array <- outer(outer(rep.int(0, p1), 0:m, "+"), 0:m, "+")
    zeros <- rep.int(0, n)
    Gc <- list()
    gc <- list()
    Gc[1L:p1] <- list(zeros)
    gc[1L:p1] <- list(zeros)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] * mu + L1 * gc[[i]]
        dks[i + 1L, 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
        gc[[1L]] <- (tG + (Go[[1L]][[1L]] + dks[1L, k, 1L])) * mu + L2 * go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:p) {
            tG <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                  L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG + (Go[[1L]][[i + 1L]]
                                + dks[i + 1L, k, 1L])) * mu +
                               L1 * gc[[i]] + L2 * go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                gc[[1L]] <- (tG + (Go[[j + 1L]][[1L]] + Go[[j]][[1L]] + dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j])) * mu + L2 * go[[j + 1L]][[1L]] + L3 * go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
                for(i in 1L:p) {
                    tG <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                          L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                          L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG + (Go[[j + 1L]][[i + 1L]] + Go[[j]][[i + 1L]]
                                        + dks[i + 1L, k - j, j + 1L]
                                        + dks[i + 1L, k - j + 1L, j])) * mu +
                                       L1 * gc[[i]] + L2 * go[[j + 1L]][[i + 1L]] + L3 * go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
        gc[[1L]] <- (tG + (Go[[k]][[1L]] + dks[1L, 1L, k])) * mu + L3 * go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:p) {
            tG <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                  L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG + (Go[[k]][[i + 1L]]
                             + dks[i + 1L, 1L, k])) * mu +
                            L1 * gc[[i]] + L3 * go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gn)) > thr || max(unlist(gn)) > thr) {
            ind_dks <- which(order_array == k)
            ind_lscf <- which(order_array >= k)
            dks[ind_dks] <- dks[ind_dks] / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            lscf[ind_lscf] <- lscf[ind_lscf] - log(1e10)
        }
    }
    attr(dks, "logscale") <- lscf
    return(dks)
}
