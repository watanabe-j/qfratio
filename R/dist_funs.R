##### pqfr_A1B1 #####
#' Distribution Function of Ratio of Quadratic Forms
#'
#' \sQuote{Exact} distribution function of the ratio of quadratic forms in
#' normal variables,
#' \eqn{\frac{ (\mathbf{x^\mathit{T} A x}) }{ (\mathbf{x^\mathit{T} B x}) }
#'      }{ (x^T A x) / (x^T B x) }, where
#' \eqn{\mathbf{x} \sim N_n(\bm{0}_n, \mathbf{\Sigma})}{x ~ N_n(0_n, \Sigma)},
#' evaluated as a partial sum of infinite series including top-order zonal
#' polynomials (Forchini 2002).
#'
#' \code{quantile} should be of length one, because the cumulative distribution
#' is a complicated function of it and there is no practical way for
#' vectorization.
#'
#' @inheritParams qfrm
#'
#' @param quantile
#'   Quantile \eqn{q}; numeric of length one
#' @param Sigma
#'   Covariance matrix \eqn{\mathbf{\Sigma}}{\Sigma} for
#'   \eqn{\mathbf{x}}{x}
#' @param check_convergence
#'   Specifies how numerical convergence is checked (see \dQuote{Details} in
#'   \code{\link{qfrm}}).
#'   Options:
#'   \itemize{
#'     \item{\code{"relative"}: }{default; magnitude of the last term of
#'           the series relative to the sum is compared with \code{tol_conv}}
#'     \item{\code{"strict_relative"} or \code{TRUE}: }{same, but stricter than
#'           default by setting \code{tol_conv = .Machine$double.eps}
#'           (unless a smaller value is specified by the user)}
#'     \item{\code{"absolute"}: }{absolute magnitude of the last term is
#'           compared with \code{tol_conv}}
#'     \item{\code{"none"} or \code{FALSE}: }{skips convergence check}
#'   }
#'
#' @references
#' Forchini, G. (2002) The exact cumulative distribution function of
#'   a ratio of quadratic forms in normal variables, with application to
#'   the AR(1) model. *Econometric Theory*, **18**, 823--852.
#'   \doi{10.1017/S0266466602184015}.
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
pqfr_A1B1 <- function(quantile, A, B, m = 100L,
                      # mu = rep.int(0, n),
                      Sigma = diag(n),
                      check_convergence = c("relative", "strict_relative",
                                            "absolute", "none"),
                      use_cpp = TRUE,
                      tol_conv = .Machine$double.eps ^ (1/4),
                      tol_zero = .Machine$double.eps * 100,
                      tol_sing = .Machine$double.eps * 100,
                      thr_margin = 100) {
    if(isTRUE(check_convergence)) check_convergence <- "strict_relative"
    if(isFALSE(check_convergence)) check_convergence <- "none"
    check_convergence <- match.arg(check_convergence)
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
    # zeros <- rep.int(0, n)
    ## If Sigma is given, transform A, B, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma) && !iseq(Sigma, In, tol_zero)) {
        KiKS <- KiK(Sigma, tol_sing)
        K <- KiKS$K
        # iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        KtBK <- t(K) %*% B %*% K
        # iKmu <- iK %*% mu
        # ## If Sigma is singular, check conditions for A, B, mu, and Sigma
        # if(ncol(K) != n) {
        #     okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
        #             (iseq(A %*% mu, zeros, tol_zero) &&
        #              iseq(B %*% mu, zeros, tol_zero)) ||
        #             (iseq(crossprod(iK, KtAK %*% iK), A) &&
        #              iseq(crossprod(iK, KtBK %*% iK), B))
        #     if(!okay) {
        #         stop("For singular Sigma, certain condition must be met ",
        #              "for A, B, mu.\n  ",
        #              "Function for situations not satisfying this has not ",
        #              "developed.\n  See documentation for details")
        #     }
        # }
        return(pqfr_A1B1(quantile, KtAK, KtBK, m = m,
                         check_convergence = check_convergence,
                         tol_conv = tol_conv, tol_zero = tol_zero,
                         tol_sing = tol_sing))
    }
    ## Check basic requirements for arguments
    stopifnot(
        "A and B must be square matrices" = all(c(dim(A), dim(B)) == n),
        "quantile must be a length-one numeric" =
            is.numeric(quantile) && (length(quantile) == 1)
    )
    Ds <- eigen(A - quantile * B, symmetric = TRUE)$values
    D1 <- Ds[Ds > tol_sing]
    D2 <- -Ds[Ds < -tol_sing]
    n1 <- length(D1)
    n2 <- length(D2)
    if(n2 == 0) return(list(p = 0, terms = NULL))
    if(n1 == 0) return(list(p = 1, terms = NULL))
    if(use_cpp) {
        cppres <- A1B1_E(D1, D2, m, thr_margin)
        ansseq <- cppres$ansseq
    } else {
        D1i <- 1 / D1
        D2i <- 1 / D2
        sumtrDi <- sum(D1i) + sum(D2i)
        D1d <- D1i / sumtrDi
        D2d <- D2i / sumtrDi
        dk1 <- d1_i(rep.int(sum(D1d), n1) - D1d, m, thr_margin)
        dk2 <- d1_i(rep.int(sum(D2d), n2) - D2d, m, thr_margin)
        seq0m <- seq.int(0, m)
        ordmat <- outer(seq0m, seq0m, FUN = "+") + (n1 + n2) / 2
        ansmat <- lgamma(ordmat)
        ansmat <- ansmat - lgamma(n1 / 2 + 1 + seq0m) + log(abs(dk1))
        ansmat <- t(t(ansmat) - lgamma(n2 / 2 + seq0m) + log(abs(dk2)))
        ansmat <- ansmat + (sum(log(D1d)) + sum(log(D2d))) / 2
        ansmat <- ansmat - attr(dk1, "logscale")
        ansmat <- t(t(ansmat) - attr(dk2, "logscale"))
        ansmat <- exp(ansmat)
        # ansmat <- ansmat * sign(dk1)
        # ansmat <- t(t(ansmat) * sign(dk2))
        ansmat <- ansmat * gsl::hyperg_2F1(ordmat, 1, n1 / 2 + 1 + seq0m,
                                           sum(D1d))
        ansseq <- sum_counterdiag(ansmat)
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
