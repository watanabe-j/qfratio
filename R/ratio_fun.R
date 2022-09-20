##### qfrm #####
#' Moment of ratio of quadratic forms in normal variables
#'
#' @export
#'
qfrm <- function(A, B, p = 1, q = p, m = 100L, mu, Sigma,
                 tol_zero = .Machine$double.eps * 100,
                 tol_sing = .Machine$double.eps, ...) {
    ##
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
    if(missing(mu)) mu <- zeros
    ## If Sigma is given, transform A, B, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma)) {
        KiKS <- KiK(Sigma, tol_sing) # }, check_convergence = FALSE) {
        K <- KiKS$K
        iK <- KiKS$iK
        KtAK <- t(K) %*% A %*% K
        KtBK <- t(K) %*% B %*% K
        iKmu <- iK %*% mu
        ## If Sigma is singular, check conditions for A, B, mu, and Sigma
        if(ncol(K) != n) {
            okay <- (iseq(K %*% iKmu, mu, tol_zero)) ||
                    (iseq(A %*% mu, zeros, tol_zero) && iseq(B %*% mu, zeros, tol_zero)) ||
                    (iseq(crossprod(iK, KtAK %*% iK), A) &&
                     iseq(crossprod(iK, KtBK %*% iK), B))
            if(!okay) {
                stop("For singular Sigma, certain condition need to be met ",
                     "for A, B, or mu.\n  ",
                     "Function for situations not satisfying this has not ",
                     "developed.\n  See documentation for details.")
            }
        }
        return(qfrm(KtAK, KtBK, p, q, m = m, mu = iKmu, ...))
    }
    if(iseq(B, In, tol_zero)) {
        if((p %% 1) == 0 && p > 0) {
            return(qfrm_ApIq_int(A = A, p = p, q = q, mu = mu))
        } else {
            return(qfrm_ApIq_npi(A = A, p = p, q = q, m = m, mu = mu, ...))
        }
    } else {
        if(iseq(A, In, tol_zero)) {
            return(qfrm_ApIq_npi(A = B, p = -q, q = -p, m = m, mu = mu, ...))
        }
    }
    if((p %% 1) == 0) {
        return(qfrm_ApBq_int(A = A, B = B, p = p, q = q, m = m, mu = mu, ...))
    } else {
        return(qfrm_ApBq_npi(A = A, B = B, p = p, q = q, m = m, mu = mu, ...))
    }
}
##### qfmrm #####
#' Moment of multiple ratio of quadratic forms in normal variables
#'
#' @export
#'
qfmrm <- function(A, B, D, p = 1, q = p / 2, r = p / 2, m = 100L, mu, Sigma,
                 tol_zero = .Machine$double.eps * 100,
                 tol_sing = .Machine$double.eps, ...) {
    ##
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
    if(missing(mu)) mu <- zeros
    ## If any pair of the three arguments are equal,
    ## reduce the problem to a simple ratio
    if(iseq(B, D, tol_zero)) {
        return(qfrm(A, B, p, q + r, m = m, mu = mu, Sigma = Sigma, ...))
    }
    if(iseq(A, D, tol_zero) && p >= r) {
        return(qfrm(A, B, p - r, q, m = m, mu = mu, Sigma = Sigma, ...))
    }
    if(iseq(A, B, tol_zero) && p >= q) {
        return(qfrm(A, D, p - r, r, m = m, mu = mu, Sigma = Sigma, ...))
    }
    ## If Sigma is given, transform A, B, D, and mu, and
    ## call this function recursively with new arguments
    if(!missing(Sigma)) {
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
                    (iseq(A %*% mu, zeros, tol_zero) && iseq(B %*% mu, zeros, tol_zero)
                          && iseq(D %*% mu, zeros, tol_zero)) ||
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
        return(qfmrm(KtAK, KtBK, KtDK, p, q, r, m = m, mu = iKmu, ...))
    }
    if(iseq(A, In, tol_zero)) {
        return(qfmrm_IpBDqr_gen(B = B, D = D, p = p, q = q, r = r, m = m, mu = mu, ...))
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
            return(qfmrm_ApBIqr_int(A = A, B = B, p = p, q = q, r = r, m = m, mu = mu, ...))
        } else {
            return(qfmrm_ApBIqr_npi(A = A, B = B, p = p, q = q, r = r, m = m, mu = mu, ...))
        }
    }
    if((p %% 1) == 0) {
        return(qfmrm_ApBDqr_int(A = A, B = B, D = D, p = p, q = q, r = r, m = m, mu = mu, ...))
    } else {
        return(qfmrm_ApBDqr_npi(A = A, B = B, D = D, p = p, q = q, r = r, m = m, mu = mu, ...))
    }
}


###############################
## Function for positive integer moment of a quadratic form
###############################
#' @export
#'
qfm_Ap_int <- function(A, p = 1, mu = rep.int(0, n),
                          use_cpp = FALSE, cpp_method = "Eigen",
                             tol_zero = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    n <- ncol(A)
    stopifnot(
        "A should be a square matrix" = all(c(dim(A)) == n),
        "p should be a positive integer" = {
            length(p) == 1 &&
            (p %% 1) == 0 &&
            p >= 1
        },
        "mu should be an n-vector" = length(mu) == n
    )
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
    ansseq <- ans
    errseq <- 0
    attr(errseq, "exact") <- TRUE
    errorb <- errseq
    structure(list(statistic = ans, error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}

###############################
## Function for product of quadratic forms
###############################
#' @export
#'
qfpm_ABpq_int <- function(A, B, p = 1, q = 1, mu = rep.int(0, n),
                          use_cpp = FALSE, cpp_method = "Eigen",
                          tol_zero = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
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
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero)
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
                # LA <- diag(A)
                # dpq <- d2_ij_v(LA, LB, m = p + q, p, q)[p + 1, q + 1]
                dpq <- d2_pj_v(LA, LB, m = p + q, p)[p + 1, q + 1]
            } else {
                # dpq <- d2_ij_m(A, diag(LB), m = p + q, p, q)[p + 1, q + 1]
                dpq <- d2_pj_m(A, diag(LB), m = p + q, p)[p + 1, q + 1]
            }
        } else {
            if(use_vec) {
                # LA <- diag(A)
                # dpq <- dtil2_ij_v(LA, LB, mu, m = p + q, p, q)[p + 1, q + 1]
                dpq <- dtil2_pq_v(LA, LB, mu, p, q)[p + 1, q + 1]
            } else {
                # dpq <- dtil2_ij_m(A, diag(LB), mu, m = p + q, p, q)[p + 1, q + 1]
                dpq <- dtil2_pq_m(A, diag(LB), mu, p, q)[p + 1, q + 1]
            }
        }
        ans <- 2 ^ (p + q) * factorial(p) * factorial(q) * dpq
    }
    ansseq <- ans
    errseq <- 0
    attr(errseq, "exact") <- TRUE
    errorb <- errseq
    structure(list(statistic = ans, error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}

#' @export
#'
qfpm_ABDpqr_int <- function(A, B, D, p = 1, q = 1, r = 1, mu = rep.int(0, n),
                            use_cpp = FALSE, cpp_method = "Eigen",
                             tol_zero = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
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
        "A and B should be square matrices" = all(c(dim(A), dim(B), dim(D)) == n),
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
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    ## Rotate A and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero) && is_diagonal(D, tol_zero)
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
                # LA <- diag(A)
                # LD <- diag(D)
                # dpqr <- d3_ijk_v(LA, LB, LD, m = p + q + r, p, q, r)[p + 1, q + 1, r + 1]
                dpqr <- d3_pjk_v(LA, LB, LD, q + r, p)[p + 1, q + 1, r + 1]
            } else {
                # dpqr <- d3_ijk_m(A, diag(LB), D, m = p + q + r, p, q, r)[p + 1, q + 1, r + 1]
                dpqr <- d3_pjk_m(A, diag(LB), D, q + r, p)[p + 1, q + 1, r + 1]
            }
        } else {
            if(use_vec) {
                # LA <- diag(A)
                # LD <- diag(D)
                # dpqr <- dtil3_ijk_v(LA, LB, LD, mu, m = p + q + r, p, q, r)[p + 1, q + 1, r + 1]
                dpqr <- dtil3_pqr_v(LA, LB, LD, mu, p, q, r)[p + 1, q + 1, r + 1]
            } else {
                # dpqr <- dtil3_ijk_m(A, diag(LB), D, mu, m = p + q + r, p, q, r)[p + 1, q + 1, r + 1]
                dpqr <- dtil3_pqr_m(A, diag(LB), D, mu, p, q, r)[p + 1, q + 1, r + 1]
            }
        }
        ans <- 2 ^ (p + q + r) * factorial(p) * factorial(q) * factorial(r) * dpqr
    }
    ansseq <- ans
    errseq <- 0
    attr(errseq, "exact") <- TRUE
    errorb <- errseq
    structure(list(statistic = ans, error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}




###############################
## Functions for multiple ratio
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
#' \code{qfrm_ApIq_int()}: Positive integer moment of
#' \eqn{(\mathbf{x}^T \mathbf{A} \mathbf{x}) / (\mathbf{x}^T \mathbf{x})},
#' where \eqn{\mathbf{x}} is a vector of independent standard normal variables.
#'
# #' @importFrom gsl hyperg_1F1
#'
#' @export
#'
qfrm_ApIq_int <- function(A, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                          use_cpp = FALSE, cpp_method = "Eigen",
                          tol_zero = .Machine$double.eps * 100) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
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
            ans <- exp((p - q) * log(2) + lgamma(p + 1) + lgamma(n / 2 + p - q)
                       - lgamma(n / 2 + p)) * dp
            ansseq <- ans
        } else {
            mu <- c(mu)
            if(requireNamespace("gsl", quietly = TRUE)) {
                ## This is an exact expression (Hillier et al. 2014, (58))
                D <- c(crossprod(eigA$vectors, mu)) ^ 2
                aps <- arl(LA, D, m = p)[p + 1, ]
                ls <- 0:p
                ansseq <- exp((p - q) * log(2) + lgamma(p + 1)
                            + lgamma(n / 2 + p - q + ls) - ls * log(2)
                            - lgamma(ls + 1) - lgamma(n / 2 + p + ls)) *
                          gsl::hyperg_1F1(q, n / 2 + p + ls, -crossprod(mu) / 2) * aps
            } else {
                ## This is a recursive alternative (Hillier et al. 2014, (53))
                ## which is less accurate (by truncation) and slower
                warning("This expression is suboptimal with a nonzero mean.\n  ",
                        "An exact alternative is available with ",
                        "the package \"gsl\" (recommended).")
                dks <- d2_ij_m(A, tcrossprod(mu), m, m1 = p)[p + 1, ]
                ansseq <- exp((p - q) * log(2) + lgamma(1 + p) - c(crossprod(mu)) / 2
                              + lgamma(n / 2 + p - q + 0:m) - 0:m * log(2)
                              - lgamma(1/2 + 0:m) + lgamma(1/2) - lgamma(n / 2 + p + 0:m)
                              + log(dks))
                exact <- FALSE
            }
            ans <- sum(ansseq)
        }
    }
    if(exact) {
        errseq <- ans - cumsum(ansseq)
        errorb <- errseq[length(errseq)]
        attr(errorb, "exact") <- TRUE
        attr(errseq, "exact") <- TRUE
    } else {
        errorb <- NA_real_
        errseq <- NULL
    }
    structure(list(statistic = ans, error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}

##### qfrm_ApIq_npi #####
#' Non-positive-integer moment of ratio of quadratic forms in normal variables
#'
#' \code{qfrm_ApIq_npi()}: Non-positive-integer moment of
#' \eqn{(\mathbf{x}^T \mathbf{A} \mathbf{x}) / (\mathbf{x}^T \mathbf{x})},
#' where \eqn{\mathbf{x}} is a vector of independent standard normal variables.
#'
#' @export
#'
qfrm_ApIq_npi <- function(A, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                    alpha1 = 1,
                    error_bound = TRUE, check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    n <- ncol(A)
    stopifnot(
        "A should be a square matrix" = all(c(dim(A)) == n),
        # "mu should be an n-vector" = length(mu) == n,
        "p should be a real number" = length(p) == 1,
        "q should be a real number" = length(q) == 1
    )
    if((p %% 1) == 0 && p > 0) {
        warning("For integral p, qfrm_ApIq_int() works better")
    }
    A <- (A + t(A)) / 2
    eigA <- eigen(A, symmetric = TRUE)
    LA <- eigA$values
    UA <- eigA$vectors
    nndefA <- all(LA >= 0)
    if(!nndefA && ((p %% 1) != 0 || p < 0)) {
        warning("(Numerically) negative eigenvalue(s) detected.\n  ",
                "Ensure A is nonnegative definite, as non-positive-integer ",
                "moment of\n  the quadratic form is not well defined otherwise.")
    }
    ## Check condition for existence of moment (Bao & Kan, 2013, prop. 1)
    cond_exist <- n / 2 + p > q ## condition(1)
    stopifnot("Moment does not exist in this combination of p, q, and rank(B)" = cond_exist)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    b1 <- alpha1 / max(abs(LA))
    if(use_cpp) {
        if(central) {
            cppres <- ApIq_npi_cvE(LA, b1, p, q, m, error_bound)
        } else {
            cppres <- ApIq_npi_nvE(LA, UA, b1, mu, p, q, m)
        }
        ansseq <- cppres$ansseq
    } else {
        LAh <- rep.int(1, n) - b1 * LA
        if(central) {
            dks <- d1_i(LAh, m = m)
            ansseq <- hgs_1d(dks, -p, n / 2, ((p - q) * log(2) - p * log(b1)
                             + lgamma(n / 2 + p - q) - lgamma(n / 2)))
        } else {
            mu <- c(crossprod(UA, c(mu)))
            # ## This is based on recursion for d as in Hillier et al. (2009)
            # dks <- d2_ij_m(diag(LAh), tcrossprod(mu), m)
            # ansmat <- hgs_dmu_2d(dks, -p, n / 2 + p - q, n / 2,
            #                      ((p - q) * log(2) - c(crossprod(mu)) / 2
            #                      - p * log(b1) + lgamma(n / 2 + p - q) - lgamma(n / 2)))
            ## This is based on recursion for h as in Hillier et al. (2014)
            dks <- h2_ij_v(LAh, rep.int(0, n), mu, m)
            ansmat <- hgs_2d(dks, -p, q, n / 2, ((p - q) * log(2) - p * log(b1)
                             + lgamma(n / 2 + p - q) - lgamma(n / 2)))
            ansseq <- sum_counterdiag(ansmat)
        }
        scf <- attr(dks, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence) {
        # try(plot(seq_along(ansseq) - 1L, cumsum(ansseq), type = "l", col = "royalblue4",
        #      ylim = sum(ansseq) * c(0.9, 1.1), xlab = "Order of polynomials",
        #      ylab = "Moment of ratio"))
        if(abs(ansseq[length(ansseq)]) > tol_conv) {
            warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
        }
    }
    if(error_bound && central) {
        if(use_cpp) {
            errseq <- cppres$errseq
            twosided <- FALSE
        } else {
            # if(!nndefA && (p %% 2) == 1) {
            #     Lp <- abs(LAh)
            #     dkst <- d1_i(Lp, m = m)
            #     twosided <- TRUE
            # } else {
            Lp <- LAh
            dkst <- dks
            twosided <- FALSE
            # }
            lcoefe <- (lgamma(-p + 0:m + 1) - lgamma(-p)
                       - lgamma(n / 2 + 0:m + 1) + lgamma(n / 2 + p - q)
                       + (p - q) * log(2) - p * log(b1))
            errseq <- exp(lcoefe - sum(log(1 - Lp)) / 2) - exp((lcoefe + log(cumsum(dkst[1:(m + 1)]))) - log(scf))
            errseq <- errseq * cumprod(sign(-p + 0:m))
        }
        errorb <- errseq[length(errseq)]
        attr(errseq, "twosided") <- twosided
        attr(errorb, "twosided") <- twosided
        if(any(LA < tol_sing)) {
            warning("Argument matrix is numerically close to singular.\n  ",
            "If it is singular, this error bound is invalid.")
            attr(errseq, "singular") <- TRUE
            attr(errorb, "singular") <- TRUE
        }
        if(alpha1 > 1) {
            warning("Error bound is unreliable when alpha1 > 1\n  ",
            "It is returned purely for heuristic purpose")
            attr(errseq, "alphaout") <- TRUE
            attr(errorb, "alphaout") <- TRUE
        }
        # if(check_convergence) {
        #     lines(errseq + cumsum(ansseq), col = "tomato", lty = 2)
        # }
    } else {
        if(error_bound) {
            warning("Error bound is unavailable for qfrm_ApIq_npi ",
                    "when mu is nonzero")
        }
        errseq <- NULL
        errorb <- NULL
    }
    structure(list(statistic = sum(ansseq), error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}


#'
#'
#' @export
#'
qfrm_ApBq_int <- function(A, B, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                    alpha2 = 1,
                    # fun = c("dki1", "dk2"),
                    error_bound = TRUE, check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    ## If A or B is missing, let it be an identity matrix
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
    }
    if(missing(B)) B <- In
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
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if(iseq(B, In, tol_zero)) {
        warning("For B = I, qfrm_ApIq_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
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
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "Moment does not exist in this combination of p, q, and rank(B)" = cond_exist)
    use_vec <- is_diagonal(A, tol_zero)
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
                cppres <- ApBq_int_cvE(LA, LB, b2, p, q, m, error_bound)
            } else {
                cppres <- ApBq_int_cmE(A, LA, UA, LB, b2, p, q, m, error_bound)
            }
        } else {
            if(use_vec) {
                cppres <- ApBq_int_nvE(LA, LB, b2, mu, p, q, m, error_bound)
            } else {
                cppres <- ApBq_int_nmE(A, LA, UA, LB, b2, mu, p, q, m, error_bound)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(use_vec) {
            # UA <- In
            # LA <- diag(A)
            LBh <- rep.int(1, n) - b2 * LB
            if(central) {
                # dks <- d2_ij_v(LA, LBh, m, m1 = p)[p + 1, 1:(m + 1)]
                dksm <- d2_pj_v(LA, LBh, m, m1 = p)
            } else {
                # dks <- htil2_ij_v(LA, LBh, mu, m, m1 = p)[p + 1, 1:(m + 1)]
                dksm <- htil2_pj_v(LA, LBh, mu, m, m1 = p)
            }
        } else {
            # eigA <- eigen(A, symmetric = TRUE)
            # UA <- eigA$vectors
            # LA <- eigA$values
            Bh <- In - b2 * diag(LB)
            if(central) {
                # dks <- d2_ij_m(A, Bh, m, m1 = p)[p + 1, 1:(m + 1)]
                dksm <- d2_pj_m(A, Bh, m, m1 = p)
            } else {
                # dks <- htil2_ij_m(A, Bh, mu, m, m1 = p)[p + 1, 1:(m + 1)]
                dksm <- htil2_pj_m(A, Bh, mu, m, m1 = p)
            }
        }
        dks <- dksm[p + 1, 1:(m + 1)]
        ansseq <- hgs_1d(dks, q, n / 2 + p, ((p - q) * log(2) + q * log(b2)
                         + lgamma(p + 1) + lgamma(n / 2 + p - q) - lgamma(n / 2 + p)))
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence) {
        # try(plot(seq_along(ansseq) - 1L, cumsum(ansseq), type = "l", col = "royalblue4",
        #      ylim = sum(ansseq) * c(0.9, 1.1), xlab = "Order of polynomials",
        #      ylab = "Moment of ratio"))
        if(abs(ansseq[length(ansseq)]) > tol_conv) {
            warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
        }
    }
    if(error_bound) {
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
                        dkstm <- d2_ij_v(LAp, LBh, m, m1 = p)
                    } else {
                        Ap <- S_fromUL(UA, LAp)
                        dkstm <- d2_ij_m(Ap, Bh, m, m1 = p)
                    }
                } else {
                    # LAp <- LA
                    Ap <- A
                    dkstm <- dksm
                }
                if(use_vec) {
                    dp <- d1_i(LAp / LB / b2, p)[p + 1]
                } else {
                    # Bisqr <- S_fromUL(In, 1 / sqrt(LB))
                    # dp <- d1_i(eigen(Bisqr %*% Ap %*% Bisqr, symmetric = TRUE)$values / b2, p)[p + 1]
                    Bisqr <- 1 / sqrt(LB)
                    dp <- d1_i(eigen(t(t(Ap * Bisqr) * Bisqr), symmetric = TRUE)$values / b2, p)[p + 1]
                }
            } else {
                twosided <- TRUE
                mub <- sqrt(2 / b2) * mu / sqrt(LB)
                deldif2 <- (sum(mub ^ 2) - sum(mu ^ 2)) / 2
                if(use_vec) {
                    # dkst <- hhat2_ij_v(LAp, LBh, mu, m, m1 = p)[p + 1, 1:(m + 1)]
                    dkstm <- hhat2_pj_v(LAp, LBh, mu, m, m1 = p)
                    dp <- dtil1_i_v(LAp / LB / b2, mub, p)[p + 1]
                } else {
                    Ap <- S_fromUL(UA, LAp)
                    # dkst <- hhat2_ij_m(Ap, Bh, mu, m, m1 = p)[p + 1, 1:(m + 1)]
                    dkstm <- hhat2_pj_m(Ap, Bh, mu, m, m1 = p)
                    # Bisqr <- S_fromUL(In, 1 / sqrt(LB))
                    # dp <- dtil1_i_m(Bisqr %*% Ap %*% Bisqr / b2, mub, p)[p + 1]
                    Bisqr <- 1 / sqrt(LB)
                    dp <- dtil1_i_m(t(t(Ap * Bisqr) * Bisqr) / b2, mub, p)[p + 1]
                }
            }
            dkst <- dkstm[p + 1, 1:(m + 1)]
            scft <- attr(dkstm, "scale")
            lBdet <- sum(log(LB * b2))
            lcoefe <- (lgamma(q + 0:m + 1) - lgamma(q)
                        - lgamma(n / 2 + p + 0:m + 1) + lgamma(n / 2 + p - q)
                        + (p - q) * log(2) + q * log(b2) + lgamma(p + 1))
            errseq <- exp(lcoefe + (deldif2 + log(dp) - lBdet / 2)) -
            exp(lcoefe + log(cumsum(dkst)) - log(scft))
        }
        errorb <- errseq[length(errseq)]
        attr(errseq, "twosided") <- twosided
        attr(errorb, "twosided") <- twosided
        if(any(LB < tol_sing)) {
            warning("Argument matrix B is numerically close to singular.\n  ",
                    "If it is singular, this error bound is invalid.")
            attr(errseq, "singular") <- TRUE
            attr(errorb, "singular") <- TRUE
        }
        if(alpha2 > 1) {
            warning("Error bound is unreliable ",
                    "when alpha2 > 1\n  ",
                    "It is returned purely for heuristic purpose")
            attr(errseq, "alphaout") <- TRUE
            attr(errorb, "alphaout") <- TRUE
        }
        # if(check_convergence) {
        #     lines(errseq + cumsum(ansseq), col = "tomato", lty = 2)
        # }
    } else {
        errseq <- NULL
        errorb <- NULL
    }
    structure(list(statistic = sum(ansseq), error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}

#' @export
#'
qfrm_ApBq_npi <- function(A, B, p = 1, q = p, m = 100L, mu = rep.int(0, n),
                    alpha1 = 1, alpha2 = 1,
                    # fun = c("dki1", "dk2"),
                    check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    ## If A or B is missing, let it be an identity matrix
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
    }
    if(missing(B)) B <- In
    ## Check basic requirements for arguments
    stopifnot(
        "A and B should be square matrices" = all(c(dim(A), dim(B)) == n),
        "p and q should be nonnegative real numbers" = {
            length(p) == 1 &&
            length(q) == 1 &&
            p >= 0 &&
            q >= 0
        },
        "alpha1 should be a scalar with 0 < alpha1 < 2" = {
            is.numeric(alpha1) &&
            length(alpha1) == 1 &&
            alpha1 > 0 &&
            alpha1 < 2
        },
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if((p %% 1) == 0) {
        warning("For integral p, qfrm_ApBq_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
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
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "Moment does not exist in this combination of p, q, and rank(B)" = cond_exist)
    use_vec <- is_diagonal(A, tol_zero)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    LA <- if(use_vec) diag(A) else eigen(A, symmetric = TRUE)$values
    b1 <- alpha1 / max(abs(LA))
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ApBq_npi_cvE(LA, LB, b1, b2, p, q, m)
            } else {
                cppres <- ApBq_npi_cmE(A, LB, b1, b2, p, q, m)
            }
        } else {
            if(use_vec) {
                cppres <- ApBq_npi_nvE(LA, LB, b1, b2, mu, p, q, m)
            } else {
                cppres <- ApBq_npi_nmE(A, LB, b1, b2, mu, p, q, m)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(use_vec) {
            # LA <- diag(A)
            # b1 <- alpha1 / max(abs(LA))
            LAh <- rep.int(1, n) - b1 * LA
            LBh <- rep.int(1, n) - b2 * LB
            if(central) {
                dksm <- d2_ij_v(LAh, LBh, m)
            } else {
                dksm <- h2_ij_v(LAh, LBh, mu, m)
            }
        } else {
            # eigA <- eigen(A, symmetric = TRUE)
            # LA <- eigA$values
            # b1 <- alpha1 / max(abs(LA))
            Ah <- In - b1 * A
            Bh <- In - b2 * diag(LB)
            if(central) {
                dksm <- d2_ij_m(Ah, Bh, m)
            } else {
                dksm <- h2_ij_m(Ah, Bh, mu, m)
            }
        }
        ansmat <- hgs_2d(dksm, -p, q, n / 2, ((p - q) * log(2)
                         - p * log(b1) + q * log(b2) + lgamma(n / 2 + p - q) -
                         lgamma(n / 2)))
        ansseq <- sum_counterdiag(ansmat)
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence) {
        # try(plot(seq_along(ansseq) - 1L, cumsum(ansseq), type = "l", col = "royalblue4",
        #      ylim = sum(ansseq) * c(0.9, 1.1), xlab = "Order of polynomials",
        #      ylab = "Moment of ratio"))
        if(abs(ansseq[length(ansseq)]) > tol_conv) {
            warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
        }
    }
    structure(list(statistic = sum(ansseq), error_bound = NA_real_,
                   res_seq = ansseq, err_seq = NULL), class = "qfrm")
}



###############################
## Functions for multiple ratio
###############################

#'
#'
#' @export
#'
qfmrm_ApBIqr_int <- function(A, B, p = 1, q = 1, r = 1, m = 100L,
                    mu = rep.int(0, n), alpha2 = 1,
                    # fun = c("dki1", "dk2"),
                    error_bound = TRUE, check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    ## If A or B is missing, let it be an identity matrix
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
    }
    if(missing(B)) B <- In
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
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
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
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "Moment does not exist in this combination of p, q, r, and rank(B)" = cond_exist)
    use_vec <- is_diagonal(A, tol_zero)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    if(use_vec) {
        LA <- diag(A)
        LBh <- rep.int(1, n) - b2 * LB
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        UA <- eigA$vectors
        LA <- eigA$values
        Bh <- In - b2 * diag(LB)
    }
    if(use_cpp) {
        # if(cpp_method == "arma") {
        #     ansseq <- ApBDqr_int_vA(LA, LBh, b2, p, q, r, m)
        # } else {
        # }
        if(central) {
            if(use_vec) {
                cppres <- ApBIqr_int_cvE(LA, LB, b2, p, q, r, m, error_bound)
            } else {
                cppres <- ApBIqr_int_cmE(A, LA, UA, LB, b2, p, q, r, m, error_bound)
            }
        } else {
            if(use_vec) {
                cppres <- ApBIqr_int_nvE(LA, LB, b2, mu, p, q, r, m, error_bound)
            } else {
                cppres <- ApBIqr_int_nmE(A, LA, UA, LB, b2, mu, p, q, r, m, error_bound)
            }
        }
        # browser()
        ansseq <- cppres$ansseq
    } else {
        if(central) {
            if(use_vec) {
                # dks <- d2_ij_v(LA, LBh, m, m1 = p)[p + 1, 1:(m + 1)]
                dksm <- d2_pj_v(LA, LBh, m, m1 = p)
            } else {
                # dks <- d2_ij_m(A, Bh, m, m1 = p)[p + 1, 1:(m + 1)]
                dksm <- d2_pj_m(A, Bh, m, m1 = p)
            }
            dks <- dksm[p + 1, 1:(m + 1)]
            ansseq <- hgs_1d(dks, q, n / 2 + p, ((p - q - r) * log(2) + q * log(b2)
                    + lgamma(p + 1) + lgamma(n / 2 + p - q - r) - lgamma(n / 2 + p)))
        } else {
            if(use_vec) {
                # dksm <- htil3_ijk_v(LA, LBh, rep.int(0, n), mu, m, m1 = p)[p + 1, , ]
                dksm <- htil3_pjk_v(LA, LBh, rep.int(0, n), mu, m, m1 = p)
            } else {
                # dksm <- htil3_ijk_m(A, Bh, matrix(0, n, n), mu, m, m1 = p)[p + 1, , ]
                dksm <- htil3_pjk_m(A, Bh, matrix(0, n, n), mu, m, m1 = p)
            }
            dks <- dksm[p + 1, , ]
            ansmat <- hgs_2d(dks, q, r, n / 2 + p, ((p - q - r) * log(2)
                            + q * log(b2) + lgamma(p + 1)
                            + lgamma(n / 2 + p - q - r) - lgamma(n / 2 + p)))
            ansseq <- sum_counterdiag(ansmat)
        }
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence && abs(ansseq[length(ansseq)]) > tol_conv) {
        warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
    }
    if(error_bound) {
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
                        # dkst <- d2_ij_v(LAp, LBh, m, m1 = p)[p + 1, 1:(m + 1)]
                        dkstm <- d2_pj_v(LAp, LBh, m, m1 = p)
                    } else {
                        Ap <- S_fromUL(UA, LAp)
                        # dkst <- d2_ij_m(Ap, Bh, m, m1 = p)[p + 1, 1:(m + 1)]
                        dkstm <- d2_pj_m(Ap, Bh, m, m1 = p)
                    }
                } else {
                    # LAp <- LA
                    Ap <- A
                    dkstm <- dksm
                }
                dkst <- dkstm[p + 1, 1:(m + 1)]
                if(use_vec) {
                    dp <- d1_i(LAp / LB / b2, p)[p + 1]
                } else {
                    # Bisqr <- S_fromUL(In, 1 / sqrt(LB))
                    # dp <- d1_i(eigen(Bisqr %*% Ap %*% Bisqr, symmetric = TRUE)$values / b2, p)[p + 1]
                    Bisqr <- 1 / sqrt(LB)
                    dp <- d1_i(eigen(t(t(Ap * Bisqr) * Bisqr), symmetric = TRUE)$values / b2, p)[p + 1]
                }
            } else {
                twosided <- TRUE
                mub <- sqrt(3 / b2) * mu / sqrt(LB)
                deldif2 <- (sum(mub ^ 2) - sum(mu ^ 2)) / 2
                s <- max(q, r)
                if(use_vec) {
                    # dksmt <- hhat3_ijk_v(LAp, LBh, rep.int(0, n), mu, m, m1 = p)[p + 1, , ]
                    dkstm <- hhat3_pjk_v(LAp, LBh, rep.int(0, n), mu, m, m1 = p)
                    dp <- dtil1_i_v(LAp / LB / b2, mub, p)[p + 1]
                } else {
                    Ap <- S_fromUL(UA, LAp)
                    # dksmt <- hhat3_ijk_m(Ap, Bh, matrix(0, n, n), mu, m, m1 = p)[p + 1, , ]
                    dkstm <- hhat3_pjk_m(Ap, Bh, matrix(0, n, n), mu, m, m1 = p)
                    # Bisqr <- S_fromUL(In, 1 / sqrt(LB))
                    # dp <- dtil1_i_m(Bisqr %*% Ap %*% Bisqr / b2, mub, p)[p + 1]
                    Bisqr <- 1 / sqrt(LB)
                    dp <- dtil1_i_m(t(t(Ap * Bisqr) * Bisqr) / b2, mub, p)[p + 1]
                }
                dkst <- sum_counterdiag(dkstm[p + 1, , ])
            }
            scft <- attr(dkstm, "scale")
            lBdet <- sum(log(LB * b2))
            lcoefe <- (lgamma(s + 0:m + 1) - lgamma(s)
                       - lgamma(n / 2 + p + 0:m + 1) + lgamma(n / 2 + p - q - r)
                       + (p - q - r) * log(2) + q * log(b2) + lgamma(p + 1))
            errseq <- exp(lcoefe + (deldif2 + log(dp) - lBdet / 2)) -
                      exp(lcoefe + log(cumsum(dkst)) - log(scft))
        }
        errorb <- errseq[length(errseq)]
        attr(errseq, "twosided") <- twosided
        attr(errorb, "twosided") <- twosided
        if(any(LB < tol_sing)) {
            warning("Argument matrix B is numerically close to singular.\n  ",
                    "If it is singular, this error bound is invalid.")
            attr(errseq, "singular") <- TRUE
            attr(errorb, "singular") <- TRUE
        }
        if(alpha2 > 1) {
            warning("Error bound is unreliable ",
                    "when alpha2 > 1\n  ",
                    "It is returned purely for heuristic purpose")
            attr(errseq, "alphaout") <- TRUE
            attr(errorb, "alphaout") <- TRUE
        }
    } else {
        errseq <- NULL
        errorb <- NULL
    }
    structure(list(statistic = sum(ansseq), error_bound = errorb,
                   res_seq = ansseq, err_seq = errseq), class = "qfrm")
}

#'
#'
#' @export
#'
qfmrm_ApBIqr_npi <- function(A, B, p = 1, q = 1, r = 1, m = 100L,
                    mu = rep.int(0, n), alpha1 = 1, alpha2 = 1,
                    # fun = c("dki1", "dk2"),
                    use_cpp = FALSE, cpp_method = "Eigen",
                    check_convergence = TRUE,
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    ## If A or B is missing, let it be an identity matrix
    if(missing(A)) {
        if(missing(B)) stop("Provide at least one of A and B")
        n <- dim(B)[1L]
        In <- diag(n)
        A <- In
    } else {
        n <- dim(A)[1L]
        In <- diag(n)
    }
    if(missing(B)) B <- In
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
        "alpha1 should be a scalar with 0 < alpha1 < 2" = {
            is.numeric(alpha1) &&
            length(alpha1) == 1 &&
            alpha1 > 0 &&
            alpha1 < 2
        },
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if((p %% 1) == 0) {
        warning("For integral p, qfmrm_ApBIqr_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
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
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "Moment does not exist in this combination of p, q, r, and rank(B)" = cond_exist)
    use_vec <- is_diagonal(A, tol_zero)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    if(use_vec) {
        LA <- diag(A)
        b1 <- alpha1 / max(abs(LA))
        LAh <- rep.int(1, n) - b1 * LA
        LBh <- rep.int(1, n) - b2 * LB
        # if(central) {
        #     dksm <- d2_ij_v(LAh, LBh, m)
        # } else {
        #     dksm <- h2_ij_v(LAh, LBh, mu, m)
        # }
    } else {
        eigA <- eigen(A, symmetric = TRUE)
        LA <- eigA$values
        b1 <- alpha1 / max(abs(LA))
        Ah <- In - b1 * A
        Bh <- In - b2 * diag(LB)
        # if(central) {
        #     dksm <- d2_ij_m(Ah, Bh, m)
        # } else {
        #     dksm <- h2_ij_m(Ah, Bh, mu, m)
        # }
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ApBIqr_npi_cvE(LA, LB, b1, b2, p, q, r, m)
            } else {
                cppres <- ApBIqr_npi_cmE(A, LB, b1, b2, p, q, r, m)
            }
        } else {
            if(use_vec) {
                cppres <- ApBIqr_npi_nvE(LA, LB, b1, b2, mu, p, q, r, m)
            } else {
                cppres <- ApBIqr_npi_nmE(A, LB, b1, b2, mu, p, q, r, m)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(central) {
            if(use_vec) {
                dksm <- d2_ij_v(LAh, LBh, m)
            } else {
                dksm <- d2_ij_m(Ah, Bh, m)
            }
            ansmat <- hgs_2d(dksm, -p, q, n / 2, ((p - q - r) * log(2)
                        - p * log(b1) + q * log(b2) + lgamma(n / 2 + p - q - r)
                        - lgamma(n / 2)))
            ansseq <- sum_counterdiag(ansmat)
        } else {
            if(use_vec) {
                dksm <- h3_ijk_v(LAh, LBh, rep.int(0, n), mu, m)
            } else {
                dksm <- h3_ijk_m(Ah, Bh, matrix(0, n, n), mu, m)
            }
            ansarr <- hgs_3d(dksm, -p, q, r, n / 2, ((p - q - r) * log(2)
                        - p * log(b1) + q * log(b2) + lgamma(n / 2 + p - q - r)
                        - lgamma(n / 2)))
            ansseq <- sum_counterdiag3D(ansarr)
        }
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence && abs(ansseq[length(ansseq)]) > tol_conv) {
        warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
    }
    structure(list(statistic = sum(ansseq), error_bound = NA_real_,
                   res_seq = ansseq, err_seq = NULL), class = "qfrm")
}

## WIP: Condition(s) for existence of moment
#'
#'
#' @export
#'
qfmrm_IpBDqr_gen <- function(B, D, p = 1, q = 1, r = 1, mu = rep.int(0, n),
                    m = 100L, alpha2 = 1, alpha3 = 1,
                    # fun = c("dki1", "dk2"),
                    check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
    ## If A or B is missing, let it be an identity matrix
    if(missing(B)) {
        if(missing(D)) stop("Provide at least one of B and D")
        n <- dim(D)[1L]
        In <- diag(n)
        B <- In
    } else {
        n <- dim(B)[1L]
        In <- diag(n)
    }
    if(missing(D)) D <- In
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
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "alpha3 should be a scalar with 0 < alpha3 < 2" = {
            is.numeric(alpha3) &&
            length(alpha3) == 1 &&
            alpha3 > 0 &&
            alpha3 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
    ## Rotate D and mu with eigenvectors of B
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(D, tol_zero)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    if(use_vec) {
        LD <- diag(D)
        b3 <- alpha3 / max(LD)
        LBh <- rep.int(1, n) - b2 * LB
        LDh <- rep.int(1, n) - b3 * LD
        # if(central) {
        #     dksm <- d2_ij_v(LBh, LDh, m)
        # } else {
        #     dksm <- h2_ij_v(LBh, LDh, mu, m)
        # }
    } else {
        eigD <- eigen(D, symmetric = TRUE)
        LD <- eigD$values
        b3 <- alpha3 / max(LD)
        Bh <- In - b2 * diag(LB)
        Dh <- In - b3 * D
        # if(central) {
        #     dksm <- d2_ij_m(Bh, Dh, m)
        # } else {
        #     dksm <- h2_ij_m(Bh, Dh, mu, m)
        # }
    }
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- IpBDqr_gen_cvE(LB, LD, b2, b3, p, q, r, m)
            } else {
                cppres <- IpBDqr_gen_cmE(LB, D, b2, b3, p, q, r, m)
            }
        } else {
            if(use_vec) {
                cppres <- IpBDqr_gen_nvE(LB, LD, b2, b3, mu, p, q, r, m)
            } else {
                cppres <- IpBDqr_gen_nmE(LB, D, b2, b3, mu, p, q, r, m)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(central) {
            if(use_vec) {
                dksm <- d2_ij_v(LBh, LDh, m)
            } else {
                dksm <- d2_ij_m(Bh, Dh, m)
            }
            ansmat <- hgs_2d(dksm, q, r, n / 2, ((p - q - r) * log(2)
                             + q * log(b2) + r * log(b3) + lgamma(n / 2 + p - q - r)
                             - lgamma(n / 2)))
            ansseq <- sum_counterdiag(ansmat)
        } else {
            if(use_vec) {
                dksm <- h3_ijk_v(rep.int(0, n), LBh, LDh, mu, m)
            } else {
                dksm <- h3_ijk_m(matrix(0, n, n), Bh, Dh, mu, m)
            }
            ansarr <- hgs_3d(dksm, -p, q, r, n / 2, ((p - q - r) * log(2)
                            + q * log(b2) + r * log(b3) + lgamma(n / 2 + p - q - r)
                            - lgamma(n / 2)))
            ansseq <- sum_counterdiag3D(ansarr)
        }
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
    }
## WIP: Condition(s) for existence of moment
    ## Check condition for existence of moment
    rB <- sum(LB > tol_sing)
    ## When A == I, the condition simplifies as A12 == A22 == 0
    cond_exist <- rB / 2 + p > q + r ## condition(1)
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "D should be nonnegative definite" = all(LD >= 0),
              "Moment does not exist in this combination of p, q, r, and rank(B)" = cond_exist)
##
    # ansmat <- hgs_2d(dksm, q, r, n / 2, ((p - q - r) * log(2)
    #                  + q * log(b2) + r * log(b3) + lgamma(n / 2 + p - q - r)
    #                  - lgamma(n / 2)))
    # ansseq <- sum_counterdiag(ansmat)
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
    if(check_convergence && abs(ansseq[length(ansseq)]) > tol_conv) {
        warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
    }
    structure(list(statistic = sum(ansseq), error_bound = NA_real_,
                   res_seq = ansseq, err_seq = NULL), class = "qfrm")
}

## WIP: Condition(s) for existence of moment
#'
#'
#' @export
#'
qfmrm_ApBDqr_int <- function(A, B, D, p = 1, q = 1, r = 1, m = 100L,
                    mu = rep.int(0, n), alpha2 = 1, alpha3 = 1,
                    # fun = c("dki1", "dk2"),
                    check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
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
    }
    if(missing(B)) B <- In
    if(missing(D)) D <- In
    ## Check basic requirements for arguments
    stopifnot(
        "A, B and D should be square matrices" = all(c(dim(A), dim(B), dim(D)) == n),
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
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "alpha3 should be a scalar with 0 < alpha3 < 2" = {
            is.numeric(alpha3) &&
            length(alpha3) == 1 &&
            alpha3 > 0 &&
            alpha3 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
    ## Rotate A, D, and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero) && is_diagonal(D, tol_zero)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    LA <- if(use_vec) diag(A) else eigen(A, symmetric = TRUE)$values
    LD <- if(use_vec) diag(D) else eigen(D, symmetric = TRUE)$values
    b3 <- alpha3 / max(LD)
## WIP: Condition(s) for existence of moment
    ## Check condition for existence of moment
    rB <- sum(LB > tol_sing)
    ## When A == I, the condition simplifies as A12 == A22 == 0
    cond_exist <- rB / 2 + p > q + r ## condition(1)
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "D should be nonnegative definite" = all(LD >= 0),
              "Moment does not exist in this combination of p, q, r, and rank(B)" = cond_exist)
##
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ApBDqr_int_cvE(LA, LB, LD, b2, b3, p, q, r, m)
            } else {
                cppres <- ApBDqr_int_cmE(A, LB, D, b2, b3, p, q, r, m)
            }
        } else {
            if(use_vec) {
                cppres <- ApBDqr_int_nvE(LA, LB, LD, b2, b3, mu, p, q, r, m)
            } else {
                cppres <- ApBDqr_int_nmE(A, LB, D, b2, b3, mu, p, q, r, m)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(use_vec) {
            # LA <- diag(A)
            # LD <- diag(D)
            # b3 <- alpha3 / max(LD)
            LBh <- rep.int(1, n) - b2 * LB
            LDh <- rep.int(1, n) - b3 * LD
            if(central) {
                # dksm <- d3_ijk_v(LA, LBh, LDh, m, m1 = p)[p + 1, , ]
                dksm <- d3_pjk_v(LA, LBh, LDh, m, m1 = p)
            } else {
                # dksm <- htil3_ijk_v(LA, LBh, LDh, mu, m, m1 = p)[p + 1, , ]
                dksm <- htil3_pjk_v(LA, LBh, LDh, mu, m, m1 = p)
            }
        } else {
            # eigD <- eigen(D, symmetric = TRUE)
            # LD <- eigD$values
            # b3 <- alpha3 / max(LD)
            Bh <- In - b2 * diag(LB)
            Dh <- In - b3 * D
            if(central) {
                # dksm <- d3_ijk_m(A, Bh, Dh, m, m1 = p)[p + 1, , ]
                dksm <- d3_pjk_m(A, Bh, Dh, m, m1 = p)
            } else {
                # dksm <- htil3_ijk_m(A, Bh, Dh, mu, m, m1 = p)[p + 1, , ]
                dksm <- htil3_pjk_m(A, Bh, Dh, mu, m, m1 = p)
            }
        }
        dks <- dksm[p + 1, , ]
        ansmat <- hgs_2d(dks, q, r, n / 2 + p, ((p - q - r) * log(2)
                         + q * log(b2) + r * log(b3) + lgamma(p + 1)
                         + lgamma(n / 2 + p - q - r) - lgamma(n / 2 + p)))
        ansseq <- sum_counterdiag(ansmat)
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence && abs(ansseq[length(ansseq)]) > tol_conv) {
        warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
    }
    structure(list(statistic = sum(ansseq), error_bound = NA_real_,
                   res_seq = ansseq, err_seq = NULL), class = "qfrm")
}

## WIP: Condition(s) for existence of moment
#'
#'
#' @export
#'
qfmrm_ApBDqr_npi <- function(A, B, D, p = 1, q = 1, r = 1,
                    m = 100L, mu = rep.int(0, n),
                    alpha1 = 1, alpha2 = 1, alpha3 = 1,
                    # fun = c("dki1", "dk2"),
                    check_convergence = TRUE,
                    use_cpp = FALSE, cpp_method = "Eigen",
                    tol_conv = .Machine$double.eps ^ (1/4),
                    tol_zero = .Machine$double.eps * 100,
                    tol_sing = .Machine$double.eps) {
    if(!missing(cpp_method)) use_cpp <- TRUE
    # cpp_method <- match.arg(cpp_method)
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
    }
    if(missing(B)) B <- In
    if(missing(D)) D <- In
    ## Check basic requirements for arguments
    stopifnot(
        "A, B and D should be square matrices" = all(c(dim(A), dim(B), dim(D)) == n),
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
        "alpha1 should be a scalar with 0 < alpha1 < 2" = {
            is.numeric(alpha1) &&
            length(alpha1) == 1 &&
            alpha1 > 0 &&
            alpha1 < 2
        },
        "alpha2 should be a scalar with 0 < alpha2 < 2" = {
            is.numeric(alpha2) &&
            length(alpha2) == 1 &&
            alpha2 > 0 &&
            alpha2 < 2
        },
        "alpha3 should be a scalar with 0 < alpha3 < 2" = {
            is.numeric(alpha3) &&
            length(alpha3) == 1 &&
            alpha3 > 0 &&
            alpha3 < 2
        },
        "mu should be an n-vector" = length(mu) == n
    )
    if((p %% 1) == 0) {
        warning("For integral p, qfmrm_ApBDqr_int() works better")
    }
    eigB <- eigen(B, symmetric = TRUE)
    LB <- eigB$values
    b2 <- alpha2 / max(LB)
    ## Rotate A, D, and mu with eigenvectors of B
    A <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    D <- with(eigB, crossprod(crossprod(D, vectors), vectors))
    mu <- c(crossprod(eigB$vectors, c(mu)))
    use_vec <- is_diagonal(A, tol_zero) && is_diagonal(D, tol_zero)
    central <- iseq(mu, rep.int(0, n), tol_zero)
    LA <- if(use_vec) diag(A) else eigen(A, symmetric = TRUE)$values
    LD <- if(use_vec) diag(D) else eigen(D, symmetric = TRUE)$values
    b1 <- alpha1 / max(abs(LA))
    b3 <- alpha3 / max(LD)
## WIP: Condition(s) for existence of moment
    ## Check condition for existence of moment
    rB <- sum(LB > tol_sing)
    ## When A == I, the condition simplifies as A12 == A22 == 0
    cond_exist <- rB / 2 + p > q + r ## condition(1)
    stopifnot("B should be nonnegative definite" = all(LB >= 0),
              "D should be nonnegative definite" = all(LD >= 0),
              "Moment does not exist in this combination of p, q, r, and rank(B)" = cond_exist)
##
    if(use_cpp) {
        if(central) {
            if(use_vec) {
                cppres <- ApBDqr_npi_cvE(LA, LB, LD, b1, b2, b3, p, q, r, m)
            } else {
                cppres <- ApBDqr_npi_cmE(A, LB, D, b1, b2, b3, p, q, r, m)
            }
        } else {
            if(use_vec) {
                cppres <- ApBDqr_npi_nvE(LA, LB, LD, b1, b2, b3, mu, p, q, r, m)
            } else {
                cppres <- ApBDqr_npi_nmE(A, LB, D, b1, b2, b3, mu, p, q, r, m)
            }
        }
        ansseq <- cppres$ansseq
    } else {
        if(use_vec) {
            # LA <- diag(A)
            # LD <- diag(D)
            # b1 <- alpha1 / max(abs(LA))
            # b3 <- alpha3 / max(LD)
            LAh <- rep.int(1, n) - b1 * LA
            LBh <- rep.int(1, n) - b2 * LB
            LDh <- rep.int(1, n) - b3 * LD
            if(central) {
                dksm <- d3_ijk_v(LAh, LBh, LDh, m)
            } else {
                dksm <- h3_ijk_v(LAh, LBh, LDh, mu, m)
            }
        } else {
            # eigA <- eigen(A, symmetric = TRUE)
            # eigD <- eigen(D, symmetric = TRUE)
            # LA <- eigA$values
            # LD <- eigD$values
            # b1 <- alpha1 / max(abs(LA))
            # b3 <- alpha3 / max(LD)
            Ah <- In - b1 * A
            Bh <- In - b2 * diag(LB)
            Dh <- In - b3 * D
            if(central) {
                dksm <- d3_ijk_m(Ah, Bh, Dh, m)
            } else {
                dksm <- h3_ijk_m(Ah, Bh, Dh, mu, m)
            }
        }
        ansarr <- hgs_3d(dksm, -p, q, r, n / 2, ((p - q - r) * log(2)
                         - p * log(b1) + q * log(b2) + r * log(b3)
                         + lgamma(n / 2 + p - q - r) - lgamma(n / 2)))
        ansseq <- sum_counterdiag3D(ansarr)
        scf <- attr(dksm, "scale")
        ansseq <- ansseq / scf
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
    if(check_convergence && abs(ansseq[length(ansseq)]) > tol_conv) {
        warning("Last term of the series is larger than specified tolerance, ",
            tol_conv, "\n  Check sensitivity with different m?")
    }
    structure(list(statistic = sum(ansseq), error_bound = NA_real_,
                   res_seq = ansseq, err_seq = NULL), class = "qfrm")
}
