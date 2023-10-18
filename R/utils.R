##### S_fromUL #####
#' Make covariance matrix from eigenstructure
#'
#' This is an internal utility function to make covariance matrix from
#' eigenvectors and eigenvalues.  Symmetry is assumed for the original matrix.
#'
#' @param evec
#'   Matrix whose columns are eigenvectors
#' @param evalues
#'   Vector of eigenvalues
#'
S_fromUL <- function(evec, evalues) {
    te <- t(evec)
    crossprod(te * evalues, te)
}

##### KiK #####
#' Matrix square root and generalized inverse
#'
#' This internal function calculates the decomposition
#' \eqn{\mathbf{S} = \mathbf{K} \mathbf{K}^T}{S = K K^T} for an
#' \eqn{n \times n}{n x n} covariance matrix \eqn{\mathbf{S}}{S}, so that
#' \eqn{\mathbf{K}}{K} is an \eqn{n \times m}{n x m} matrix with \eqn{m} being
#' the rank of \eqn{\mathbf{S}}{S}.  Returns this
#' \eqn{\mathbf{K}}{K} and its generalized inverse,
#' \eqn{\mathbf{K}^-}{K^-}, in a list.
#'
#' At present, this utilizes \code{svd()},
#' although there may be better alternatives.
#'
#' @param S
#'   Covariance matrix.  Symmetry and positive (semi-)definiteness are checked.
#' @param tol
#'   Tolerance to determine the rank of \eqn{\mathbf{S}}{S}.  Eigenvalues
#'   smaller than this value are considered zero.
#'
#' @return
#'   List with \code{K} and \code{iK}, with the latter being
#'   \eqn{\mathbf{K}^-}{K^-}
#'
KiK <- function(S, tol = .Machine$double.eps * 100) {
    if(!isSymmetric(S)) stop("Covariance matrix must be symmetric")
    svdS <- svd(S, nv = 0)
    d <- svdS$d
    u <- svdS$u
    if(any(d < 0)) stop("Covariance matrix must be nonnegative definite")
    pos <- d > tol
    K <- u[, pos] %*% diag(sqrt(d[pos]), nrow = sum(pos))
    iK <- diag(1 / sqrt(d[pos]), nrow = sum(pos)) %*% t(u[, pos])
    list(K = K, iK = iK)
}

##### tr #####
#' Matrix trace function
#'
#' This is an internal function.  No check is done on the structure of \code{X}.
#'
#' @param X
#'   Square matrix whose trace is to be calculated
#'
tr <- function(X) sum(diag(X))

##### sum_counterdiag #####
#' Summing up counter-diagonal elements
#'
#' \code{sum_counterdiag()} sums up counter-diagonal elements of
#' a square matrix from the upper-left; i.e.,
#' \code{c(X[1, 1], X[1, 2] + X[2, 1], X[1, 3] + X[2, 2] + X[3, 1], ...)}
#' (or a sequence of \eqn{\sum_{i=1}^k x_{i, k - i + 1}} for
#' \eqn{k = 1, 2, ...}).  \code{sum_counterdiag3D()} does a comparable in
#' a 3D cubic array.  No check is done on the structure of \code{X}.
#'
#' @param X
#'   Square matrix or cubic array
#'
#' @name sum_counterdiag
#'
sum_counterdiag <- function(X) {
    n <- nrow(X)
    ans <- rep.int(0, n)
    for(i in 1:n) {
        for(j in 1:i) {
            x <- X[i - j + 1, j]
            if(!is.na(x)) ans[i] <- ans[i] + x
        }
    }
    return(ans)
}

#' @rdname sum_counterdiag
sum_counterdiag3D <- function(X) {
    n <- dim(X)[1]
    ans <- rep.int(0, n)
    for(i in 1:n) {
        for(j in 1:i) {
            for(k in 1:(i - j + 1)) {
                x <- X[i - j - k + 2, j, k]
                if(!is.na(x)) ans[i] <- ans[i] + x
            }
        }
    }
    return(ans)
}

##### iseq #####
#' Are these vectors equal?
#'
#' This internal function is used to determine whether two vectors/matrices have
#' the same elements (or, a vector/matrix is all equal to 0)
#' using \code{all.equal()}.  Attributes and dimensions are ignored as
#' they are passed as vectors using \code{c()}.
#'
#' @param x
#'   Main \code{target} vector/matrix in \code{all.equal()}
#' @param y
#'   \code{current} in \code{all.equal()}.  Default zero vector.
#' @param tol
#'   Numeric to specify \code{tolerance} in \code{all.equal()}
#'
#' @seealso \code{\link[base]{all.equal}}
#'
iseq <- function(x, y = rep.int(0, length(x)),
                 tol = .Machine$double.eps * 100) {
    isTRUE(all.equal(c(x), c(y), tol = tol, check.attributes = FALSE))
}

##### is_diagonal #####
#' Is this matrix diagonal?
#'
#' This internal function is used to determine whether a square matrix
#' is diagonal (within a specified tolerance).  Returns \code{TRUE}
#' when the absolute values of all off-diagonal elements
#' are below \code{tol}, using \code{all.equal()}.
#'
#' @inheritParams iseq
#'
#' @param A
#'   Square matrix.  No check is done.
#' @param symmetric
#'   If \code{FALSE} (default), sum of absolute values of the corresponding
#'   lower and upper triangular elements are examined with a doubled
#'   \code{tol}.  If \code{TRUE}, only the lower triangular elements are
#'   examined assuming symmetry.
#'
#' @seealso \code{\link[base]{all.equal}}
#'
is_diagonal <- function(A, tol = .Machine$double.eps * 100, symmetric = FALSE) {
    n <- dim(A)[1]
    if(symmetric) {
        return(isTRUE(all.equal(A[lower.tri(A)],
                                rep.int(0, n * (n - 1) / 2), tol)))
    } else {
        A <- abs(A)
        return(isTRUE(all.equal((A + t(A))[lower.tri(A)],
                                rep.int(0, n * (n - 1) / 2), tol * 2)))
    }
}

##### range_qfr #####
#' Get range of ratio of quadratic forms
#'
#' \code{range_qfr()}: internal function to obtain the possible range of
#' a ratio of quadratic forms,
#' \eqn{\frac{ \mathbf{x^{\mathit{T}} A x} }{ \mathbf{x^{\mathit{T}} B x} }
#' }{ (x^T A x) / (x^T B x) }.
#'
#' @param A,B
#'   Square matrices.  No check is done.
#' @param eigB
#'   Result of \code{eigen(B)} can be passed when already computed
#' @param tol
#'   Tolerance to determine numerical zero
#'
range_qfr <- function(A, B, eigB = eigen(B, symmetric = TRUE),
                      tol = .Machine$double.eps * 100, t = 1e-3) {
    LB <- eigB$values
    rB <- sum(LB > tol)
    n <- length(LB)
    Ad <- with(eigB, crossprod(crossprod(A, vectors), vectors))
    if(rB == n) {
        LBiA <- gen_eig(A, B, eigB, Ad, tol = tol, t = t)
    } else {
        A11 <- Ad[1:rB, 1:rB]
        A12 <- Ad[1:rB, (rB + 1):n]
        A22 <- Ad[(rB + 1):n, (rB + 1):n]
        A12_is_zero <- all(abs(A12) < tol)
        A22_is_zero <- all(abs(A22) < tol)
        ## If A12 and A22 are zero, B's null space is within A's
        ## It suffices to consider the non-null space of B
        if(A12_is_zero && A22_is_zero) {
            LBiA <- gen_eig(A11, diag(LB[1:rB]), tol = tol, t = t)
        } else {
            ## This is not ideal but try anyway
            LBiA <- try(gen_eig(A, B, eigB, Ad, tol = tol, t = t), TRUE)
            if(inherits(LBiA, "try-error")) return(c(-Inf, Inf))
        }
    }
    ## NaN in LBiA usually corresponds to common null space so is negligible
    res <- range(LBiA[!is.nan(LBiA)])
    return(res)
}

##### gen_eig #####
#' Get generalized eigenvalues
#'
#' \code{gen_eig()} is an internal function to obtain generalized eigenvalues,
#' i.e., roots of
#' \eqn{\det{\mathbf{A} - \lambda \mathbf{B}} = 0}{A - \lambda B = 0},
#' which are the eigenvalues of \eqn{\mathbf{B}^{-1} \mathbf{A}} if
#' \eqn{\mathbf{B}} is nonsingular.
#'
#' \code{gen_eig()} solves the generalized eigenvalue problem with
#' Jennings et al.'s (1978) algorithm.  The sign of infinite eigenvalue
#' (when present) cannot be determined from this algorithm, so is deduced
#' as follows: (1) \eqn{\mathbf{A}}{A} and \eqn{\mathbf{B}}{B} are rotated by
#' the eigenvectors of \eqn{\mathbf{B}}{B}; (2) the submatrix of rotated
#' \eqn{\mathbf{A}}{A} corresponding to the null space of \eqn{\mathbf{B}}{B}
#' is examined; (3) if this is nonnegative (nonpositive) definite, the result
#' must have positive (negative, resp.) infinity; if this is indefinite,
#' the result must have both positive and negative infinities;
#' if this is (numerically) zero, the result must have \code{NaN}.  The last
#' case is expeted to happen very rarely, as in this case Jennings algorithm
#' would fail.  This is where the null space of \eqn{\mathbf{B}}{B} is
#' a subspace of that of \eqn{\mathbf{A}}{A}, so that the range of ratio of
#' quadratic forms can be well-behaved.  \code{range_qfr()} tries to detect
#' this case and handle the range accordingly, but if that is infeasible
#' it returns \code{c(-Inf, Inf)}.
#'
#' @param Ad
#'   \code{A} rotated with eigenvectors of \code{B} can be passed
#'   when already computed
#' @param t
#'   Tolerance used to determine whether estimates are numerically stable;
#'   \eqn{t} in Jennings et al. (1978).
#'
#' @references
#' Jennings, A., Halliday, J. and Cole, M. J. (1978) Solution of linear
#'   generalized eigenvalue problems containing singular matrices.
#'   *Journal of the Institute of Mathematics and Its Applications*, **22**,
#'   401--410.
#'   \doi{10.1093/imamat/22.4.401}.
#'
#' @rdname range_qfr
#'
gen_eig <- function(A, B, eigB = eigen(B, symmetric = TRUE),
                    Ad = with(eigB, crossprod(crossprod(A, vectors), vectors)),
                    tol = .Machine$double.eps * 100, t = 1e-3) {
    LB <- eigB$values
    rB <- sum(LB > tol)
    n <- length(LB)
    if(rB == n) {
        BiA <- t(Ad / sqrt(LB)) / sqrt(LB)
        LBiA <- eigen(BiA, symmetric = TRUE, only.values = TRUE)$values
    } else {
        ## Use Jennings algorithm
        sqnorm_B <- sum(B ^ 2)
        alpha <- sum(A * B) / sqnorm_B
        Abar <- A - alpha * B
        i <- 1
        fail_prev <- FALSE
        while(TRUE) {
            LM <- try(eigen(solve(Abar, B), only.values = TRUE)$values, TRUE)
            fail_curr <- inherits(LM, "try-error")
            if(fail_prev && fail_curr) stop("problem looks ill-conditioned")
            if(!fail_curr) {
                if(all(abs(LM - alpha) > t * abs(alpha))) break
            }
            ## If above didn't work, failure; update alpha
            fail_prev <- TRUE
            alpha <- alpha + sqrt(sum(A ^ 2) / sqnorm_B)
            Abar <- A - alpha * B
            ## To avoid infinite loop; this should not happen
            i <- i + 1
            if(i > 5) stop("unexpected error: max iteration reached; ",
                           "contact maintainer")
        }
        LBiA <- 1 / LM + alpha
        ## If any of LM is ~0, the corresponding LBiA should be +/-Inf
        ## but its sign from LM can be inaccurate.
        ## Determine the sign(s) by looking at the submatrix of Ad
        ## corresponding to the null space of B
        inf_inds <- abs(LM) < tol
        if(any(inf_inds)) {
            A22 <- Ad[(rB + 1):n, (rB + 1):n]
            LA22 <- eigen(A22, symmetric = TRUE, only.values = TRUE)$values
            LA22min <- min(LA22)
            LA22max <- max(LA22)
            ## If A22 is indefinite, LBiA should have both +Inf and -Inf
            ##           nonnegative,                 +Inf only
            ##           nonpositive,                 -Inf only
            ##           not clearly any of the above, use NaN (should be rare)
            if(LA22min < -tol && LA22max > tol) {
                LBiA[inf_inds] <- rep_len(c(Inf, -Inf), sum(inf_inds))
            } else if(LA22min > -tol && LA22max > tol) {
                LBiA[inf_inds] <-  Inf
            } else if(LA22min < -tol && LA22max < tol) {
                LBiA[inf_inds] <- -Inf
            } else {
                LBiA[inf_inds] <- NaN
            }
        }
    }
    ## Mathematically, the generalized eigenvalues should be real
    ## But LBiA or LM may be complex for numerically infinite ones,
    ## so Re() is taken for safeguarding
    return(Re(LBiA))
}
