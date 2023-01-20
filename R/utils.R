##### S_fromUL #####
#' Make covariance matrix from eigenstructure
#'
#' This is an internal utility function to make covariance matrix from
#' eigenvectors and eigenvalues. Symmetry is assumed for the original matrix.
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
#' \eqn{\mathbf{S} = \mathbf{K} \mathbf{K}^T} for an \eqn{n \times n}
#' covariance matrix \eqn{\mathbf{S}}, so that \eqn{\mathbf{K}} is an
#' \eqn{n \times m} matrix with \eqn{m} being the rank of \eqn{\mathbf{S}}.
#' Returns this \eqn{\mathbf{K}} and its generalized inverse,
#' \eqn{\mathbf{K}^-}, in a list.
#'
#' At present, this utilizes \code{svd()},
#' although there may be better alternatives.
#'
#' @param S
#'   Covariance matrix. Symmetry and positive (semi-)definiteness are checked.
#' @param tol
#'   Tolerance to determine the rank of \eqn{\mathbf{S}}.
#'   Eigenvalues smaller than this value is considered zero.
#'
#' @return
#'   List with \code{K} and \code{iK}, with the latter being \eqn{\mathbf{K}^-}
#'
#' @examples
#' (S1 <- diag(c(3, 2, 1, 0)))
#' qfratio:::KiK(S1)
#'
#' (S2 <- diag(1e-4^(1:5)))
#' qfratio:::KiK(S2)
#' ## Note that the smallest(s) of the eigenvalues was considered zero above
#' ## Specify tol if this is to be avoided
#' qfratio:::KiK(S2, tol = 1e-20)
#'
KiK <- function(S, tol = .Machine$double.eps * 100) {
    if(!isSymmetric(S)) stop("Covariance matrix should be symmetric")
    svdS <- svd(S, nv = 0)
    d <- svdS$d
    u <- svdS$u
    if(any(d < 0)) stop("Covariance matrix should be nonnegative definite")
    pos <- d > tol
    K <- u[, pos] %*% diag(sqrt(d[pos]), nrow = sum(pos))
    iK <- diag(1 / sqrt(d[pos]), nrow = sum(pos)) %*% t(u[, pos])
    list(K = K, iK = iK)
}

##### tr #####
#' Matrix trace function
#'
#' This is an internal function. No check is done on the structure of \code{X}.
#'
#' @param X
#'   Square matrix whose trace is to be calculated
#'
tr <- function(X) sum(diag(X))

##### sum_counterdiag #####
#' Summing up counter-diagonal elements
#'
#' sum_counterdiag() sums up counter-diagonal elements of a square matrix from
#' the upper-left; i.e.,
#' \code{c(X[1, 1], X[1, 2] + X[2, 1], X[1, 3] + X[2, 2] + X[3, 1], ...)}
#' (or a sequence of
#' \eqn{\sum_{i=1}^k x_{i, k - i + 1}} for \eqn{k = 1, 2, ...}).
#' sum_counterdiag3D() does a comparable in a 3D cubic array.
#' No check is done on the structure of \code{X}.
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
#' using \code{all.equal()}. Attributes and dimensions are ignored as
#' they are passed as vectors using \code{c()}.
#'
#' @param x
#'   Main \code{target} vector/matrix in \code{all.equal()}
#' @param y
#'   \code{current} in \code{all.equal()}. Default zero vector.
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
#' is diagonal (within a specified tolerance).
#' Returns \code{TRUE} when the absolute values of all off-diagonal elements
#' are below \code{tol}, using \code{all.equal()}.
#'
#' @inheritParams iseq
#'
#' @param A
#'   Square matrix. No check is done.
#' @param symmetric
#'   If \code{FALSE} (default), sum of absolute values of the corresponding
#'   lower and upper triangular elements are examined with a doubled \code{tol}.
#'   If \code{TRUE}, only the lower triangular elements are examined
#'   assuming symmetry.
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
