
S_fromUL <- function(evec = evec, evalues = evalues) {
    te <- t(evec)
    crossprod(te * evalues, te)
}

# matsqrt_K <- function(S, tol = .Machine$double.eps * 100) {
#     if(!isSymmetric(S)) stop("Covariance matrix should be symmetric")
#     svdS <- svd(S, nv = 0)
#     d <- svdS$d
#     u <- svdS$u
#     if(any(d < 0)) stop("Covariance matrix should be nonnegative definite")
#     pos <- d > tol
#     u[, pos] %*% diag(sqrt(d[pos]))
#     # t(t(u[, pos]) * sqrt(d[pos]))
# }

KiK <- function(S, tol = .Machine$double.eps * 100) {
    if(!isSymmetric(S)) stop("Covariance matrix should be symmetric")
    svdS <- svd(S, nv = 0)
    d <- svdS$d
    u <- svdS$u
    if(any(d < 0)) stop("Covariance matrix should be nonnegative definite")
    pos <- d > tol
    K <- u[, pos] %*% diag(sqrt(d[pos]))
    iK <- diag(1 / sqrt(d[pos])) %*% t(u[, pos])
    list(K = K, iK = iK)
}

tr <- function(X) sum(diag(X))

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

iseq <- function(x, y = rep.int(0, length(x)),
                 tol = .Machine$double.eps * 100) {
    isTRUE(all.equal(c(x), c(y), tol = tol, check.attributes = FALSE))
}

is_diagonal <- function(A, tol = .Machine$double.eps * 100) {
    n <- dim(A)[1]
    In <- diag(n)
    isTRUE(all.equal(as.numeric(abs(A) > tol), c(In)))
}
