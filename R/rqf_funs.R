#' @export
#'
rqfr <- function(nit = 1000L, A, B, p = 1, q = p, mu, Sigma, use_cpp = FALSE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package \"mvtnorm\" is required to use this function.")
    }
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
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    # X <- eigvaldisp::rmvn(nit, mean = mu, Sigma = Sigma)
    if(use_cpp) {
        return(rqfrE(nit, A, B, p, q, mu, Sigma))
    }
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    diag(tcrossprod(tcrossprod(X, A), X)) ^ p / diag(tcrossprod(tcrossprod(X, B), X)) ^ q
}

#' @export
#'
rqfmr <- function(nit = 1000L, A, B, D, p = 1, q = p, r = q, mu, Sigma, use_cpp = FALSE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package \"mvtnorm\" is required to use this function.")
    }
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
    if(missing(p) && !missing(q)) p <- q + r
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    if(use_cpp) {
        return(rqfmrE(nit, A, B, D, p, q, r, mu, Sigma))
    }
    # X <- eigvaldisp::rmvn(nit, mean = mu, Sigma = Sigma)
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    diag(tcrossprod(tcrossprod(X, A), X)) ^ p / diag(tcrossprod(tcrossprod(X, B), X)) ^ q / diag(tcrossprod(tcrossprod(X, D), X)) ^ r
}

#' @export
#'
rqfp <- function(nit = 1000L, A, B, D, p = 1, q = p, r = q, mu, Sigma, use_cpp = FALSE) {
    if(!requireNamespace("mvtnorm", quietly = TRUE) && !use_cpp) {
        stop("Package \"mvtnorm\" is required to use this function.")
    }
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
    if(missing(p) && !missing(q)) p <- q + r
    if(missing(mu)) mu <- rep.int(0, n)
    if(missing(Sigma)) Sigma <- In
    if(use_cpp) {
        return(rqfpE(nit, A, B, D, p, q, r, mu, Sigma))
    }
    # X <- eigvaldisp::rmvn(nit, mean = mu, Sigma = Sigma)
    X <- mvtnorm::rmvnorm(nit, mean = mu, sigma = Sigma)
    diag(tcrossprod(tcrossprod(X, A), X)) ^ p * diag(tcrossprod(tcrossprod(X, B), X)) ^ q * diag(tcrossprod(tcrossprod(X, D), X)) ^ r
}
