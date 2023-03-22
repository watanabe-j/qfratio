test_that("Expect average eigenvalues when p = q = 1", {
    nvs <- 2:10
    for(nv in nvs) {
        L1 <- 1:nv
        A1 <- diag(L1)

        trA1 <- tr(A1)
        A1I1 <- qfrm(A1, p = 1)$statistic

        expect_equal(trA1 / nv, A1I1)
    }
})

test_that("Expect ordinary positive moments when q = 0", {
    nvs <- 2:10
    ks <- 1:5
    for(nv in nvs) {
        L1 <- 1:nv
        A1 <- diag(L1)
        for(k in ks) {
            Ak <- qfm_Ap_int(A1, p = k)$statistic
            AkI0 <- qfrm(A1, p = k, q = 0)$statistic

            expect_equal(Ak, AkI0)
        }
    }
})

if(requireNamespace("stats", quietly = TRUE)) {
    test_that("Expect identical results for simultaneously rotated matrices", {
        nvs <- 2:10
        ps <- 1:3
        for(nv in nvs) {
            L1 <- 1:nv
            L2 <- nv:1
            L3 <- sqrt(nv:1)
            A1 <- diag(L1)
            A2 <- diag(L2)
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1r <- Q %*% A1 %*% t(Q)
            A2r <- Q %*% A2 %*% t(Q)

            for(p in ps) {
                expect_equal(qfrm(A1, A2, p), qfrm(A1r, A2r, p))
            }
        }
    })
}

# This yields a message once per session, which is to be ignored
suppressMessages(
    qfrm(diag(4), p = 1/2, mu = rep.int(1, 4))
)

test_that("Expect silence or warning around error bound", {
    nvs <- 4:6
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv / nv
        A1s <- diag(L1 - 1)

        ## Acceptable parameter values
        expect_silent(qfrm(A1, p = 1/2, alphaA = 0.8,           m = m, check_convergence = FALSE))
        expect_silent(qfrm(A1, p = 1/2, mu = mu, alphaA = 0.8,  m = m, check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, alphaB = 0.8,                m = m, check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, mu = mu, alphaB = 0.8,       m = m, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, alphaB = 0.8,               m = m, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, mu = mu, alphaB = 0.8,      m = m, check_convergence = FALSE))

        ## Unacceptable value / singular argument causing warning in calculating error bound
        expect_warning(qfrm(A1, p = 1/2, alphaA = 1.5,      m = m, check_convergence = FALSE))
        expect_warning(qfrm(A1s, p = 1/2,                   m = m, check_convergence = FALSE))
        expect_warning(qfrm(A1, A2, alphaB = 1.5,           m = m, check_convergence = FALSE))
        expect_warning(qfrm(A1, A2, mu = mu, alphaB = 1.5,  m = m, check_convergence = FALSE))
        expect_warning(qfrm(A1, A1s,                        m = m, check_convergence = FALSE))
        expect_warning(qfrm(A1, A1s, mu = mu,               m = m, check_convergence = FALSE))
        expect_warning(qfmrm(A1, A1s,                       m = m, check_convergence = FALSE))
        expect_warning(qfmrm(A1, A1s, mu = mu,              m = m, check_convergence = FALSE))
        expect_warning(qfmrm(A1, A2, alphaB = 1.5,          m = m, check_convergence = FALSE))
        expect_warning(qfmrm(A1, A2, mu = mu, alphaB = 1.5, m = m, check_convergence = FALSE))

        ## Above is dismissed when no error bound is returned
        expect_silent(qfrm(A1, p = 1/2, alphaA = 1.5,      m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1s, p = 1/2,                   m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, alphaB = 1.5,           m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, mu = mu, alphaB = 1.5,  m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A1s,                        m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A1s, mu = mu,               m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A1s,                       m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A1s, mu = mu,              m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, alphaB = 1.5,          m = m, error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, mu = mu, alphaB = 1.5, m = m, error_bound = FALSE, check_convergence = FALSE))

        ## Does not matter when error bound is unavailable
        expect_silent(qfrm(A1, p = 1/2, mu = mu, alphaA = 1.5,  m = m, check_convergence = FALSE))
        expect_silent(qfrm(A1s, p = 1/2, mu = mu,               m = m, check_convergence = FALSE))

    }
})


test_that("Existence conditions: qfrm, nonsingular", {
    nvs <- 2:4
    ks <- c(1:3, 1/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        I <- diag(nv)
        mu <- 1:nv / nv
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
        }

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks) {
                if(nv / 2 + p <= q) {
                    # expect_error(qfrm(A1, I,  p, q, m = m))
                    # expect_error(suppressWarnings(qfrm(A1, I,  p, q, m = m, mu = mu)))
                    expect_error(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, I,  p, q, m = m))
                    if(requireNamespace("gsl", quietly = TRUE)) expect_silent(qfrm(A1, I,  p, q, m = m, mu = mu))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks) {
                if(nv / 2 + p <= q) {
                    # expect_error(qfrm(A1, I,  p, q, m = m, check_convergence = FALSE))
                    # expect_error(suppressWarnings(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE)))
                    expect_error(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, I,  p, q, m = m, check_convergence = FALSE))
                    expect_silent(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                }
            }
        }
    }
})

test_that("Existence conditions: qfmrm, nonsingular", {
    nvs <- 2:4
    ks <- c(1, 2, 1/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        I <- diag(nv)
        mu <- 1:nv / nv
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
            A3 <- Q %*% A3 %*% t(Q)
        }

        for(p in ks) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(nv / 2 + p <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
    }
})

test_that("Existence conditions: qfrm, singular A12 = A22 = 0", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L1[nv] <- 0
        L2[nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks) {
                if((nv - 1) / 2 + p <= q) {
                    expect_error(qfrm(A1, A2,  p, q, m = m,          error_bound = FALSE, check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, error_bound = FALSE, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, A2,  p, q, m = m,          error_bound = FALSE, check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, error_bound = FALSE, check_convergence = FALSE))
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks) {
                if((nv - 1) / 2 + p <= q) {
                    expect_error(qfrm(A1, A2,  p, q, m = m,          check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, A2,  p, q, m = m,          check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                }
            }
        }
    }
})

test_that("Existence conditions: qfrm, singular, A22 = 0, A12 != 0", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L1[nv] <- 0
        L2[nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        A1[nv, 1] <- A1[1, nv] <- 0.5
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks) {
                if(((nv - 1) + p) / 2 <= q) {
                    expect_error(qfrm(A1, A2,  p, q, m = m,          error_bound = FALSE, check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, error_bound = FALSE, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, A2,  p, q, m = m,          error_bound = FALSE, check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, error_bound = FALSE, check_convergence = FALSE))
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks) {
                if(((nv - 1) + p) / 2 <= q) {
                    expect_error(qfrm(A1, A2,  p, q, m = m,          check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                } else {
                    # ## Here A1 is indefinite and moment is undefined
                    # expect_silent(qfrm(A1, A2,  p, q, m = m,          check_convergence = FALSE))
                    # expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                }
            }
        }
    }
})

test_that("Existence conditions: qfrm, singular, A22 != 0", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L1[nv] <- 0
        L2[nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        A1[nv, nv] <- A1[nv, 1] <- A1[1, nv] <- 0.5
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks) {
                if((nv - 1) / 2 <= q) {
                    expect_error(qfrm(A1, A2,  p, q, m = m,          error_bound = FALSE, check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, error_bound = FALSE, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, A2,  p, q, m = m,          error_bound = FALSE, check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, error_bound = FALSE, check_convergence = FALSE))
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks) {
                if((nv - 1) / 2 <= q) {
                    expect_error(qfrm(A1, A2,  p, q, m = m,          check_convergence = FALSE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                } else {
                    expect_silent(qfrm(A1, A2,  p, q, m = m,          check_convergence = FALSE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE))
                }
            }
        }
    }
})


test_that("Existence conditions: qfmrm, singular, range identical, A12 = A22 = 0", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        L1[nv] <- 0
        L2[nv] <- 0
        L3[nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        I <- diag(nv)
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
            A3 <- Q %*% A3 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if((nv - 1) / 2 + p <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m,           error_bound = FALSE, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  error_bound = FALSE, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m,           error_bound = FALSE, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  error_bound = FALSE, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if((nv - 1) / 2 + p <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
    }
})

test_that("Existence conditions: qfmrm, singular, range identical, A22 = 0, A12 != 0", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 1
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        L1[nv] <- 0
        L2[nv] <- 0
        L3[nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        I <- diag(nv)
        A1[nv, 1] <- A1[1, nv] <- 0.5
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
            A3 <- Q %*% A3 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(((nv - 1) + p) / 2 <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m,           error_bound = FALSE, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  error_bound = FALSE, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m,           error_bound = FALSE, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  error_bound = FALSE, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(((nv - 1) + p) / 2 <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        # ## Here A1 is indefinite and moment is undefined
                        # expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m,           check_convergence = FALSE))
                        # expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        # expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        # expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
    }
})

test_that("Existence conditions: qfmrm, singular, A22 != 0", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 1

    ## A1 has nonzero elements in the null space of A2/A3, whose rank is (nv - 1)
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        L2[nv] <- 0
        L3[nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        I <- diag(nv)
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
            A3 <- Q %*% A3 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if((nv - 1) / 2 <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m,           error_bound = FALSE, check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  error_bound = FALSE, check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m,           error_bound = FALSE, check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  error_bound = FALSE, check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if((nv - 1) / 2 <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
    }

    ## A1 has nonzero elements in the null space of A3, whose rank is (nv - 2)
    for(nv in nvs[nvs != 2]) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        L1[nv] <- 0
        L2[nv] <- 0
        L3[(nv - 1):nv] <- 0
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        I <- diag(nv)
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A1 <- Q %*% A1 %*% t(Q)
            A2 <- Q %*% A2 %*% t(Q)
            A3 <- Q %*% A3 %*% t(Q)
        }
        mu <- 1:nv / nv

        for(p in ks[(ks %% 1) == 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if((nv - 2) / 2 <= q + r) {
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
        for(p in ks[(ks %% 1) != 0]) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if((nv - 2) / 2 <= q + r) {
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_error(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_error(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    } else {
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m,           check_convergence = FALSE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu,  check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m,          check_convergence = FALSE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE))
                    }
                }
            }
        }
    }
})
