test_that("Expect average eigenvalues when nv = q = 1", {
    nvs <- 2:10
    # ks <- 1:5
    for(nv in nvs) {
        L1 <- 1:nv
        # L2 <- nv:1
        # L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        # A2 <- diag(L2)
        # A3 <- diag(L3)
        # mu <- 1:nv / nv

        trA1 <- tr(A1)
        A1I1 <- qfrm(A1, nv = 1)$statistic

        expect_equal(trA1 / nv, A1I1)
    }
})

test_that("Expect ordinary positive moments when q = 0", {
    nvs <- 2:10
    ks <- 1:5
    for(nv in nvs) {
        L1 <- 1:nv
        # L2 <- nv:1
        # L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        # A2 <- diag(L2)
        # A3 <- diag(L3)
        # mu <- 1:nv / nv
        for(k in ks) {
            Ak <- qfm_Ap_int(A1, p = k)$statistic
            AkI0 <- qfrm(A1, p = k, q = 0)$statistic
            AkII00 <- qfmrm(A1, p = k, q = 0, r = 0)$statistic

            expect_equal(Ak, AkI0)
            expect_equal(Ak, AkII00)
        }
    }
})

test_that("Expect identical results for simultaneously rotated matrices", {
    nvs <- 2:10
    ps <- 1:3
    # ks <- 1:5
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        # A3 <- diag(L3)
        # mu <- 1:nv / nv
        Q <- qr.Q(qr(matrix(rnorm(nv^2), nv, nv)))
        A1r <- Q %*% A1 %*% t(Q)
        A2r <- Q %*% A2 %*% t(Q)
        # A3r <- Q %*% A3 %*% t(Q)

        for(p in ps) {
            expect_equal(qfrm(A1, A2, p), qfrm(A1r, A2r, p))
        }
    }
})

# This yields a warning once per session
qfrm(diag(4), p = 1/2, mu = rep.int(1, 4), check_convergence = FALSE)

test_that("Expect silence or warning around error bound", {
    nvs <- 4:6
    # ps <- 1:3
    # ks <- 1:5
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        # L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        # A3 <- diag(L3)
        mu <- 1:nv / nv
        # Q <- qr.Q(qr(matrix(rnorm(nv^2), nv, nv)))
        # A1r <- Q %*% A1 %*% t(Q)
        # A2r <- Q %*% A2 %*% t(Q)
        # A3r <- Q %*% A3 %*% t(Q)
        A1s <- diag(L1 - 1)

        ## Acceptable parameter values
        expect_silent(qfrm(A1, p = 1/2, alphaA = 0.8,           check_convergence = FALSE))
        expect_silent(qfrm(A1, p = 1/2, mu = mu, alphaA = 0.8,  check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, alphaB = 0.8,                check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, mu = mu, alphaB = 0.8,       check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, alphaB = 0.8,               check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, mu = mu, alphaB = 0.8,      check_convergence = FALSE))

        ## Unacceptable value / singular argument causing warning in calculating error bound
        expect_warning(qfrm(A1, p = 1/2, alphaA = 1.5,      check_convergence = FALSE))
        expect_warning(qfrm(A1s, p = 1/2,                   check_convergence = FALSE))
        expect_warning(qfrm(A1, A2, alphaB = 1.5,           check_convergence = FALSE))
        expect_warning(qfrm(A1, A2, mu = mu, alphaB = 1.5,  check_convergence = FALSE))
        expect_warning(qfrm(A1, A1s,                        check_convergence = FALSE))
        expect_warning(qfrm(A1, A1s, mu = mu,               check_convergence = FALSE))
        expect_warning(qfmrm(A1, A1s,                       check_convergence = FALSE))
        expect_warning(qfmrm(A1, A1s, mu = mu,              check_convergence = FALSE))
        expect_warning(qfmrm(A1, A2, alphaB = 1.5,          check_convergence = FALSE))
        expect_warning(qfmrm(A1, A2, mu = mu, alphaB = 1.5, check_convergence = FALSE))

        ## Above is dismissed when no error bound is returned
        expect_silent(qfrm(A1, p = 1/2, alphaA = 1.5,      error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1s, p = 1/2,                   error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, alphaB = 1.5,           error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A2, mu = mu, alphaB = 1.5,  error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A1s,                        error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfrm(A1, A1s, mu = mu,               error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A1s,                       error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A1s, mu = mu,              error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, alphaB = 1.5,          error_bound = FALSE, check_convergence = FALSE))
        expect_silent(qfmrm(A1, A2, mu = mu, alphaB = 1.5, error_bound = FALSE, check_convergence = FALSE))

        ## Does not matter when error bound is unavailable
        expect_silent(qfrm(A1, p = 1/2, mu = mu, alphaA = 1.5,  check_convergence = FALSE))
        expect_silent(qfrm(A1s, p = 1/2, mu = mu,               check_convergence = FALSE))

    }
})

test_that("Expect silence or error for existence conditions: qfrm", {
    nvs <- 2:4
    ks <- c(1:3, 1/2)
    m <- 50
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        # L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        # A3 <- diag(L3)
        I <- diag(nv)
        mu <- 1:nv / nv
        # Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
        # A1r <- Q %*% A1 %*% t(Q)
        # A2r <- Q %*% A2 %*% t(Q)
        # A3r <- Q %*% A3 %*% t(Q)

        for(p in ks) {
            for(q in ks) {
                if(nv / 2 + p <= q) {
                    # expect_error(qfrm(A1, I,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE))
                    # expect_error(suppressWarnings(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE)))
                    expect_error(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE))
                    expect_error(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                } else {
                    expect_silent(qfrm(A1, I,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE))
                    expect_silent(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE))
                    expect_silent(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                }
            }
        }
    }
})

test_that("Expect silence, warning, or error for existence conditions: qfmrm", {
    nvs <- 3:4
    ks <- c(1, 2, 1/2)
    m <- 25
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        I <- diag(nv)
        mu <- 1:nv / nv
        # Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
        # A1r <- Q %*% A1 %*% t(Q)
        # A2r <- Q %*% A2 %*% t(Q)
        # A3r <- Q %*% A3 %*% t(Q)

        for(p in ks) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(nv / 2 + p <= q + r) {
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE))
                        expect_error(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                        ## The following seem to yield NaN in this condition, unless nv = 2, against which a warning is thrown
                        expect_warning(expect_warning(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE)))
                        expect_warning(expect_warning(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE)))
                        expect_warning(expect_warning(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE)))
                        expect_warning(expect_warning(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE)))
                    } else {
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE))
                        expect_silent(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE))
                        expect_silent(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE))
                        expect_silent(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                    }
                }
            }
        }
    }
})
