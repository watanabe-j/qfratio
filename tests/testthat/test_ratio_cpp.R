tol <- 5e-7

test_that("Expect equal from R and Cpp methods: qfrm", {
    nvs <- 2:4
    ks <- c(1:3, 1/2)
    m <- 3
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        I <- diag(nv)
        mu <- 1:nv / nv
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A2r <- Q %*% A2 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks[ks %% 1 == 0]) {
            for(q in ks) {
                if(nv / 2 + p <= q) next
                expect_equal(qfrm(A1, I,  p, q, m = m, use_cpp = FALSE),
                             qfrm(A1, I,  p, q, m = m, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                if(requireNamespace("gsl", quietly = TRUE)) {
                    expect_equal(suppressMessages(qfrm(A1, I,  p, q, m = m, mu = mu, use_cpp = FALSE)),
                                 suppressMessages(qfrm(A1, I,  p, q, m = m, mu = mu, use_cpp = TRUE, nthreads = 1)), tolerance = tol)
                }
                expect_equal(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
            }
        }
        for(p in ks[ks %% 1 != 0]) {
            for(q in ks) {
                if(nv / 2 + p <= q) next
                expect_equal(qfrm(A1, I,  p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, I,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(suppressMessages(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE)),
                             suppressMessages(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1)), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
            }
        }
    }
})


test_that("Expect equal from R and Cpp methods: qfmrm", {
    nvs <- 2:4
    ks <- c(1, 2, 1/2)
    m <- 3
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
            A2r <- Q %*% A2 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(nv / 2 + p <= q + r) next
                    expect_equal(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, I, p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2r, I, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, I, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2r, I, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(I, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(I, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(I, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(I, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1), tolerance = tol)
                }
            }
        }
    }
})

test_that("Expect equal from R and Cpp methods: qfm", {
    nvs <- 2:4
    ks <- c(1, 2)
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
            A2r <- Q %*% A2 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks) {
            expect_equal(qfm_Ap_int(A1, p,          use_cpp = FALSE),
                         qfm_Ap_int(A1, p,          use_cpp = TRUE), tolerance = tol)
            expect_equal(qfm_Ap_int(A1, p, mu = mu, use_cpp = FALSE),
                         qfm_Ap_int(A1, p, mu = mu, use_cpp = TRUE), tolerance = tol)
            for(q in ks) {
                expect_equal(qfpm_ABpq_int(A1, A2,  p, q,          use_cpp = FALSE),
                             qfpm_ABpq_int(A1, A2,  p, q,          use_cpp = TRUE), tolerance = tol)
                expect_equal(qfpm_ABpq_int(A1, A2,  p, q, mu = mu, use_cpp = FALSE),
                             qfpm_ABpq_int(A1, A2,  p, q, mu = mu, use_cpp = TRUE), tolerance = tol)
                expect_equal(qfpm_ABpq_int(A1, A2r, p, q,          use_cpp = FALSE),
                             qfpm_ABpq_int(A1, A2r, p, q,          use_cpp = TRUE), tolerance = tol)
                expect_equal(qfpm_ABpq_int(A1, A2r, p, q, mu = mu, use_cpp = FALSE),
                             qfpm_ABpq_int(A1, A2r, p, q, mu = mu, use_cpp = TRUE), tolerance = tol)
                for(r in ks) {
                    expect_equal(qfpm_ABDpqr_int(A1, A2,  A3, p, q, r,          use_cpp = FALSE),
                                 qfpm_ABDpqr_int(A1, A2,  A3, p, q, r,          use_cpp = TRUE), tolerance = tol)
                    expect_equal(qfpm_ABDpqr_int(A1, A2,  A3, p, q, r, mu = mu, use_cpp = FALSE),
                                 qfpm_ABDpqr_int(A1, A2,  A3, p, q, r, mu = mu, use_cpp = TRUE), tolerance = tol)
                    expect_equal(qfpm_ABDpqr_int(A1, A2r, A3, p, q, r,          use_cpp = FALSE),
                                 qfpm_ABDpqr_int(A1, A2r, A3, p, q, r,          use_cpp = TRUE), tolerance = tol)
                    expect_equal(qfpm_ABDpqr_int(A1, A2r, A3, p, q, r, mu = mu, use_cpp = FALSE),
                                 qfpm_ABDpqr_int(A1, A2r, A3, p, q, r, mu = mu, use_cpp = TRUE), tolerance = tol)
                }
            }
        }
    }
})

test_that("Expect equal from double and long double in C++: qfrm", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 3
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        I <- diag(nv)
        mu <- 1:nv / nv
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A2r <- Q %*% A2 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks[ks %% 1 != 0]) {
            for(q in ks) {
                if(nv / 2 + p <= q) next
                expect_equal(suppressMessages(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double")),
                             suppressMessages(qfrm(A1, I,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double")), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
            }
        }
    }
})

test_that("Expect equal from double and long double in C++: qfmrm", {
    nvs <- 2:4
    ks <- c(1, 2, 1/2)
    m <- 3
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
            A2r <- Q %*% A2 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(nv / 2 + p <= q + r) next
                    expect_equal(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, I, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, I, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, I, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, I, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(I, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(I, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "long_double"), tolerance = tol)
                }
            }
        }
    }
})

test_that("Expect equal from double and coef_wise in C++: qfrm", {
    nvs <- 2:4
    ks <- c(1:3, 1/2, 3/2)
    m <- 3
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        I <- diag(nv)
        mu <- 1:nv / nv
        if(requireNamespace("stats", quietly = TRUE)) {
            Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
            A2r <- Q %*% A2 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks[ks %% 1 != 0]) {
            for(q in ks) {
                if(nv / 2 + p <= q) next
                expect_equal(suppressWarnings(qfrm(A1, I,  p, q, m = m, mu = mu, use_cpp = TRUE, nthreads = 1, cpp_method = "double")),
                             suppressWarnings(qfrm(A1, I,  p, q, m = m, mu = mu, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise")), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                expect_equal(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                expect_equal(qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                             qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
            }
        }
    }
})

test_that("Expect equal from double and coef_wise in C++: qfmrm", {
    nvs <- 2:4
    ks <- c(1, 2, 1/2)
    m <- 3
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
            # A1r <- Q %*% A1 %*% t(Q)
            A2r <- Q %*% A2 %*% t(Q)
            # A3r <- Q %*% A3 %*% t(Q)
        } else {
            A2r <- A2
        }

        for(p in ks) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(nv / 2 + p <= q + r) next
                    expect_equal(qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, I,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, I,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, I, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, I, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, I, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, I, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(I, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(I, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(I, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "double"),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE, nthreads = 1, cpp_method = "coef_wise"), tolerance = tol)
                }
            }
        }
    }
})
