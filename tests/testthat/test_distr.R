tol <- 5e-7
tol2 <- 5e-4
tol3 <- 5e-3
CQF_available <- requireNamespace("CompQuadForm", quietly = TRUE)

test_that("Expect identical results for simultaneously rotated matrices", {
    nvs <- 2:4
    for(nv in nvs) {
        qs <- 0:nv + 0.5
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2
        Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
        A1r <- Q %*% A1 %*% t(Q)
        A2r <- Q %*% A2 %*% t(Q)
        mur <- Q %*% mu

        expect_equal(pqfr(qs, A1,  A2,  mu = mu,  method = "imhof"),
                     pqfr(qs, A1r, A2r, mu = mur, method = "imhof"), tolerance = tol)
        expect_equal(pqfr(qs, A1,  A2,  mu = mu,  method = "forchini", m = 5, check_convergence = FALSE),
                     pqfr(qs, A1r, A2r, mu = mur, method = "forchini", m = 5, check_convergence = FALSE), tolerance = tol)
        expect_equal(pqfr(qs, A1,  A2,  mu = mu,  method = "butler"),
                     pqfr(qs, A1r, A2r, mu = mur, method = "butler"), tolerance = tol)
        expect_equal(dqfr(qs, A1,  A2,  mu = mu,  method = "broda"),
                     dqfr(qs, A1r, A2r, mu = mur, method = "broda"), tolerance = tol)
        expect_equal(dqfr(qs, A1,                 method = "hillier", m = 5, check_convergence = FALSE),
                     dqfr(qs, A1r,                method = "hillier", m = 5, check_convergence = FALSE), tolerance = tol)
        expect_equal(dqfr(qs, A1,  A2,  mu = mu,  method = "butler"),
                     dqfr(qs, A1r, A2r, mu = mur, method = "butler"), tolerance = tol)
    }
})

test_that("Expect equal from R and C++ methods", {
    nvs <- 2:4
    for(nv in nvs) {
        qs <- 0:nv + 0.5
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2

        if(CQF_available) {
            expect_equal(pqfr(qs, A1, A2, mu = mu, method = "imhof", return_abserr_attr = FALSE, use_cpp = TRUE),
                         pqfr(qs, A1, A2, mu = mu, method = "imhof", return_abserr_attr = FALSE, use_cpp = FALSE), tolerance = tol)
        }
        expect_equal(pqfr(qs, A1, A2, mu = mu, method = "forchini", m = 5, check_convergence = FALSE, use_cpp = TRUE),
                     pqfr(qs, A1, A2, mu = mu, method = "forchini", m = 5, check_convergence = FALSE, use_cpp = FALSE), tolerance = tol)
        expect_equal(pqfr(qs, A1, A2, mu = mu, method = "butler", use_cpp = TRUE),
                     pqfr(qs, A1, A2, mu = mu, method = "butler", use_cpp = FALSE), tolerance = tol)
        expect_equal(pqfr(qs, A1, A2, mu = mu, method = "butler", order_spa = 1, use_cpp = TRUE),
                     pqfr(qs, A1, A2, mu = mu, method = "butler", order_spa = 1, use_cpp = FALSE), tolerance = tol)
        expect_equal(dqfr(qs, A1, A2, mu = mu, method = "broda", return_abserr_attr = FALSE, use_cpp = TRUE),
                     dqfr(qs, A1, A2, mu = mu, method = "broda", return_abserr_attr = FALSE, use_cpp = FALSE), tolerance = tol)
        expect_equal(dqfr(qs, A1,              method = "hillier", m = 5, check_convergence = FALSE, use_cpp = TRUE),
                     dqfr(qs, A1,              method = "hillier", m = 5, check_convergence = FALSE, use_cpp = FALSE), tolerance = tol)
        expect_equal(dqfr(qs, A1, A2, mu = mu, method = "butler", use_cpp = TRUE),
                     dqfr(qs, A1, A2, mu = mu, method = "butler", use_cpp = FALSE), tolerance = tol)
        expect_equal(dqfr(qs, A1, A2, mu = mu, method = "butler", order_spa = 1, use_cpp = TRUE),
                     dqfr(qs, A1, A2, mu = mu, method = "butler", order_spa = 1, use_cpp = FALSE), tolerance = tol)
    }
})

test_that("Expect similar results between different methods", {
    nvs <- 2:4
    for(nv in nvs) {
        qs <- 0:nv + 0.5
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2
        pseq_imhof    <- pqfr(qs, A1, A2, mu = mu, method = "imhof", return_abserr_attr = FALSE)
        pseq_forchini <- pqfr(qs, A1, A2, mu = mu, method = "forchini", m = 50, check_convergence = FALSE)
        dseq_broda    <- dqfr(qs, A1, method = "broda", return_abserr_attr = FALSE)
        dseq_hillier  <- dqfr(qs, A1, method = "hillier", m = 50, check_convergence = FALSE)

        expect_equal(pseq_imhof, pseq_forchini, tolerance = tol2)
        if(CQF_available) {
            pseq_davies <- pqfr(qs, A1, A2, mu = mu, method = "davies")
            expect_equal(pseq_imhof, pseq_davies,   tolerance = tol2)
        }
        expect_equal(dseq_broda, dseq_hillier,  tolerance = tol2)
    }
})

test_that("Expect equal p-values with exponents with nnd matrices", {
    nvs <- 2:4
    for(nv in nvs) {
        qs <- 0:nv + 0.2
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2

        res_p1 <- pqfr(qs,   A1, A2, 1, mu = mu, return_abserr_attr = FALSE)
        res_p2 <- pqfr(qs^2, A1, A2, 2, mu = mu, return_abserr_attr = FALSE)
        res_p3 <- pqfr(qs^3, A1, A2, 3, mu = mu, return_abserr_attr = FALSE)
        res_p05 <- pqfr(qs^(1/2), A1, A2, 1/2, mu = mu, return_abserr_attr = FALSE)

        expect_equal(res_p1, res_p2, tolerance = tol)
        expect_equal(res_p1, res_p3, tolerance = tol)
        expect_equal(res_p1, res_p05, tolerance = tol)
    }
})

test_that("Expect equal p-values with odd exponents with indefinite matrices", {
    nvs <- 2:4
    for(nv in nvs) {
        qs <- (-2):nv + 0.2
        L1 <- (1:nv - 1.5) * 2
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2

        res_p1 <- pqfr(qs,   A1, A2, 1, mu = mu, return_abserr_attr = FALSE)
        res_p3 <- pqfr(qs^3, A1, A2, 3, mu = mu, return_abserr_attr = FALSE)

        expect_equal(res_p1, res_p3, tolerance = tol)
    }
})

test_that("Expect unity when density is integrated with nnd matrices", {
    ## Tests with nv = 2, 3 often hit singularity, causing error
    nvs <- 4:5
    for(nv in nvs) {
        qs <- 0:nv + 0.2
        L1 <- 1:nv
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2
        qmin <- min(L1 / L2)
        qmax <- max(L1 / L2)
        
        res_p1  <- stats::integrate(dqfr, qmin^1, qmax^1, A1, A2, p = 1, mu = mu)$value
        res_p2  <- stats::integrate(dqfr, qmin^2, qmax^2, A1, A2, p = 2, mu = mu)$value
        res_p3  <- stats::integrate(dqfr, qmin^3, qmax^3, A1, A2, p = 3, mu = mu)$value
        res_p05 <- stats::integrate(dqfr, qmin^(1/2), qmax^(1/2), A1, A2, p = 0.5, mu = mu)$value

        expect_equal(res_p1,  1, tolerance = tol3)
        expect_equal(res_p2,  1, tolerance = tol3)
        expect_equal(res_p3,  1, tolerance = tol3)
        expect_equal(res_p05, 1, tolerance = tol3)
    }
})

test_that("Expect unity when density is integrated with indefinite matrices", {
    nvs <- 4:5
    for(nv in nvs) {
        qs <- 0:nv + 0.2
        L1 <- (1:nv - 1.5) * 2
        L2 <- nv:1
        A1 <- diag(L1)
        A2 <- diag(L2)
        mu <- 1:nv * 0.2
        qmin <- min(L1 / L2)
        qmax <- max(L1 / L2)

        res_p1 <- stats::integrate(dqfr, qmin^1, qmax^1, A1, A2, p = 1, mu = mu)$value
        res_p2 <- stats::integrate(dqfr, 0,      qmax^2, A1, A2, p = 2, mu = mu)$value
        res_p3 <- stats::integrate(dqfr, qmin^3, qmax^3, A1, A2, p = 3, mu = mu)$value

        expect_equal(res_p1, 1, tolerance = tol3)
        expect_equal(res_p2, 1, tolerance = tol3)
        expect_equal(res_p3, 1, tolerance = tol3)
    }
})
