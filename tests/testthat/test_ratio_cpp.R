test_that("Expect equal from R anc Cpp methods: qfrm", {
    nvs <- 2:5
    ks <- c(1:3, 1/2)
    m <- 100
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        # L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        # A3 <- diag(L3)
        mu <- 1:nv / nv
        Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
        # A1r <- Q %*% A1 %*% t(Q)
        A2r <- Q %*% A2 %*% t(Q)
        # A3r <- Q %*% A3 %*% t(Q)

        for(p in ks) {
            for(q in ks) {
                if(nv / 2 + p <= q) next
                expect_equal(qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2,  p, q, m = m, check_convergence = FALSE, use_cpp = TRUE))
                expect_equal(qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2,  p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                expect_equal(qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2r, p, q, m = m, check_convergence = FALSE, use_cpp = TRUE))
                expect_equal(qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                             qfrm(A1, A2r, p, q, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
            }
        }
    }
})


test_that("Expect equal from R anc Cpp methods: qfmrm", {
    nvs <- 2:5
    ks <- c(1:3, 1/2)
    m <- 50
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        mu <- 1:nv / nv
        Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
        # A1r <- Q %*% A1 %*% t(Q)
        A2r <- Q %*% A2 %*% t(Q)
        # A3r <- Q %*% A3 %*% t(Q)

        for(p in ks) {
            for(q in ks/2) {
                for(r in ks/2) {
                    if(nv / 2 + p <= q + r) next
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE))
                    expect_equal(qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2, A3,  p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, check_convergence = FALSE, use_cpp = TRUE))
                    expect_equal(qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
                                 qfmrm(A1, A2r, A3, p, q, r, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
                }
            }
        }
    }
})

# test_that("Expect equal from R anc Cpp methods: misc", {
#     nvs <- 2:5
#     m <- 60
#     for(nv in nvs) {
#         L1 <- 1:nv
#         L2 <- nv:1
#         L3 <- sqrt(nv:1)
#         A1 <- diag(L1)
#         A2 <- diag(L2)
#         A3 <- diag(L3)
#         mu <- 1:nv / nv
#         Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
#         A1r <- Q %*% A1 %*% t(Q)
#         Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
#         A2r <- Q %*% A2 %*% t(Q)
#         Q <- qr.Q(qr(matrix(stats::rnorm(nv^2), nv, nv)))
#         A3r <- Q %*% A3 %*% t(Q)
#
#         expect_equal(qfm_Ap_int(A1, 3, check_convergence = FALSE, use_cpp = FALSE),
#                      qfm_Ap_int(A1, 3, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfm_Ap_int(A1r, 3, check_convergence = FALSE, use_cpp = FALSE),
#                      qfm_Ap_int(A1r, 3, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfm_Ap_int(A1, 4, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfm_Ap_int(A1, 4, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfm_Ap_int(A1r, 3, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfm_Ap_int(A1r, 3, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfpm_ABpq_int(A1, A2, 2, 3, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABpq_int(A1, A2, 2, 3, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfpm_ABpq_int(A1r, A2, 2, 3, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABpq_int(A1r, A2, 2, 3, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfpm_ABpq_int(A1, A2, 3, 4, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABpq_int(A1, A2, 3, 4, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfpm_ABpq_int(A1r, A2, 2, 3, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABpq_int(A1r, A2, 2, 3, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfpm_ABDpqr_int(A1, A2, A3, 2, 3, 4, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABDpqr_int(A1, A2, A3, 2, 3, 4, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfpm_ABDpqr_int(A1r, A2, A3, 2, 2, 1, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABDpqr_int(A1r, A2, A3, 2, 2, 1, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfpm_ABDpqr_int(A1, A2, A3, 2, 3, 4, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABDpqr_int(A1, A2, A3, 2, 3, 4, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfpm_ABDpqr_int(A1r, A2, A3, 2, 3, 4, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfpm_ABDpqr_int(A1r, A2, A3, 2, 3, 4, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfrm_ApIq_int(A1r, 2, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApIq_int(A1r, 2, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfrm_ApIq_int(A1r, 2, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApIq_int(A1r, 2, mu = mu, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#
#         expect_equal(qfrm_ApIq_npi(A1, -1, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApIq_npi(A1, -1, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApIq_npi(A1, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApIq_npi(A1, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(suppressWarnings(qfrm_ApIq_npi(A1, -1, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE)),
#                      suppressWarnings(qfrm_ApIq_npi(A1, -1, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE)),
#                      ignore_attr = TRUE)
#         expect_equal(suppressWarnings(qfrm_ApIq_npi(A1, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE)),
#                      suppressWarnings(qfrm_ApIq_npi(A1, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE)),
#                      ignore_attr = TRUE)
#
#         expect_equal(qfrm_ApBq_int(A1, A3, 2, 2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_int(A1, A3, 2, 2, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApBq_int(A1, A3r, 2, 2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_int(A1, A3r, 2, 2, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApBq_int(A1, A3, 2, 2, mu = mu, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_int(A1, A3, 2, 2, mu = mu, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApBq_int(A1, A3r, 2, 2, mu = mu, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_int(A1, A3r, 2, 2, mu = mu, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#
#         expect_equal(qfrm_ApBq_npi(A1, crossprod(A1), 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_npi(A1, crossprod(A1), 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApBq_npi(A1, A3r, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_npi(A1, A3r, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApBq_npi(A1, crossprod(A1), 1/2, 1/2, mu = mu, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_npi(A1, crossprod(A1), 1/2, 1/2, mu = mu, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#         expect_equal(qfrm_ApBq_npi(A1, A3r, 1/2, 1/2, mu = mu, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfrm_ApBq_npi(A1, A3r, 1/2, 1/2, mu = mu, m = m, check_convergence = FALSE, use_cpp = TRUE),
#                      ignore_attr = TRUE)
#
#         expect_equal(qfmrm_ApBIqr_int(A1, crossprod(A1), 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_int(A1, crossprod(A1), 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBIqr_int(A1, A3r, 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_int(A1, A3r, 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBIqr_int(A1, crossprod(A1), 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_int(A1, crossprod(A1), 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBIqr_int(A1, A3r, 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_int(A1, A3r, 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfmrm_ApBIqr_npi(A1, crossprod(A1), 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_npi(A1, crossprod(A1), 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBIqr_npi(A1, A3r, 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_npi(A1, A3r, 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBIqr_npi(A1, crossprod(A1), 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_npi(A1, crossprod(A1), 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBIqr_npi(A1, A3r, 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBIqr_npi(A1, A3r, 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfmrm_IpBDqr_gen(A1, diag(1 / L1), 2, 1, 1, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_IpBDqr_gen(A1, diag(1 / L1), 2, 1, 1, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_IpBDqr_gen(A1, A3r, 2, 1, 1, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_IpBDqr_gen(A1, A3r, 2, 1, 1, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_IpBDqr_gen(A1, diag(1 / L1), 2, 1, 1, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_IpBDqr_gen(A1, diag(1 / L1), 2, 1, 1, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_IpBDqr_gen(A1, A3r, 2, 1, 1, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_IpBDqr_gen(A1, A3r, 2, 1, 1, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfmrm_ApBDqr_int(A1 * A3, A1^2, A3^2, 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_int(A1 * A3, A1^2, A3^2, 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBDqr_int(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_int(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1, 1/2, 1/2, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBDqr_int(A1 * A3, A1^2, A3^2, 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_int(A1 * A3, A1^2, A3^2, 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBDqr_int(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_int(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1, 1/2, 1/2, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#         expect_equal(qfmrm_ApBDqr_npi(A1 * A3, A1^2, A3^2, 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_npi(A1 * A3, A1^2, A3^2, 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBDqr_npi(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_npi(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1/2, 1/4, 1/4, m = m, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBDqr_npi(A1 * A3, A1^2, A3^2, 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_npi(A1 * A3, A1^2, A3^2, 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#         expect_equal(qfmrm_ApBDqr_npi(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = FALSE),
#                      qfmrm_ApBDqr_npi(A1 %*% A3r, A1 ^ 2, crossprod(A3r), 1/2, 1/4, 1/4, m = m, mu = mu, check_convergence = FALSE, use_cpp = TRUE))
#
#     }
# })
