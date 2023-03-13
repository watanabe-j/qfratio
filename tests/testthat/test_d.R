test_that("Expect d_1(A) as tr(A) / 2", {
    nvs <- 2:10
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)

        trA1 <- tr(A1)
        trA2 <- tr(A2)
        trA3 <- tr(A3)

        expect_equal(trA1 / 2, d1_i(L1, 3)[2])
        expect_equal(trA2 / 2, d1_i(L2, 3)[2])
        expect_equal(trA3 / 2, d1_i(L3, 3)[2])
        expect_equal(trA1 / 2, d2_ij_m(A1, A2, 3)[2, 1])
        expect_equal(trA2 / 2, d2_ij_m(A1, A2, 3)[1, 2])
        expect_equal(trA1 / 2, d2_pj_m(A1, A2, 2, 1)[2, 1])
        expect_equal(trA2 / 2, d2_pj_m(A2, A1, 2, 1)[2, 1])
        expect_equal(trA1 / 2, d2_ij_v(L1, L2, 3)[2, 1])
        expect_equal(trA2 / 2, d2_ij_v(L1, L2, 3)[1, 2])
        expect_equal(trA1 / 2, d2_pj_v(L1, L2, 2, 1)[2, 1])
        expect_equal(trA2 / 2, d2_pj_v(L2, L1, 2, 1)[2, 1])
        expect_equal(trA1 / 2, d3_ijk_m(A1, A2, A3, 3)[2, 1, 1])
        expect_equal(trA2 / 2, d3_ijk_m(A1, A2, A3, 3)[1, 2, 1])
        expect_equal(trA3 / 2, d3_ijk_m(A1, A2, A3, 3)[1, 1, 2])
        expect_equal(trA1 / 2, d3_ijk_v(L1, L2, L3, 3)[2, 1, 1])
        expect_equal(trA2 / 2, d3_ijk_v(L1, L2, L3, 3)[1, 2, 1])
        expect_equal(trA3 / 2, d3_ijk_v(L1, L2, L3, 3)[1, 1, 2])
        expect_equal(trA1 / 2, d3_pjk_m(A1, A2, A3, 3, 1)[2, 1, 1])
        expect_equal(trA2 / 2, d3_pjk_m(A2, A1, A3, 3, 1)[2, 1, 1])
        expect_equal(trA3 / 2, d3_pjk_m(A3, A1, A2, 3, 1)[2, 1, 1])

    }
})

test_that("Consistency between d, h, etc.", {
    nvs <- 2:5
    m <- 20
    for(nv in nvs) {
        L1 <- 1:nv
        L2 <- nv:1
        L3 <- sqrt(nv:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        mu <- 1:nv / nv

        expect_equal(dtil2_pq_v(L1, L2, mu, m, m), dtil2_pq_m(A1, A2, mu, m, m))
        expect_equal(dtil3_pqr_v(L1, L2, L3, mu, 2, m, m), dtil3_pqr_m(A1, A2, A3, mu, 2, m, m))

        expect_equal(d2_ij_v(L1, L2, m, 2), d2_pj_v(L1, L2, m, 2))
        expect_equal(d2_ij_v(L1, L2, m, 1), d2_1j_v(L1, L2, m))
        expect_equal(htil2_pj_v(L1, L2, mu = mu, m, 1), htil2_pj_m(A1, A2, mu = mu, m))
        expect_equal(hhat2_pj_v(L1, L2, mu = mu, m, 1), hhat2_pj_m(A1, A2, mu = mu, m))

        expect_equal(d3_ijk_m(A1, A2, A3, m, 2), d3_pjk_m(A1, A2, A3, m, 2))
        expect_equal(d3_ijk_v(L1, L2, L3, m, 2), d3_pjk_v(L1, L2, L3, m, 2))
        expect_equal(d3_pjk_v(L1, L2, L3, m, 2), d3_pjk_m(A1, A2, A3, m, 2))
        expect_equal(htil3_pjk_v(L1, L2, L3, mu, m, 2), htil3_pjk_m(A1, A2, A3, mu, m, 2))
        expect_equal(hhat3_pjk_v(L1, L2, L3, mu, m, 2), hhat3_pjk_m(A1, A2, A3, mu, m, 2))

        expect_equal(d3_pjk_m(A1, A2, A3, m, 2), htil3_pjk_m(A1, A2, A3, rep.int(0, nv), m, 2))
        expect_equal(d3_pjk_v(L1, L2, L3, m, 2), htil3_pjk_v(L1, L2, L3, rep.int(0, nv), m, 2))
        expect_equal(d3_pjk_v(L1, L2, L3, m, 2), htil3_pjk_m(A1, A2, A3, rep.int(0, nv), m, 2))

        expect_equal(d3_pjk_m(A1, A2, A3, m, 2), hhat3_pjk_m(A1, A2, A3, rep.int(0, nv), m, 2))
        expect_equal(d3_pjk_v(L1, L2, L3, m, 2), hhat3_pjk_v(L1, L2, L3, rep.int(0, nv), m, 2))
        expect_equal(d3_pjk_v(L1, L2, L3, m, 2), hhat3_pjk_m(A1, A2, A3, rep.int(0, nv), m, 2))

    }
})
