test_that("Expect d_1(A) as tr(A) / 2", {
    ps <- 2:10
    # ks <- 1:5
    for(p in ps) {
        L1 <- 1:p
        L2 <- p:1
        L3 <- sqrt(p:1)
        A1 <- diag(L1)
        A2 <- diag(L2)
        A3 <- diag(L3)
        # mu <- 1:p / p

        trA1 <- tr(A1)
        trA2 <- tr(A2)
        trA3 <- tr(A3)

        expect_equal(trA1 / 2, d1_i(L1)[2])
        expect_equal(trA2 / 2, d1_i(L2)[2])
        expect_equal(trA3 / 2, d1_i(L3)[2])
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
