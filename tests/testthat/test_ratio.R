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
