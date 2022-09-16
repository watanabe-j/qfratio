d1_i <- function(L, m = 100L) {
    n <- length(L)
    dks <- rep.int(c(1, 0), c(1L, m))
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    uk <- rep.int(0, n)
    for(k in seq_len(m)) {
        uk <- L * (dks[k] + uk)
        dks[k + 1L] <- sum(uk) / (2 * k)
        if(max(uk) > thr) {
            dks <- dks / 1e10
            uk <- uk / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

#' Recursion for dtilde_k
#'
#' \code{dtil1_i()} calculates \eqn{\tilde{d}_k} by the recursion in
#' eqs. 28--30 in Hillier et al. (2014)
#'
dtil1_i_v <- function(L, mu = rep.int(0, n), m = 100L) {
    n <- length(L)
    D <- mu ^ 2
    dks <- rep.int(c(1, 0), c(1L, m))
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    uk <- rep.int(0, n)
    vk <- rep.int(0, n)
    for(k in seq_len(m)) {
        uk <- L * (dks[k] + uk)
        vk <- D * uk + L * vk
        dks[k + 1L] <- sum(uk + vk) / (2 * k)
        if(max(uk) > thr || max(vk) > thr) {
            dks <- dks / 1e10
            uk <- uk / 1e10
            vk <- vk / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

#' Recursion for dtilde_k
#'
#' \code{dtil1_i()} calculates \eqn{\tilde{d}_k} by the recursion in
#' eqs. 28--30 in Hillier et al. (2014)
#'
dtil1_i_m <- function(A, mu = rep.int(0, n), m = 100L) {
    eigA <- eigen(A, symmetric = TRUE)
    L <- eigA$values
    D <- c(crossprod(eigA$vectors, mu)) ^ 2
    n <- length(L)
    dks <- rep.int(c(1, 0), c(1L, m))
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    uk <- rep.int(0, n)
    vk <- rep.int(0, n)
    for(k in seq_len(m)) {
        uk <- L * (dks[k] + uk)
        vk <- D * uk + L * vk
        dks[k + 1L] <- sum(uk + vk) / (2 * k)
        if(max(uk) > thr || max(vk) > thr) {
            dks <- dks / 1e10
            uk <- uk / 1e10
            vk <- vk / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

dtil2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m,
                         fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    hs <- list(zerovec)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        hs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            ht <- tG %*% mu
            if(i1 >= 1L) ht <- ht + A1 %*% hs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) ht <- ht + A2 %*% hs[[il2(i1, i2 - 1L)]]
            hs[[il2(i1, i2)]] <- ht
            dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, ht))) / (2 * k)
        }
        if(max(tG) > thr || max(th)) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            hs <- lapply(hs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

dtil2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m,
                         fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    hs <- list(zeros)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        hs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            ht <- tG * mu
            if(i1 >= 1L) ht <- ht + L1 * hs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) ht <- ht + L2 * hs[[il2(i1, i2 - 1L)]]
            hs[[il2(i1, i2)]] <- ht
            dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * ht)) / (2 * k)
        }
        if(max(tG) > thr || max(ht) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            hs <- lapply(hs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}


dtil3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    hs <- list(zerovec)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        hs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                th <- tG %*% mu
                if(i1 >= 1L) th <- th + A1 %*% hs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) th <- th + A2 %*% hs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) th <- th + A3 %*% hs[[il3(i1, i2, i3 - 1L)]]
                hs[[il3(i1, i2, i3)]] <- th
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, th))) / (2 * k)
            }
        }
        if(max(tG) > thr || max(th) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            hs <- lapply(hs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

dtil3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    hs <- list(zeros)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        hs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                th <- tG * mu
                if(i1 >= 1L) th <- th + L1 * hs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) th <- th + L2 * hs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) th <- th + L3 * hs[[il3(i1, i2, i3 - 1L)]]
                hs[[il3(i1, i2, i3)]] <- th
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * th)) / (2 * k)
            }
        }
        if(max(tG) > thr || max(th) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            hs <- lapply(hs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

dtil2_pq_m <- function(A1, A2, mu = rep.int(0, n), p = 1L, q = 1L) {
    if(p == 1L) return(dtil2_1q_m(A1, A2, mu, q))
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, p + 1L, q + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(p + 1L)] <- list(zeromat)
    g_k_i[1L:(p + 1L)] <- list(zerovec)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu + A1 %*% g_k_i[[i]]
        dks[i + 1L, 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * i)
    }
    for(k in 1L:q) {
        G_k_i[[1L]] <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        g_k_i[[1L]] <- G_k_i[[1L]] %*% mu + A2 %*% g_k_i[[1L]]
        dks[1L, k + 1L] <- (tr(G_k_i[[1L]]) + c(crossprod(mu, g_k_i[[1L]]))) / (2 * k)
        for(i in 1L:p) {
            G_k_i[[i + 1L]] <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu +
                               A1 %*% g_k_i[[i]] + A2 %*% g_k_i[[i + 1L]]
            dks[i + 1L, k + 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

dtil2_1q_m <- function(A1, A2, mu = rep.int(0, n), q = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, 2L, q + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1 # * (dks[1L, 1L] + G_k_0)
    g_k_0 <- matrix(0, n, 1)
    g_k_1 <- G_k_1 %*% mu # + A1 %*% g_k_0
    dks[2L, 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / 2
    for(k in 1L:q) {
        G_k_0 <- A2 %*% (dks[1L, k] * In + G_k_0)
        g_k_0 <- G_k_0 %*% mu + A2 %*% g_k_0
        dks[1L, k + 1L] <- (tr(G_k_0) + c(crossprod(mu, g_k_0))) / (2 * k)
        G_k_1 <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        g_k_1 <- G_k_1 %*% mu + A1 %*% g_k_0 + A2 %*% g_k_1
        dks[2L, k + 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / (2 * (k + 1))
        if(max(G_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

dtil2_pq_v <- function(L1, L2, mu = rep.int(0, n), p = 1L, q = 1L) {
    if(p == 1L) return(dtil2_1j_v(L1, L2, mu, q))
    n <- length(L1)
    dks <- matrix(0, p + 1L, q + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(p + 1L)] <- list(zeros)
    g_k_i[1L:(p + 1L)] <- list(zeros)
    for(i in 1L:p) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu + L1 * g_k_i[[i]]
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:q) {
        G_k_i[[1L]] <- L2 * (dks[1L, k] + G_k_i[[1L]])
        g_k_i[[1L]] <- G_k_i[[1L]] * mu + L2 * g_k_i[[1L]]
        dks[1L, k + 1L] <- sum(G_k_i[[1L]] + mu * g_k_i[[1L]]) / (2 * k)
        for(i in 1L:p) {
            G_k_i[[i + 1L]] <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu +
                               L1 * g_k_i[[i]] + L2 * g_k_i[[i + 1L]]
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

dtil2_1q_v <- function(L1, L2, mu = rep.int(0, n), q = 1L) {
    n <- length(L1)
    dks <- matrix(0, 2L, q + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1 # * (dks[1L, 1L] + G_k_0)
    g_k_0 <- rep.int(0, n)
    g_k_1 <- G_k_1 * mu # + L1 %*% g_k_0
    dks[2L, 1L] <- sum(G_k_1 + mu * g_k_1) / 2
    for(k in 1L:m) {
        G_k_0 <- L2 * (dks[1L, k] + G_k_0)
        g_k_0 <- G_k_0 * mu + L2 * g_k_0
        dks[1L, k + 1L] <- sum(G_k_0 + mu * g_k_0) / (2 * k)
        G_k_1 <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        g_k_1 <- G_k_1 * mu + L1 * g_k_0 + L2 * g_k_1
        dks[2L, k + 1L] <- sum(G_k_1 + mu * g_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}


dtil3_pqr_m <- function(A1, A2, A3, mu = rep.int(0, n), p = 1L, q = 1L, r = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    m <- q + r
    dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    Gc <- list()
    gc <- list()
    Gc[1L:(p + 1L)] <- list(zeromat)
    gc[1L:(p + 1L)] <- list(zerovec)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu + A1 %*% gc[[i]]
        dks[i + 1L, 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        if(k <= q) {
            Gc[[1L]] <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
            gc[[1L]] <- Gc[[1L]] %*% mu + A2 %*% go[[1L]][[1L]]
            dks[1L, k + 1L, 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                                A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu +
                                   A1 %*% gc[[i]] + A2 %*% go[[1L]][[i + 1L]]
                dks[i + 1L, k + 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
            }
            Gn <- list(Gc)
            gn <- list(gc)
        } else {
            Gn <- list(0)
            gn <- list(0)
        }
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                if(k - j <= q && j <= r) {
                    Gc[[1L]] <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                    gc[[1L]] <- Gc[[1L]] %*% mu + A2 %*% go[[j + 1L]][[1L]] + A3 %*% go[[j]][[1L]]
                    dks[1L, k - j + 1L, j + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
                    for(i in 1L:p) {
                        Gc[[i + 1L]] <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                                        A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                                        A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu +
                                        A1 %*% gc[[i]] + A2 %*% go[[j + 1L]][[i + 1L]] + A3 %*% go[[j]][[i + 1L]]
                        dks[i + 1L, k - j + 1L, j + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
                    }
                    Gn <- c(Gn, list(Gc))
                    gn <- c(gn, list(gc))
                } else {
                    Gn <- c(Gn, list(0))
                    gn <- c(gn, list(0))
                }
            }
        }
        if(k <= r) {
            Gc[[1L]] <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
            gc[[1L]] <- Gc[[1L]] %*% mu + A3 %*% go[[k]][[1L]]
            dks[1L, 1L, k + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                                A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu +
                                A1 %*% gc[[i]] + A3 %*% go[[k]][[i + 1L]]
                dks[i + 1L, 1L, k + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
            }
            Gn <- c(Gn, list(Gc))
            gn <- c(gn, list(gc))
        } else {
            Gn <- c(Gn, list(0))
            gn <- c(gn, list(0))
        }
        if(max(unlist(Gc)) > thr || max(unlist(gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

dtil3_pqr_v <- function(L1, L2, L3, mu = rep.int(0, n), p = 1L, q = 1L, r = 1L) {
    n <- length(L1)
    m <- q + r
    dks <- array(0, dim = c(p + 1L, q + 1L, r + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    Gc <- list()
    gc <- list()
    Gc[1L:(p + 1L)] <- list(zeros)
    gc[1L:(p + 1L)] <- list(zeros)
    for(i in 1L:p) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] * mu + L1 * gc[[i]]
        dks[i + 1L, 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        if(k <= q) {
            Gc[[1L]] <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
            gc[[1L]] <- Gc[[1L]] * mu + L2 * go[[1L]][[1L]]
            dks[1L, k + 1L, 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                                L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] * mu +
                                   L1 * gc[[i]] + L2 * go[[1L]][[i + 1L]]
                dks[i + 1L, k + 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
            }
            Gn <- list(Gc)
            gn <- list(gc)
        } else {
            Gn <- list(0)
            gn <- list(0)
        }
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                if(k - j <= q && j <= r) {
                    Gc[[1L]] <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                    gc[[1L]] <- Gc[[1L]] * mu + L2 * go[[j + 1L]][[1L]] + L3 * go[[j]][[1L]]
                    dks[1L, k - j + 1L, j + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
                    for(i in 1L:p) {
                        Gc[[i + 1L]] <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                                        L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                                        L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                        gc[[i + 1L]] <- Gc[[i + 1L]] * mu +
                                        L1 * gc[[i]] + L2 * go[[j + 1L]][[i + 1L]] + L3 * go[[j]][[i + 1L]]
                        dks[i + 1L, k - j + 1L, j + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
                    }
                    Gn <- c(Gn, list(Gc))
                    gn <- c(gn, list(gc))
                } else {
                    Gn <- c(Gn, list(0))
                    gn <- c(gn, list(0))
                }
            }
        }
        if(k <= r) {
            Gc[[1L]] <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
            gc[[1L]] <- Gc[[1L]] * mu + L3 * go[[k]][[1L]]
            dks[1L, 1L, k + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
            for(i in 1L:p) {
                Gc[[i + 1L]] <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                                L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
                gc[[i + 1L]] <- Gc[[i + 1L]] * mu +
                                L1 * gc[[i]] + L3 * go[[k]][[i + 1L]]
                dks[i + 1L, 1L, k + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
            }
            Gn <- c(Gn, list(Gc))
            gn <- c(gn, list(gc))
        } else {
            Gn <- c(Gn, list(0))
            gn <- c(gn, list(0))
        }
        if(max(unlist(Gc)) > thr || max(unlist(gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}



#' Recursion for a_{r,l}
#'
#' \code{arl()} calculates \eqn{a_{k,l}} by the recursion in
#' eqs. 31--32 in Hillier et al. (2014)
#'
#' Note that \eqn{w_{r,i}} there should be understood as \eqn{w_{r,l,i}} with
#' the index \eqn{l} fixed for each \eqn{a_{r,l}}.
#'
#' @param L
#'   Eigenvalues of the argument matrix; vectors of \eqn{\lambda_i}
#' @param D
#'   Squared norm of the mean vector projected on the eigenvalues of
#'   the argument matrix: vectors of \eqn{\delta_i}
#' @param m
#'   Scalar to specify the desired order
#'
arl <- function(L, ...) {
    UseMethod("arl", L)
}

#' @exportS3Method
#'
arl.matrix <- function(A, mu, m = 10L) {
    if(any(dim(A) == 1L)) return(arl.default(A, D, m))
    eigA <- eigen(A, symmetric = TRUE)
    L <- eigA$values
    D <- c(crossprod(eigA$vectors, mu)) ^ 2
    return(arl.default(L, D, m))
}

#' @exportS3Method
#'
arl.default <-  function(L, D, m = 10L) {
    n <- length(L)
    arls <- matrix(0, m + 1L, m + 1L)
    arls[, 1L] <- d1_i(L, m)
    wrls <- matrix(0, n, m)
    for(k in 1:m) {
        for(l in 1:k) {
            wrls[, l] <- L * (arls[k, l] + wrls[, l])
            arls[k + 1L, l + 1L] <- sum(D * wrls[, l])
        }
    }
    return(arls)
}


# d2_ij_mb <- function(A1, A2, m = 100L, m1 = m, m2 = m, fill_all = !missing(m1) || !missing(m2)) {
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
#     Gs <- array(0, dim = c(n, n, dim(dks)))
#     kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
#     for(k in 1L:kmax) {
#         for(i1 in min(k, m1):max(k - m2, 0L)) {
#             i2 <- k - i1
#             tG <- zeromat
#             if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[, , i1, i2 + 1L])
#             if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[, , i1 + 1L, i2])
#             Gs[, , i1 + 1L, i2 + 1L] <- tG
#             dks[i1 + 1L, i2 + 1L] <- sum(tG) / (2 * k)
#         }
#     }
#     return(dks)
# }

d2_ij_vb <- function(L1, L2, m = 100L, m1 = m, m2 = m, fill_all = !missing(m1) || !missing(m2)) {
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- array(0, dim = c(n, dim(dks)))
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[, i1, i2 + 1L])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[, i1 + 1L, i2])
            Gs[, i1 + 1L, i2 + 1L] <- tG
            dks[i1 + 1L, i2 + 1L] <- sum(tG) / (2 * k)
        }
        if(max(tG) > thr) {
            dks <- dks / 1e10
            Gs <- Gs / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d2_ij_m <- function(A1, A2, m = 100L, m1 = m, m2 = m, fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    # tr <- function(X) sum(diag(X))
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            dks[i1 + 1L, i2 + 1L] <- tr(tG) / (2 * k)
        }
        if(max(tG) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d2_ij_v <- function(L1, L2, m = 100L, m1 = m, m2 = m, fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            dks[i1 + 1L, i2 + 1L] <- sum(tG) / (2 * k)
        }
        if(max(tG) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d2_pj_m <- function(A1, A2, m = 100L, m1 = 1L) {
    if(m1 == 1L) return(d2_1j_m(A1, A2, m))
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, m1 + 1L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    G_k_i <- list()
    G_k_i[1L:(m1 + 1L)] <- list(zeromat)
    for(i in 1L:m1) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        dks[i + 1L, 1L] <- tr(G_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        G_k_i[[1L]] <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        dks[1L, k + 1L] <- tr(G_k_i[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            G_k_i[[i + 1L]] <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            dks[i + 1L, k + 1L] <- tr(G_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(G_k_i[[i + 1L]]) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d2_1j_m <- function(A1, A2, m = 100L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1 # * (dks[1L, 1L] + G_k_0)
    dks[2L, 1L] <- tr(G_k_1) / 2
    for(k in 1L:m) {
        G_k_0 <- A2 %*% (dks[1L, k] * In + G_k_0)
        dks[1L, k + 1L] <- tr(G_k_0) / (2 * k)
        G_k_1 <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        dks[2L, k + 1L] <- tr(G_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d2_pj_v <- function(L1, L2, m = 100L, m1 = 1L) {
    if(m1 == 1L) return(d2_1j_v(L1, L2, m))
    n <- length(L1)
    dks <- matrix(0, m1 + 1L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    G_k_i[1L:(m1 + 1L)] <- list(zeros)
    for(i in 1L:m1) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        G_k_i[[1L]] <- L2 * (dks[1L, k] + G_k_i[[1L]])
        dks[1L, k + 1L] <- sum(G_k_i[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            G_k_i[[i + 1L]] <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(G_k_i[[i + 1L]]) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d2_1j_v <- function(L1, L2, m = 100L) {
    n <- length(L1)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1 # * (dks[1L, 1L] + G_k_0)
    dks[2L, 1L] <- sum(G_k_1) / 2
    for(k in 1L:m) {
        G_k_0 <- L2 * (dks[1L, k] + G_k_0)
        dks[1L, k + 1L] <- sum(G_k_0) / (2 * k)
        G_k_1 <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        dks[2L, k + 1L] <- sum(G_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}


# d2_i1_m <- function(A1, A2, m = 101L) {
#     # tr <- function(X) sum(diag(X))
#     n <- ncol(A1)
#     In <- diag(n)
#     dks <- matrix(0, m + 1L, 2L)
#     dks[1L, 1L] <- 1
#     dks[1L, 2L] <- tr(A2) / 2
#     dks[2L, 1L] <- tr(A1) / 2
#     G_k2_1 <- A2
#     G_k1_0 <- A1
#     for(k in 2L:m) {
#         G_k1_1 <- A1 %*% (dks[k - 1L, 2L] * In + G_k2_1) +
#                   A2 %*% (dks[k, 1L] * In + G_k1_0)
#         dks[k, 2L] <- tr(G_k1_1) / (2 * k)
#         G_k_0 <- A1 %*% (dks[k, 1L] * In + G_k1_0)
#         dks[k + 1L, 1L] <- tr(G_k_0) / (2 * k)
#         G_k2_1 <- G_k1_1
#         G_k1_0 <- G_k_0
#     }
#     return(dks)
# }
#
# d2_i1_v <- function(L1, L2, m = 101L) {
#     n <- length(L1)
#     dks <- matrix(0, m + 1L, 2L)
#     dks[1L, 1L] <- 1
#     dks[1L, 2L] <- sum(L2) / 2
#     dks[2L, 1L] <- sum(L1) / 2
#     G_k2_1 <- L2
#     G_k1_0 <- L1
#     for(k in 2L:m) {
#         G_k1_1 <- L1 * (dks[k - 1L, 2L] + G_k2_1) +
#                   L2 * (dks[k, 1L] + G_k1_0)
#         dks[k, 2L] <- sum(G_k1_1) / (2 * k)
#         G_k_0 <- L1 * (dks[k, 1L] + G_k1_0)
#         dks[k + 1L, 1L] <- sum(G_k_0) / (2 * k)
#         G_k2_1 <- G_k1_1
#         G_k1_0 <- G_k_0
#     }
#     return(dks)
# }

# d2_ip_v <- function(L1, L2, m = 100L, m1 = 1L) {
#     n <- length(L1)
#     dks <- matrix(0, m + 1L, m1 + 1L)
#     dks[1L, 1L] <- 1
#     zeros <- rep.int(0, n)
#     G_k_i <- matrix(0, n, m1 + 1L)
#     for(i in 1L:m1) {
#         G_k_i[, i + 1L] <- L2 * (dks[1L, i] + G_k_i[, i])
#         dks[1L, i + 1L] <- sum(G_k_i[, i + 1L]) / (2 * i)
#     }
#     for(k in 1L:m) {
#         G_k_i[, 1L] <- L1 * (dks[k, 1L] + G_k_i[, 1L])
#         dks[k + 1L, 1L] <- sum(G_k_i[, 1L]) / (2 * k)
#         for(i in 1L:m1) {
#             G_k_i[, i + 1L] <- L1 * (dks[k, i + 1L] + G_k_i[, i + 1L]) +
#                                L2 * (dks[k + 1L, i] + G_k_i[, i])
#             dks[k + 1L, i + 1L] <- sum(G_k_i[, i + 1L]) / (2 * (k + i))
#         }
#     }
#     return(dks)
# }
#
# d2_ip_v <- function(L1, L2, m = 100L, m1 = 1L) {
#     n <- length(L1)
#     dks <- matrix(0, m + 1L, m1 + 1L)
#     dks[1L, 1L] <- 1
#     zeros <- rep.int(0, n)
#     G_k_i <- list()
#     G_k_i[1L:(m1 + 1L)] <- list(zeros)
#     for(i in 1L:m1) {
#         G_k_i[[i + 1L]] <- L2 * (dks[1L, i] + G_k_i[[i]])
#         dks[1L, i + 1L] <- sum(G_k_i[[i + 1L]]) / (2 * i)
#     }
#     for(k in 1L:m) {
#         G_k_i[[1L]] <- L1 * (dks[k, 1L] + G_k_i[[1L]])
#         dks[k + 1L, 1L] <- sum(G_k_i[[1L]]) / (2 * k)
#         for(i in 1L:m1) {
#             G_k_i[[i + 1L]] <- L1 * (dks[k, i + 1L] + G_k_i[[i + 1L]]) +
#                                L2 * (dks[k + 1L, i] + G_k_i[[i]])
#             dks[k + 1L, i + 1L] <- sum(G_k_i[[i + 1L]]) / (2 * (k + i))
#         }
#     }
#     return(dks)
# }


#' This function makes many temporary matrices which are used recursively
#' to calculate polynomials.
#'
d3_ijk_m <- function(A1, A2, A3, m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- tr(tG) / (2 * k)
            }
        }
        if(max(tG) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

d3_ijk_v <- function(L1, L2, L3, m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- sum(tG) / (2 * k)
            }
        }
        if(max(tG) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

d3_pjk_m <- function(A1, A2, A3, m = 100L, m1 = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- array(0, dim = c(m1 + 1L, m + 1L, m + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    Gc <- list()
    Gc[1L:(m1 + 1L)] <- list(zeromat)
    for(i in 1L:m1) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        dks[i + 1L, 1L, 1L] <- tr(Gc[[i + 1L]]) / (2 * i)
    }
    Gn <- list(Gc)
    for(k in 1L:m) {
        Go <- Gn
        Gc[[1L]] <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
        dks[1L, k + 1L, 1L] <- tr(Gc[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            Gc[[i + 1L]] <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                            A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
            dks[i + 1L, k + 1L, 1L] <- tr(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- list(Gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                Gc[[1L]] <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                dks[1L, k - j + 1L, j + 1L] <- tr(Gc[[1L]]) / (2 * k)
                for(i in 1L:m1) {
                    Gc[[i + 1L]] <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                                    A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                                    A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                    dks[i + 1L, k - j + 1L, j + 1L] <- tr(Gc[[i + 1L]]) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
            }
        }
        Gc[[1L]] <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
        dks[1L, 1L, k + 1L] <- tr(Gc[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                            A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
            dks[i + 1L, 1L, k + 1L] <- tr(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        if(max(unlist(Gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

d3_pjk_v <- function(L1, L2, L3, m = 100L, m1 = 1L) {
    n <- length(L1)
    dks <- array(0, dim = c(m1 + 1L, m + 1L, m + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    Gc <- list()
    Gc[1L:(m1 + 1L)] <- list(zeros)
    for(i in 1L:m1) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        dks[i + 1L, 1L, 1L] <- sum(Gc[[i + 1L]]) / (2 * i)
    }
    Gn <- list(Gc)
    for(k in 1L:m) {
        Go <- Gn
        Gc[[1L]] <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
        dks[1L, k + 1L, 1L] <- sum(Gc[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            Gc[[i + 1L]] <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                            L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
            dks[i + 1L, k + 1L, 1L] <- sum(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- list(Gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                Gc[[1L]] <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                dks[1L, k - j + 1L, j + 1L] <- sum(Gc[[1L]]) / (2 * k)
                for(i in 1L:m1) {
                    Gc[[i + 1L]] <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                                    L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                                    L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                    dks[i + 1L, k - j + 1L, j + 1L] <- sum(Gc[[i + 1L]]) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
            }
        }
        Gc[[1L]] <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
        dks[1L, 1L, k + 1L] <- sum(Gc[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            Gc[[i + 1L]] <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                            L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
            dks[i + 1L, 1L, k + 1L] <- sum(Gc[[i + 1L]]) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        if(max(unlist(Gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}


# d3_ij1_m <- function(A1, A2, A3, m = 101L) {
#     n <- ncol(A1)
#     In <- diag(n)
#     zeromat <- matrix(0, n, n)
#     dks <- array(0, dim = c(m + 1L, m + 1L, 2L))
#     dks[1L] <- 1
#     Gs <- array(0, dim = c(n, n, dim(dks)))
#     for(k in 1L:m) {
#         for(i3 in 0L:1L) {
#             for(i1 in (k - i3):0L) {
#                 i2 <- k - i1 - i3
#                 tG <- zeromat
#                 if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[, , i1, i2 + 1L, i3 + 1L])
#                 if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[, , i1 + 1L, i2, i3 + 1L])
#                 if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[, , i1 + 1L, i2 + 1L, i3])
#                 Gs[, , i1 + 1L, i2 + 1L, i3 + 1L] <- tG
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- tr(tG) / (2 * k)
#             }
#         }
#     }
#     return(dks)
# }
#
# d3_ij1_v <- function(A1, A2, A3, m = 101L) {
#     n <- ncol(A1)
#     zeros <- rep.int(0, n)
#     dks <- array(0, dim = c(m + 1L, m + 1L, 2L))
#     dks[1L] <- 1
#     Gs <- array(0, dim = c(n, n, dim(dks)))
#     for(k in 1L:m) {
#         for(i3 in 0L:1L) {
#             for(i1 in (k - i3):0L) {
#                 i2 <- k - i1 - i3
#                 tG <- zeros
#                 if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[, , i1, i2 + 1L, i3 + 1L])
#                 if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[, , i1 + 1L, i2, i3 + 1L])
#                 if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[, , i1 + 1L, i2 + 1L, i3])
#                 Gs[, , i1 + 1L, i2 + 1L, i3 + 1L] <- tG
#                 dks[i1 + 1L, i2 + 1L, i3 + 1L] <- sum(tG) / (2 * k)
#             }
#         }
#     }
#     return(dks)
# }

# h1_i_m <- function(A1, mu = rep.int(0, n), m = 100L) {
#     n <- ncol(A1)
#     In <- diag(n)
#     dks <- dks <- rep.int(c(1, 0), c(1L, m))
#     Go <- matrix(0, n, n)
#     go <- matrix(0, n, 1)
#     for(k in 1L:m) {
#         tG <- A1 %*% (dks[k] * In + Go)
#         tg <- (tG - Go) %*% mu - dks[k] * mu + A1 %*% go
#         dks[k + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
#         Go <- tG
#         go <- tg
#     }
#     return(dks)
# }

h2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m,
                      fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG %*% mu
            if(i1 >= 1L) tg <- tg - Gs[[il2(i1 - 1L, i2)]] %*% mu - dks[i1, i2 + 1L] * mu + A1 %*% gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg - Gs[[il2(i1, i2 - 1L)]] %*% mu - dks[i1 + 1L, i2] * mu + A2 %*% gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

h2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m,
                      fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG * mu
            if(i1 >= 1L) tg <- tg - (Gs[[il2(i1 - 1L, i2)]] + dks[i1, i2 + 1L]) * mu + L1 * gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg - (Gs[[il2(i1, i2 - 1L)]] + dks[i1 + 1L, i2]) * mu + L2 * gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m,
                         fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG %*% mu
            if(i1 >= 1L) tg <- tg + A1 %*% gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg - Gs[[il2(i1, i2 - 1L)]] %*% mu - dks[i1 + 1L, i2] * mu + A2 %*% gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m,
                         fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG * mu
            if(i1 >= 1L) tg <- tg + L1 * gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg - (Gs[[il2(i1, i2 - 1L)]] + dks[i1 + 1L, i2]) * mu + L2 * gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil2_pj_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    if(m1 == 1L) return(htil2_1j_m(A1, A2, mu, m))
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, m1 + 1L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(m1 + 1L)] <- list(zeromat)
    g_k_i[1L:(m1 + 1L)] <- list(zerovec)
    for(i in 1L:m1) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu + A1 %*% g_k_i[[i]]
        dks[i + 1L, 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG - G_k_i[[1L]] - (dks[1L, k] * In)) %*% mu + A2 %*% g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- (tr(G_k_i[[1L]]) + c(crossprod(mu, g_k_i[[1L]]))) / (2 * k)
        for(i in 1L:m1) {
            tG <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG - G_k_i[[i + 1L]]
                                - (dks[i + 1L, k] * In)) %*% mu +
                               A1 %*% g_k_i[[i]] + A2 %*% g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil2_1j_m <- function(A1, A2, mu = rep.int(0, n), m = 100L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1 # * (dks[1L, 1L] + G_k_0)
    g_k_0 <- matrix(0, n, 1)
    g_k_1 <- G_k_1 %*% mu # + A1 %*% g_k_0
    dks[2L, 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / 2
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_0)
        g_k_0 <- (tG - G_k_0 - (dks[1L, k] * In)) %*% mu + A2 %*% g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- (tr(G_k_0) + c(crossprod(mu, g_k_0))) / (2 * k)
        tG <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        g_k_1 <- (tG - G_k_1 - (dks[2L, k] * In)) %*% mu +
                 A1 %*% g_k_0 + A2 %*% g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil2_pj_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    if(m1 == 1L) return(htil2_1j_v(L1, L2, mu, m))
    n <- length(L1)
    dks <- matrix(0, m1 + 1L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(m1 + 1L)] <- list(zeros)
    g_k_i[1L:(m1 + 1L)] <- list(zeros)
    for(i in 1L:m1) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu + L1 * g_k_i[[i]]
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG - G_k_i[[1L]] - dks[1L, k]) * mu + L2 * g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- sum(G_k_i[[1L]] + mu * g_k_i[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            tG <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG - G_k_i[[i + 1L]] - dks[i + 1L, k]) * mu +
                               L1 * g_k_i[[i]] + L2 * g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil2_1j_v <- function(L1, L2, mu = rep.int(0, n), m = 100L) {
    n <- length(L1)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1 # * (dks[1L, 1L] + G_k_0)
    g_k_0 <- rep.int(0, n)
    g_k_1 <- G_k_1 * mu # + L1 %*% g_k_0
    dks[2L, 1L] <- sum(G_k_1 + mu * g_k_1) / 2
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_0)
        g_k_0 <- (tG - G_k_0 - dks[1L, k]) * mu + L2 * g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- sum(G_k_0 + mu * g_k_0) / (2 * k)
        tG <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        g_k_1 <- (tG - G_k_1 - dks[2L, k]) * mu + L1 * g_k_0 + L2 * g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- sum(G_k_1 + mu * g_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}


h3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG %*% mu
                if(i1 >= 1L) tg <- tg - Gs[[il3(i1 - 1L, i2, i3)]] %*% mu - dks[i1, i2 + 1L, i3 + 1L] * mu + A1 %*% gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg - Gs[[il3(i1, i2 - 1L, i3)]] %*% mu - dks[i1 + 1L, i2, i3 + 1L] * mu + A2 %*% gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg - Gs[[il3(i1, i2, i3 - 1L)]] %*% mu - dks[i1 + 1L, i2 + 1L, i3] * mu + A3 %*% gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
            }
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

h3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG * mu
                if(i1 >= 1L) tg <- tg - (Gs[[il3(i1 - 1L, i2, i3)]] + dks[i1, i2 + 1L, i3 + 1L]) * mu + L1 * gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg - (Gs[[il3(i1, i2 - 1L, i3)]] + dks[i1 + 1L, i2, i3 + 1L]) * mu + L2 * gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg - (Gs[[il3(i1, i2, i3 - 1L)]] + dks[i1 + 1L, i2 + 1L, i3]) * mu + L3 * gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
            }
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

htil3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG %*% mu
                if(i1 >= 1L) tg <- tg + A1 %*% gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg - Gs[[il3(i1, i2 - 1L, i3)]] %*% mu - dks[i1 + 1L, i2, i3 + 1L] * mu + A2 %*% gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg - Gs[[il3(i1, i2, i3 - 1L)]] %*% mu - dks[i1 + 1L, i2 + 1L, i3] * mu + A3 %*% gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
            }
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

htil3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG * mu
                if(i1 >= 1L) tg <- tg + L1 * gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg - (Gs[[il3(i1, i2 - 1L, i3)]] + dks[i1 + 1L, i2, i3 + 1L]) * mu + L2 * gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg - (Gs[[il3(i1, i2, i3 - 1L)]] + dks[i1 + 1L, i2 + 1L, i3]) * mu + L3 * gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
            }
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

htil3_pjk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- array(0, dim = c(m1 + 1L, m + 1L, m + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    Gc <- list()
    gc <- list()
    Gc[1L:(m1 + 1L)] <- list(zeromat)
    gc[1L:(m1 + 1L)] <- list(zerovec)
    for(i in 1L:m1) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu + A1 %*% gc[[i]]
        dks[i + 1L, 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
        gc[[1L]] <- (tG - Go[[1L]][[1L]] - (dks[1L, k, 1L] * In)) %*% mu + A2 %*% go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:m1) {
            tG <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                  A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG - Go[[1L]][[i + 1L]]
                                - (dks[i + 1L, k, 1L] * In)) %*% mu +
                               A1 %*% gc[[i]] + A2 %*% go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                gc[[1L]] <- (tG - Go[[j + 1L]][[1L]] - Go[[j]][[1L]] - ((dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j]) * In)) %*% mu + A2 %*% go[[j + 1L]][[1L]] + A3 %*% go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
                for(i in 1L:m1) {
                    tG <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                          A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                          A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG - Go[[j + 1L]][[i + 1L]] - Go[[j]][[i + 1L]]
                                     - ((dks[i + 1L, k - j, j + 1L] + dks[i + 1L, k - j + 1L, j]) * In)) %*% mu +
                                    A1 %*% gc[[i]] + A2 %*% go[[j + 1L]][[i + 1L]] + A3 %*% go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
        gc[[1L]] <- (tG - Go[[k]][[1L]] - (dks[1L, 1L, k] * In)) %*% mu + A3 %*% go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:m1) {
            tG <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                  A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG - Go[[k]][[i + 1L]]
                             - (dks[i + 1L, 1L, k] * In)) %*% mu +
                            A1 %*% gc[[i]] + A3 %*% go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gc)) > thr || max(unlist(gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

htil3_pjk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    n <- length(L1)
    dks <- array(0, dim = c(m1 + 1L, m + 1L, m + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    Gc <- list()
    gc <- list()
    Gc[1L:(m1 + 1L)] <- list(zeros)
    gc[1L:(m1 + 1L)] <- list(zeros)
    for(i in 1L:m1) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] * mu + L1 * gc[[i]]
        dks[i + 1L, 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
        gc[[1L]] <- (tG - (Go[[1L]][[1L]] + dks[1L, k, 1L])) * mu + L2 * go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:m1) {
            tG <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                  L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG - (Go[[1L]][[i + 1L]]
                                + dks[i + 1L, k, 1L])) * mu +
                               L1 * gc[[i]] + L2 * go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                gc[[1L]] <- (tG - (Go[[j + 1L]][[1L]] + Go[[j]][[1L]] + dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j])) * mu + L2 * go[[j + 1L]][[1L]] + L3 * go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
                for(i in 1L:m1) {
                    tG <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                          L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                          L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG - (Go[[j + 1L]][[i + 1L]] + Go[[j]][[i + 1L]]
                                        + dks[i + 1L, k - j, j + 1L]
                                        + dks[i + 1L, k - j + 1L, j])) * mu +
                                       L1 * gc[[i]] + L2 * go[[j + 1L]][[i + 1L]] + L3 * go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
        gc[[1L]] <- (tG - (Go[[k]][[1L]] + dks[1L, 1L, k])) * mu + L3 * go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:m1) {
            tG <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                  L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG - (Go[[k]][[i + 1L]]
                             + dks[i + 1L, 1L, k])) * mu +
                            L1 * gc[[i]] + L3 * go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        # if(!isFALSE(max(unlist(Gc)) > thr || max(unlist(gc)) > thr)) browser()
        if(max(unlist(Gc)) > thr || max(unlist(gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}



#' Recursion for hhat
#'
#' \code{hhat2_ij_m()} and \code{hhat2_ij_m()} are recursive algorithms
#' for truncation error in noncentral case.
#'
#' These provide recursive algorithms for \eqn{\hat{h}_{i,j}} in
#' Hillier et al. (2014, theorem 7).  This recursion is said to be very similar
#' to those for \eqn{h_{i,j}} in the note therein, but differs in the signs
#' of some terms.
#'
hhat2_ij_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeromat
            if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L] * In + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2] * In + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG %*% mu
            if(i1 >= 1L) tg <- tg + A1 %*% gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg + Gs[[il2(i1, i2 - 1L)]] %*% mu + dks[i1 + 1L, i2] * mu + A2 %*% gs[[il2(i1, i2 - 1L)]]
            gs[[il2(i1, i2)]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

hhat2_ij_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, fill_all = !missing(m1) || !missing(m2)) {
    il2 <- function(i1, i2) i1 + i2 * (m1 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- matrix(rep.int(c(1, 0), c(1L, (m1 + 1) * (m2 + 1) - 1L)), m1 + 1, m2 + 1)
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_mat <- outer(0:m1, 0:m2, "+")
    kmax <- if(fill_all) sum(dim(dks) - 1L) else max(dim(dks) - 1L)
    for(k in 1L:kmax) {
        Gs[which(order_mat == k - 2L)] <- 0
        gs[which(order_mat == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2, 0L)) {
            i2 <- k - i1
            tG <- zeros
            if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L] + Gs[[il2(i1 - 1L, i2)]])
            if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2] + Gs[[il2(i1, i2 - 1L)]])
            Gs[[il2(i1, i2)]] <- tG
            tg <- tG * mu
            if(i1 >= 1L) tg <- tg + L1 * gs[[il2(i1 - 1L, i2)]]
            if(i2 >= 1L) tg <- tg + (Gs[[il2(i1, i2 - 1L)]] + dks[i1 + 1L, i2]) * mu + L2 * gs[[il2(i1, i2 - 1L)]]
            gs[[sum(c(i1, i2) * c(1L, m1 + 1L)) + 1L]] <- tg
            dks[i1 + 1L, i2 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

hhat2_pj_m <- function(A1, A2, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    if(m1 == 1L) return(hhat2_1j_m(A1, A2, mu, m))
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, m1 + 1L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(m1 + 1L)] <- list(zeromat)
    g_k_i[1L:(m1 + 1L)] <- list(zerovec)
    for(i in 1L:m1) {
        G_k_i[[i + 1L]] <- A1 %*% (dks[i, 1L] * In + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] %*% mu + A1 %*% g_k_i[[i]]
        dks[i + 1L, 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG + G_k_i[[1L]] + (dks[1L, k] * In)) %*% mu + A2 %*% g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- (tr(G_k_i[[1L]]) + c(crossprod(mu, g_k_i[[1L]]))) / (2 * k)
        for(i in 1L:m1) {
            tG <- A1 %*% (dks[i, k + 1L] * In + G_k_i[[i]]) +
                               A2 %*% (dks[i + 1L, k] * In + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG + G_k_i[[i + 1L]]
                                + (dks[i + 1L, k] * In)) %*% mu +
                               A1 %*% g_k_i[[i]] + A2 %*% g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- (tr(G_k_i[[i + 1L]]) + c(crossprod(mu, g_k_i[[i + 1L]]))) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

hhat2_1j_m <- function(A1, A2, mu = rep.int(0, n), m = 100L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- matrix(0, n, n)
    G_k_1 <- A1 # * (dks[1L, 1L] + G_k_0)
    g_k_0 <- matrix(0, n, 1)
    g_k_1 <- G_k_1 %*% mu # + A1 %*% g_k_0
    dks[2L, 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / 2
    for(k in 1L:m) {
        tG <- A2 %*% (dks[1L, k] * In + G_k_0)
        g_k_0 <- (tG + G_k_0 + (dks[1L, k] * In)) %*% mu + A2 %*% g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- (tr(G_k_0) + c(crossprod(mu, g_k_0))) / (2 * k)
        tG <- A1 %*% (dks[1L, k + 1L] * In + G_k_0) +
                 A2 %*% (dks[2L, k] * In + G_k_1)
        g_k_1 <- (tG + G_k_1 + (dks[2L, k] * In)) %*% mu +
                 A1 %*% g_k_0 + A2 %*% g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- (tr(G_k_1) + c(crossprod(mu, g_k_1))) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

hhat2_pj_v <- function(L1, L2, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    if(m1 == 1L) return(hhat2_1j_v(L1, L2, mu, m))
    n <- length(L1)
    dks <- matrix(0, m1 + 1L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    G_k_i <- list()
    g_k_i <- list()
    G_k_i[1L:(m1 + 1L)] <- list(zeros)
    g_k_i[1L:(m1 + 1L)] <- list(zeros)
    for(i in 1L:m1) {
        G_k_i[[i + 1L]] <- L1 * (dks[i, 1L] + G_k_i[[i]])
        g_k_i[[i + 1L]] <- G_k_i[[i + 1L]] * mu + L1 * g_k_i[[i]]
        dks[i + 1L, 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * i)
    }
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_i[[1L]])
        g_k_i[[1L]] <- (tG + G_k_i[[1L]] + dks[1L, k]) * mu + L2 * g_k_i[[1L]]
        G_k_i[[1L]] <- tG
        dks[1L, k + 1L] <- sum(G_k_i[[1L]] + mu * g_k_i[[1L]]) / (2 * k)
        for(i in 1L:m1) {
            tG <- L1 * (dks[i, k + 1L] + G_k_i[[i]]) +
                               L2 * (dks[i + 1L, k] + G_k_i[[i + 1L]])
            g_k_i[[i + 1L]] <- (tG + G_k_i[[i + 1L]] + dks[i + 1L, k]) * mu +
                               L1 * g_k_i[[i]] + L2 * g_k_i[[i + 1L]]
            G_k_i[[i + 1L]] <- tG
            dks[i + 1L, k + 1L] <- sum(G_k_i[[i + 1L]] + mu * g_k_i[[i + 1L]]) / (2 * (k + i))
        }
        if(max(unlist(G_k_i)) > thr || max(unlist(g_k_i)) > thr) {
            dks <- dks / 1e10
            G_k_i <- lapply(G_k_i, function(x) x / 1e10)
            g_k_i <- lapply(g_k_i, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

hhat2_1j_v <- function(L1, L2, mu = rep.int(0, n), m = 100L) {
    n <- length(L1)
    dks <- matrix(0, 2L, m + 1L)
    dks[1L, 1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    G_k_0 <- rep.int(0, n)
    G_k_1 <- L1 # * (dks[1L, 1L] + G_k_0)
    g_k_0 <- rep.int(0, n)
    g_k_1 <- G_k_1 * mu # + L1 %*% g_k_0
    dks[2L, 1L] <- sum(G_k_1 + mu * g_k_1) / 2
    for(k in 1L:m) {
        tG <- L2 * (dks[1L, k] + G_k_0)
        g_k_0 <- (tG + G_k_0 + dks[1L, k]) * mu + L2 * g_k_0
        G_k_0 <- tG
        dks[1L, k + 1L] <- sum(G_k_0 + mu * g_k_0) / (2 * k)
        tG <- L1 * (dks[1L, k + 1L] + G_k_0) + L2 * (dks[2L, k] + G_k_1)
        g_k_1 <- (tG + G_k_1 + dks[2L, k]) * mu + L1 * g_k_0 + L2 * g_k_1
        G_k_1 <- tG
        dks[2L, k + 1L] <- sum(G_k_1 + mu * g_k_1) / (2 * (k + 1))
        if(max(G_k_1) > thr || max(g_k_1) > thr) {
            dks <- dks / 1e10
            G_k_0 <- G_k_0 / 1e10
            G_k_1 <- G_k_1 / 1e10
            g_k_0 <- g_k_0 / 1e10
            g_k_1 <- g_k_1 / 1e10
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}



hhat3_ijk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- ncol(A1)
    In <- diag(n)
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeromat)
    gs <- list(zerovec)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeromat
                if(i1 >= 1L) tG <- tG + A1 %*% (dks[i1, i2 + 1L, i3 + 1L] * In + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + A2 %*% (dks[i1 + 1L, i2, i3 + 1L] * In + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + A3 %*% (dks[i1 + 1L, i2 + 1L, i3] * In + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG %*% mu
                if(i1 >= 1L) tg <- tg + A1 %*% gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg + Gs[[il3(i1, i2 - 1L, i3)]] %*% mu + dks[i1 + 1L, i2, i3 + 1L] * mu + A2 %*% gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg + Gs[[il3(i1, i2, i3 - 1L)]] %*% mu + dks[i1 + 1L, i2 + 1L, i3] * mu + A3 %*% gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (tr(tG) + c(crossprod(mu, tg))) / (2 * k)
            }
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

hhat3_ijk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, m1 = m, m2 = m, m3 = m,
                 fill_across = c(!missing(m1), !missing(m2), !missing(m3))) { # , verbose = m > 200L) {
    il3 <- function(i1, i2, i3) i1 + i2 * (m1 + 1L) + i3 * (m1 + 1L) * (m2 + 1L) + 1L
    n <- length(L1)
    zeros <- rep.int(0, n)
    dks <- array(0, dim = c(m1 + 1L, m2 + 1L, m3 + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    Gs <- list(zeros)
    gs <- list(zeros)
    order_array <- outer(outer(0:m1, 0:m2, "+"), 0:m3, "+")
    kmax <- min(any(!fill_across) * m + sum(fill_across * c(m1, m2, m3)), sum(dim(dks) - 1L))
    for(k in 1L:kmax) {
        Gs[which(order_array == k - 2L)] <- 0
        gs[which(order_array == k - 2L)] <- 0
        for(i1 in min(k, m1):max(k - m2 - m3, 0L, fill_across[1L] * (k + m1 - kmax))) {
            for(i2 in min(k - i1, m2):max(k - i1 - m3, 0L)) {
                i3 <- k - i1 - i2
                if(fill_across[3L] && i1 + i2 + m3 > kmax) next
                if(fill_across[2L] && i1 + i3 + m2 > kmax) next
                # if(fill_across[1] && i2 + i3 + m1 > kmax) next
                tG <- zeros
                if(i1 >= 1L) tG <- tG + L1 * (dks[i1, i2 + 1L, i3 + 1L] + Gs[[il3(i1 - 1L, i2, i3)]])
                if(i2 >= 1L) tG <- tG + L2 * (dks[i1 + 1L, i2, i3 + 1L] + Gs[[il3(i1, i2 - 1L, i3)]])
                if(i3 >= 1L) tG <- tG + L3 * (dks[i1 + 1L, i2 + 1L, i3] + Gs[[il3(i1, i2, i3 - 1L)]])
                Gs[[il3(i1, i2, i3)]] <- tG
                tg <- tG * mu
                if(i1 >= 1L) tg <- tg + L1 * gs[[il3(i1 - 1L, i2, i3)]]
                if(i2 >= 1L) tg <- tg + (Gs[[il3(i1, i2 - 1L, i3)]] + dks[i1 + 1L, i2, i3 + 1L]) * mu + L2 * gs[[il3(i1, i2 - 1L, i3)]]
                if(i3 >= 1L) tg <- tg + (Gs[[il3(i1, i2, i3 - 1L)]] + dks[i1 + 1L, i2 + 1L, i3]) * mu + L3 * gs[[il3(i1, i2, i3 - 1L)]]
                gs[[il3(i1, i2, i3)]] <- tg
                dks[i1 + 1L, i2 + 1L, i3 + 1L] <- (sum(tG) + sum(mu * tg)) / (2 * k)
            }
        }
        if(max(tG) > thr || max(tg) > thr) {
            dks <- dks / 1e10
            Gs <- lapply(Gs, function(x) x / 1e10)
            gs <- lapply(gs, function(x) x / 1e10)
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
        # if((k %% 10 == 0) && verbose) cat("  k =", k, "done\n")
    }
    return(dks)
}

hhat3_pjk_m <- function(A1, A2, A3, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    n <- ncol(A1)
    In <- diag(n)
    dks <- array(0, dim = c(m1 + 1L, m + 1L, m + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeromat <- matrix(0, n, n)
    zerovec <- matrix(0, n, 1)
    Gc <- list()
    gc <- list()
    Gc[1L:(m1 + 1L)] <- list(zeromat)
    gc[1L:(m1 + 1L)] <- list(zerovec)
    for(i in 1L:m1) {
        Gc[[i + 1L]] <- A1 %*% (dks[i, 1L, 1L] * In + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] %*% mu + A1 %*% gc[[i]]
        dks[i + 1L, 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- A2 %*% (dks[1L, k, 1L] * In + Go[[1L]][[1L]])
        gc[[1L]] <- (tG + Go[[1L]][[1L]] + (dks[1L, k, 1L] * In)) %*% mu + A2 %*% go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:m1) {
            tG <- A1 %*% (dks[i, k + 1L, 1L] * In + Gc[[i]]) +
                  A2 %*% (dks[i + 1L, k, 1L] * In + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG + Go[[1L]][[i + 1L]]
                                + (dks[i + 1L, k, 1L] * In)) %*% mu +
                               A1 %*% gc[[i]] + A2 %*% go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- A2 %*% (dks[1L, k - j, j + 1L] * In + Go[[j + 1L]][[1L]]) + A3 %*% (dks[1L, k - j + 1L, j] * In + Go[[j]][[1L]])
                gc[[1L]] <- (tG + Go[[j + 1L]][[1L]] + Go[[j]][[1L]] + ((dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j]) * In)) %*% mu + A2 %*% go[[j + 1L]][[1L]] + A3 %*% go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
                for(i in 1L:m1) {
                    tG <- A1 %*% (dks[i, k - j + 1L, j + 1L] * In + Gc[[i]]) +
                          A2 %*% (dks[i + 1L, k - j, j + 1L] * In + Go[[j + 1L]][[i + 1L]]) +
                          A3 %*% (dks[i + 1L, k - j + 1L, j] * In + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG + Go[[j + 1L]][[i + 1L]] + Go[[j]][[i + 1L]]
                                     + ((dks[i + 1L, k - j, j + 1L] + dks[i + 1L, k - j + 1L, j]) * In)) %*% mu +
                                    A1 %*% gc[[i]] + A2 %*% go[[j + 1L]][[i + 1L]] + A3 %*% go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- A3 %*% (dks[1L, 1L, k] * In + Go[[k]][[1L]])
        gc[[1L]] <- (tG + Go[[k]][[1L]] + (dks[1L, 1L, k] * In)) %*% mu + A3 %*% go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (tr(Gc[[1L]]) + c(crossprod(mu, gc[[1L]]))) / (2 * k)
        for(i in 1L:m1) {
            tG <- A1 %*% (dks[i, 1L, k + 1L] * In + Gc[[i]]) +
                  A3 %*% (dks[i + 1L, 1L, k] * In + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG + Go[[k]][[i + 1L]]
                             + (dks[i + 1L, 1L, k] * In)) %*% mu +
                            A1 %*% gc[[i]] + A3 %*% go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (tr(Gc[[i + 1L]]) + c(crossprod(mu, gc[[i + 1L]]))) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gc)) > thr || max(unlist(gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}

hhat3_pjk_v <- function(L1, L2, L3, mu = rep.int(0, n), m = 100L, m1 = 1L) {
    n <- length(L1)
    dks <- array(0, dim = c(m1 + 1L, m + 1L, m + 1L))
    dks[1L] <- 1
    attr(dks, "scale") <- 1
    thr <- .Machine$double.xmax / 100 / n
    zeros <- rep.int(0, n)
    Gc <- list()
    gc <- list()
    Gc[1L:(m1 + 1L)] <- list(zeros)
    gc[1L:(m1 + 1L)] <- list(zeros)
    for(i in 1L:m1) {
        Gc[[i + 1L]] <- L1 * (dks[i, 1L, 1L] + Gc[[i]])
        gc[[i + 1L]] <- Gc[[i + 1L]] * mu + L1 * gc[[i]]
        dks[i + 1L, 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * i)
    }
    Gn <- list(Gc)
    gn <- list(gc)
    for(k in 1L:m) {
        Go <- Gn
        go <- gn
        tG <- L2 * (dks[1L, k, 1L] + Go[[1L]][[1L]])
        gc[[1L]] <- (tG + (Go[[1L]][[1L]] + dks[1L, k, 1L])) * mu + L2 * go[[1L]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, k + 1L, 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:m1) {
            tG <- L1 * (dks[i, k + 1L, 1L] + Gc[[i]]) +
                  L2 * (dks[i + 1L, k, 1L] + Go[[1L]][[i + 1L]])
            gc[[i + 1L]] <- (tG + (Go[[1L]][[i + 1L]]
                                + dks[i + 1L, k, 1L])) * mu +
                               L1 * gc[[i]] + L2 * go[[1L]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, k + 1L, 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- list(Gc)
        gn <- list(gc)
        if(k >= 2L) {
            for(j in 1L:(k - 1L)) {
                tG <- L2 * (dks[1L, k - j, j + 1L] + Go[[j + 1L]][[1L]]) + L3 * (dks[1L, k - j + 1L, j] + Go[[j]][[1L]])
                gc[[1L]] <- (tG + (Go[[j + 1L]][[1L]] + Go[[j]][[1L]] + dks[1L, k - j, j + 1L] + dks[1L, k - j + 1L, j])) * mu + L2 * go[[j + 1L]][[1L]] + L3 * go[[j]][[1L]]
                Gc[[1L]] <- tG
                dks[1L, k - j + 1L, j + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
                for(i in 1L:m1) {
                    tG <- L1 * (dks[i, k - j + 1L, j + 1L] + Gc[[i]]) +
                          L2 * (dks[i + 1L, k - j, j + 1L] + Go[[j + 1L]][[i + 1L]]) +
                          L3 * (dks[i + 1L, k - j + 1L, j] + Go[[j]][[i + 1L]])
                    gc[[i + 1L]] <- (tG + (Go[[j + 1L]][[i + 1L]] + Go[[j]][[i + 1L]]
                                        + dks[i + 1L, k - j, j + 1L]
                                        + dks[i + 1L, k - j + 1L, j])) * mu +
                                       L1 * gc[[i]] + L2 * go[[j + 1L]][[i + 1L]] + L3 * go[[j]][[i + 1L]]
                    Gc[[i + 1L]] <- tG
                    dks[i + 1L, k - j + 1L, j + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
                }
                Gn <- c(Gn, list(Gc))
                gn <- c(gn, list(gc))
            }
        }
        tG <- L3 * (dks[1L, 1L, k] + Go[[k]][[1L]])
        gc[[1L]] <- (tG + (Go[[k]][[1L]] + dks[1L, 1L, k])) * mu + L3 * go[[k]][[1L]]
        Gc[[1L]] <- tG
        dks[1L, 1L, k + 1L] <- (sum(Gc[[1L]]) + sum(mu * gc[[1L]])) / (2 * k)
        for(i in 1L:m1) {
            tG <- L1 * (dks[i, 1L, k + 1L] + Gc[[i]]) +
                  L3 * (dks[i + 1L, 1L, k] + Go[[k]][[i + 1L]])
            gc[[i + 1L]] <- (tG + (Go[[k]][[i + 1L]]
                             + dks[i + 1L, 1L, k])) * mu +
                            L1 * gc[[i]] + L3 * go[[k]][[i + 1L]]
            Gc[[i + 1L]] <- tG
            dks[i + 1L, 1L, k + 1L] <- (sum(Gc[[i + 1L]]) + sum(mu * gc[[i + 1L]])) / (2 * (k + i))
        }
        Gn <- c(Gn, list(Gc))
        gn <- c(gn, list(gc))
        if(max(unlist(Gc)) > thr || max(unlist(gc)) > thr) {
            dks <- dks / 1e10
            Gn <- lapply(Gn, function(x) lapply(x, function(y) y / 1e10))
            gn <- lapply(gn, function(x) lapply(x, function(y) y / 1e10))
            attr(dks, "scale") <- attr(dks, "scale") / 1e10
        }
    }
    return(dks)
}
