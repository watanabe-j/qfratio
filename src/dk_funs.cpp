#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <limits>
#include "dk_funs.h"

#ifdef _OPENMP
    #include <omp.h>
    // [[Rcpp::plugins(openmp)]]
#endif

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// // [[Rcpp::export]]
Eigen::ArrayXd d1_i_vE(const Eigen::ArrayXd& L, const int m, Eigen::ArrayXd& lscf) {
    int n = L.size();
    ArrayXd dks = ArrayXd::Zero(m + 1);
    dks(0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd uk = ArrayXd::Zero(n);
    for(int k = 1; k <= m; k++) {
        uk = L * (dks(k - 1) + uk);
        dks(k) = uk.sum() / (2 * k);
        if(uk.maxCoeff() > thr) {
            dks(k) /= 1e10;
            uk /= 1e10;
            lscf.tail(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXd d1_i_mE(const Eigen::MatrixXd& A, const int m, Eigen::ArrayXd& lscf) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A, Eigen::EigenvaluesOnly);
    ArrayXd L = eigA.eigenvalues();
    return d1_i_vE(L, m, lscf);
}

// // [[Rcpp::export]]
Eigen::ArrayXd dtil1_i_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& mu,
    const int m, Eigen::ArrayXd& lscf) {
    int n = L.size();
    ArrayXd D = square(mu);
    ArrayXd dks = ArrayXd::Zero(m + 1);
    dks(0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd uk = ArrayXd::Zero(n);
    ArrayXd vk = ArrayXd::Zero(n);
    for(int k = 1; k <= m; k++) {
        uk = L * (dks(k - 1) + uk);
        vk = D * uk + L * vk;
        dks(k) = (uk + vk).sum() / (2 * k);
        if(uk.maxCoeff() > thr || vk.maxCoeff() > thr) {
            dks(k) /= 1e10;
            uk /= 1e10;
            vk /= 1e10;
            lscf.tail(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXd dtil1_i_mE(const Eigen::MatrixXd& A, const Eigen::VectorXd& mu,
        const int m, Eigen::ArrayXd& lscf) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
    ArrayXd L = eigA.eigenvalues();
    ArrayXd mud = eigA.eigenvectors().transpose() * mu;
    return dtil1_i_vE(L, mud, m, lscf);
}

// // [[Rcpp::export]]
Eigen::ArrayXXd arl_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& D, const int m) {
    const int n = L.size();
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXXd arls = ArrayXXd::Zero(m + 1, m + 1);
    arls.col(0) = d1_i_vE(L, m, lscf);
    ArrayXXd wrls = ArrayXXd::Zero(n, m);
    for(int k = 0; k < m; k++) {
        for(int l = 0; l <= k; l++) {
            wrls.col(l) = L * (arls(k, l) + wrls.col(l));
            arls(k + 1, l + 1) = (D * wrls.col(l)).sum();
        }
    }
    return arls;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd arl_mE(const Eigen::MatrixXd& A, const Eigen::VectorXd& mu, const int m) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
    ArrayXd L = eigA.eigenvalues();
    ArrayXd mud = eigA.eigenvectors().transpose() * mu;
    ArrayXd D = square(mud);
    return arl_vE(L, D, m);
}


// // [[Rcpp::export]]
Eigen::ArrayXXd d2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const int m, const int p, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd G_k_i = MatrixXd::Zero(n, n * (p + 1));
    for(int i = 1; i <= p; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        dks(i, 0) = G_k_i.block(0, i * n, n, n).trace() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        G_k_i.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        dks(0, k) = G_k_i.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= p; i++) {
            G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                                          A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            dks(i, k) = G_k_i.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const int m, const int p, Eigen::ArrayXXd& lscf) {
    const int n = A1.size();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXXd G_k_i = ArrayXXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        dks(i, 0) = G_k_i.col(i).sum() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        G_k_i.col(0) = A2 * (dks(0, k - 1) + G_k_i.col(0));
        dks(0, k) = G_k_i.col(0).sum() / (2 * k);
        for(int i = 1; i <= p; i++) {
            G_k_i.col(i) = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                           A2 * (dks(i, k - 1) + G_k_i.col(i));
            dks(i, k) = G_k_i.col(i).sum() / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}


// // [[Rcpp::export]]
Eigen::ArrayXXd d2_ij_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const int m, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd G_k_i = MatrixXd::Zero(n, n * (m + 1));
    for(int i = 1; i <= m; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        dks(i, 0) = G_k_i.block(0, i * n, n, n).trace() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        G_k_i.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        dks(0, k) = G_k_i.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                                          A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            dks(i, k) = G_k_i.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d2_ij_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const int m, Eigen::ArrayXXd& lscf) {
    const int n = A1.size();
    ArrayXXd dks = ArrayXXd::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd tG(n);
    ArrayXXd G_k_i = ArrayXXd::Zero(n, m + 1);
    for(int i = 1; i <= m; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        dks(i, 0) = G_k_i.col(i).sum() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        G_k_i.col(0) = A2 * (dks(0, k - 1) + G_k_i.col(0));
        dks(0, k) = G_k_i.col(0).sum() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            G_k_i.col(i) = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                           A2 * (dks(i, k - 1) + G_k_i.col(i));
            dks(i, k) = G_k_i.col(i).sum() / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h2_ij_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd G_k_i = MatrixXd::Zero(n, n * (m + 1));
    MatrixXd g_k_i = MatrixXd::Zero(n, m + 1);
    for(int i = 1; i <= m; i++) {
        tG = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        g_k_i.col(i) = (tG - G_k_i.block(0, (i - 1) * n, n, n) - (dks(i - 1, 0) * In)) * mu + A1 * g_k_i.col(i - 1);
        G_k_i.block(0, i * n, n, n) = tG;
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        tG = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = (tG - G_k_i.block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * g_k_i.col(0);
        G_k_i.block(0, 0, n, n) = tG;
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            g_k_i.col(i) = (tG - G_k_i.block(0, (i - 1) * n, n, n) - G_k_i.block(0, i * n, n, n)
                            - ((dks(i - 1, k) + dks(i, k - 1)) * In)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.block(0, i * n, n, n) = tG;
            dks(i, k) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h2_ij_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int m, Eigen::ArrayXXd& lscf) {
    const int n = A1.size();
    ArrayXXd dks = ArrayXXd::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd tG(n);
    ArrayXXd G_k_i = ArrayXXd::Zero(n, m + 1);
    ArrayXXd g_k_i = ArrayXXd::Zero(n, m + 1);
    for(int i = 1; i <= m; i++) {
        tG = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = (tG - G_k_i.col(i - 1) - dks(i - 1, 0)) * mu + A1 * g_k_i.col(i - 1);
        G_k_i.col(i) = tG;
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        tG = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = (tG - G_k_i.col(0) - dks(0, k - 1)) * mu + A2 * g_k_i.col(0);
        G_k_i.col(0) = tG;
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                 A2 * (dks(i, k - 1) + G_k_i.col(i));
            g_k_i.col(i) = (tG - G_k_i.col(i - 1) - G_k_i.col(i)
                            - dks(i - 1, k) - dks(i, k - 1)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.col(i) = tG;
            dks(i, k) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd G_k_i = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd g_k_i = MatrixXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        g_k_i.col(i) = G_k_i.block(0, i * n, n, n) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        tG = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = (tG - G_k_i.block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * g_k_i.col(0);
        G_k_i.block(0, 0, n, n) = tG;
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            g_k_i.col(i) = (tG - G_k_i.block(0, i * n, n, n)
                            - (dks(i, k - 1) * In)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.block(0, i * n, n, n) = tG;
            dks(i, k) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf) {
    const int n = A1.size();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd tG(n);
    ArrayXXd G_k_i = ArrayXXd::Zero(n, p + 1);
    ArrayXXd g_k_i = ArrayXXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = G_k_i.col(i) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        tG = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = (tG - G_k_i.col(0) - dks(0, k - 1)) * mu + A2 * g_k_i.col(0);
        G_k_i.col(0) = tG;
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                 A2 * (dks(i, k - 1) + G_k_i.col(i));
            g_k_i.col(i) = (tG - G_k_i.col(i) - dks(i, k - 1)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.col(i) = tG;
            dks(i, k) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd G_k_i = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd g_k_i = MatrixXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        g_k_i.col(i) = G_k_i.block(0, i * n, n, n) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        tG = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = (tG + G_k_i.block(0, 0, n, n) + (dks(0, k - 1) * In)) * mu + A2 * g_k_i.col(0);
        G_k_i.block(0, 0, n, n) = tG;
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            g_k_i.col(i) = (tG + G_k_i.block(0, i * n, n, n)
                            + (dks(i, k - 1) * In)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.block(0, i * n, n, n) = tG;
            dks(i, k) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf) {
    const int n = A1.size();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd tG(n);
    ArrayXXd G_k_i = ArrayXXd::Zero(n, p + 1);
    ArrayXXd g_k_i = ArrayXXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = G_k_i.col(i) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        tG = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = (tG + G_k_i.col(0) + dks(0, k - 1)) * mu + A2 * g_k_i.col(0);
        G_k_i.col(0) = tG;
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                 A2 * (dks(i, k - 1) + G_k_i.col(i));
            g_k_i.col(i) = (tG + G_k_i.col(i) + dks(i, k - 1)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.col(i) = tG;
            dks(i, k) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf.rightCols(m + 1 - k) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil2_pq_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int p, const int q) { //, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, q + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd G_k_i = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd g_k_i = MatrixXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        g_k_i.col(i) = G_k_i.block(0, i * n, n, n) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= q; k++) {
        G_k_i.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = G_k_i.block(0, 0, n, n) * mu + A2 * g_k_i.col(0);
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                                          A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            g_k_i.col(i) = G_k_i.block(0, i * n, n, n) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            dks(i, k) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * (k + i));
        }
        // if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
        //     dks /= 1e10;
        //     G_k_i /= 1e10;
        //     g_k_i /= 1e10;
        //     lscf -= log(1e10);
        // }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil2_pq_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int p, const int q) { //, Eigen::ArrayXXd& lscf) {
    const int n = A1.size();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, q + 1);
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXXd G_k_i = ArrayXXd::Zero(n, p + 1);
    ArrayXXd g_k_i = ArrayXXd::Zero(n, p + 1);
    for(int i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = G_k_i.col(i) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(int k = 1; k <= q; k++) {
        G_k_i.col(0) = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = G_k_i.col(0) * mu + A2 * g_k_i.col(0);
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            G_k_i.col(i) = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                           A2 * (dks(i, k - 1) + G_k_i.col(i));
            g_k_i.col(i) = G_k_i.col(i) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            dks(i, k) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * (k + i));
        }
        // if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
        //     dks /= 1e10;
        //     G_k_i /= 1e10;
        //     g_k_i /= 1e10;
        //     lscf -= log(1e10);
        // }
    }
    return dks;
}


// // [[Rcpp::export]]
Eigen::ArrayXXd d3_ijk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                          const int m, Eigen::ArrayXXd& lscf, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd Go = MatrixXd::Zero(n, n * (m + 1) * m);
    MatrixXd Gn = MatrixXd::Zero(n, n * (m + 1) * (m + 1));
    for(int i = 1; i <= m; i++) {
        Gn.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        dks(i, 0) = Gn.block(0, i * n, n, n).trace() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * n * (m + 1)) = Gn.block(0, 0, n, k * n * (m + 1));
        Gn.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        dks(0, k) = Gn.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gn.block(0, i * n, n, n) =
                A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            dks(i, k) = Gn.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(int j = 1; j < k; j++) {
            Gn.block(0, j * n * (m + 1), n, n) =
                A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (m + 1), n, n)) +
                A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (m + 1), n, n));
            dks(0, (k - j) + j * (m + 1)) = Gn.block(0, j * n * (m + 1), n, n).trace() / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                Gn.block(0, j * n * (m + 1) + i * n, n, n) =
                    A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (m + 1) + (i - 1) * n, n, n)) +
                    A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (m + 1) + i * n, n, n)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (m + 1) + i * n, n, n));
                dks(i, (k - j) + j * (m + 1)) = Gn.block(0, j * n * (m + 1) + i * n, n, n).trace() / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        Gn.block(0, k * n * (m + 1), n, n) = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (m + 1), n, n));
        dks(0, k * (m + 1)) = Gn.block(0, k * n * (m + 1), n, n).trace() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gn.block(0, k * n * (m + 1) + i * n, n, n) =
                A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (m + 1) + (i - 1) * n, n, n)) +
                A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (m + 1) + i * n, n, n));
            dks(i, k * (m + 1)) = Gn.block(0, k * n * (m + 1) + i * n, n, n).trace() / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), m + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d3_ijk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                          const int m, Eigen::ArrayXXd& lscf) { // , int nthreads) {
// #ifdef _OPENMP
//     if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
//     omp_set_num_threads(nthreads);
// #endif
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXXd Go = ArrayXXd::Zero(n, (m + 1) * m);
    ArrayXXd Gn = ArrayXXd::Zero(n, (m + 1) * (m + 1));
    for(int i = 1; i <= m; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        dks(i, 0) = Gn.col(i).sum() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * (m + 1)) = Gn.block(0, 0, n, k * (m + 1));
        Gn.col(0) = A2 * (dks(0, k - 1) + Go.col(0));
        dks(0, k) = Gn.col(0).sum() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gn.col(i) = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                        A2 * (dks(i, k - 1) + Go.col(i));
            dks(i, k) = Gn.col(i).sum() / (2 * (k + i));
        }
// #ifdef _OPENMP
// #pragma omp parallel
// {
// #pragma omp for
// #endif
        for(int j = 1; j < k; j++) {
            Gn.col(j * (m + 1)) =
                A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (m + 1))) +
                A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (m + 1)));
            dks(0, (k - j) + j * (m + 1)) = Gn.col(j * (m + 1)).sum() / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                Gn.col(j * (m + 1) + i) =
                    A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gn.col(j * (m + 1) + i - 1)) +
                    A2 * (dks(i, (k - j - 1) + j * (m + 1)) + Go.col(j * (m + 1) + i)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (m + 1) + i));
                dks(i, (k - j) + j * (m + 1)) = Gn.col(j * (m + 1) + i).sum() / (2 * (k + i));
            }
        }
// #ifdef _OPENMP
// }
// #endif
        Gn.col(k * (m + 1)) = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (m + 1)));
        dks(0, k * (m + 1)) = Gn.col(k * (m + 1)).sum() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gn.col(k * (m + 1) + i) =
                A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (m + 1) + i - 1)) +
                A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (m + 1) + i));
            dks(i, k * (m + 1)) = Gn.col(k * (m + 1) + i).sum() / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), m + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                             const int m, const int p, Eigen::ArrayXXd& lscf, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd Go = MatrixXd::Zero(n, n * (p + 1) * m);
    MatrixXd Gn = MatrixXd::Zero(n, n * (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        dks(i, 0) = Gn.block(0, i * n, n, n).trace() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        Gn.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        dks(0, k) = Gn.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gn.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                                       A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            dks(i, k) = Gn.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(int j = 1; j < k; j++) {
            Gn.block(0, j * n * (p + 1), n, n) = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
            dks(0, (k - j) + j * (m + 1)) = Gn.block(0, j * n * (p + 1), n, n).trace() / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gn.block(0, j * n * (p + 1) + i * n, n, n) = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                                           A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                                           A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                dks(i, (k - j) + j * (m + 1)) = Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        Gn.block(0, k * n * (p + 1), n, n) = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
        dks(0, k * (m + 1)) = Gn.block(0, k * n * (p + 1), n, n).trace() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gn.block(0, k * n * (p + 1) + i * n, n, n) = A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                                       A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
            dks(i, k * (m + 1)) = Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const int m, const int p, Eigen::ArrayXXd& lscf) { // , int nthreads) {
// #ifdef _OPENMP
//     if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
//     omp_set_num_threads(nthreads);
// #endif
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXXd Go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd Gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        dks(i, 0) = Gn.col(i).sum() / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        Gn.col(0) = A2 * (dks(0, k - 1) + Go.col(0));
        dks(0, k) = Gn.col(0).sum() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gn.col(i) = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                        A2 * (dks(i, k - 1) + Go.col(i));
            dks(i, k) = Gn.col(i).sum() / (2 * (k + i));
        }
// #ifdef _OPENMP
// #pragma omp parallel
// {
// #pragma omp for
// #endif
        for(int j = 1; j < k; j++) {
            Gn.col(j * (p + 1)) =
                A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1))) +
                A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1)));
            dks(0, (k - j) + j * (m + 1)) = Gn.col(j * (p + 1)).sum() / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gn.col(j * (p + 1) + i) =
                    A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gn.col(j * (p + 1) + (i - 1))) +
                    A2 * (dks(i, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1) + i)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1) + i));
                dks(i, (k - j) + j * (m + 1)) = Gn.col(j * (p + 1) + i).sum() / (2 * (k + i));
            }
        }
// #ifdef _OPENMP
// }
// #endif
        Gn.col(k * (p + 1)) = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1)));
        dks(0, k * (m + 1)) = Gn.col(k * (p + 1)).sum() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gn.col(k * (p + 1) + i) =
                A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (p + 1) + i - 1)) +
                A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1) + i));
            dks(i, k * (m + 1)) = Gn.col(k * (p + 1) + i).sum() / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h3_ijk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd& mu, const int m, Eigen::ArrayXXd& lscf, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd Go = MatrixXd::Zero(n, n * (m + 1) * m);
    MatrixXd Gn = MatrixXd::Zero(n, n * (m + 1) * (m + 1));
    MatrixXd go = MatrixXd::Zero(n, (m + 1) * m);
    MatrixXd gn = MatrixXd::Zero(n, (m + 1) * (m + 1));
    for(int i = 1; i <= m; i++) {
        tG = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        gn.col(i) = (tG - Gn.block(0, (i - 1) * n, n, n) - (dks(i - 1, 0) * In)) * mu + A1 * gn.col(i - 1);
        Gn.block(0, i * n, n, n) = tG;
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * n * (m + 1)) = Gn.block(0, 0, n, k * n * (m + 1));
        go.block(0, 0, n, k * (m + 1)) = gn.block(0, 0, n, k * (m + 1));
        tG = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        gn.col(0) = (tG - Go.block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            gn.col(i) = (tG - Gn.block(0, (i - 1) * n, n, n) - Go.block(0, i * n, n, n)
                         - ((dks(i - 1, k) + dks(i, (k - 1))) * In)) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (m + 1), n, n)) +
                 A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (m + 1), n, n));
            gn.col(j * (m + 1)) =
                (tG - Go.block(0, j * n * (m + 1), n, n) - Go.block(0, (j - 1) * n * (m + 1), n, n)
                 - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                A2 * go.col(j * (m + 1)) + A3 * go.col((j - 1) * (m + 1));
            Gn.block(0, j * n * (m + 1), n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (m + 1), n, n).trace() + gn.col(j * (m + 1)).dot(mu)) / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (m + 1) + (i - 1) * n, n, n)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (m + 1) + i * n, n, n)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (m + 1) + i * n, n, n));
                gn.col(j * (m + 1) + i) =
                    (tG - Gn.block(0, j * n * (m + 1) + (i - 1) * n, n, n) - Go.block(0, j * n * (m + 1) + i * n, n, n) - Go.block(0, (j - 1) * n * (m + 1) + i * n, n, n)
                     - ((dks(i - 1, (k - j) + j * (m + 1)) + dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                    A1 * gn.col(j * (m + 1) + i - 1) + A2 * go.col(j * (m + 1) + i) + A3 * go.col((j - 1) * (m + 1) + i);
                Gn.block(0, j * n * (m + 1) + i * n, n, n) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (m + 1) + i * n, n, n).trace() + gn.col(j * (m + 1) + i).dot(mu)) / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (m + 1), n, n));
        gn.col(k * (m + 1)) = (tG - Go.block(0, (k - 1) * n * (m + 1), n, n)
                               - (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * go.col((k - 1) * (m + 1));
        Gn.block(0, k * n * (m + 1), n, n) = tG;
        dks(0, k * (m + 1)) = (Gn.block(0, k * n * (m + 1), n, n).trace() + gn.col(k * (m + 1)).dot(mu)) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (m + 1) + (i - 1) * n, n, n)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (m + 1) + i * n, n, n));
            gn.col(k * (m + 1) + i) =
                (tG - Gn.block(0, k * n * (m + 1) + (i - 1) * n, n, n) - Go.block(0, (k - 1) * n * (m + 1) + i * n, n, n)
                 - ((dks(i - 1, k * (m + 1)) + dks(i, (k - 1) * (m + 1))) * In)) * mu +
                A1 * gn.col(k * (m + 1) + i - 1) + A3 * go.col((k - 1) * (m + 1) + i);
            Gn.block(0, k * n * (m + 1) + i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gn.block(0, k * n * (m + 1) + i * n, n, n).trace() + gn.col(k * (m + 1) + i).dot(mu)) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), m + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h3_ijk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int m, Eigen::ArrayXXd& lscf) { // , int nthreads) {
// #ifdef _OPENMP
//     if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
//     omp_set_num_threads(nthreads);
// #endif
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXXd tG(n, n);
    ArrayXXd Go = ArrayXXd::Zero(n, (m + 1) * m);
    ArrayXXd Gn = ArrayXXd::Zero(n, (m + 1) * (m + 1));
    ArrayXXd go = ArrayXXd::Zero(n, (m + 1) * m);
    ArrayXXd gn = ArrayXXd::Zero(n, (m + 1) * (m + 1));
    for(int i = 1; i <= m; i++) {
        tG = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = (tG - Gn.col(i - 1) - dks(i - 1, 0)) * mu + A1 * gn.col(i - 1);
        Gn.col(i) = tG;
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * (m + 1)) = Gn.block(0, 0, n, k * (m + 1));
        go.block(0, 0, n, k * (m + 1)) = gn.block(0, 0, n, k * (m + 1));
        tG = A2 * (dks(0, k - 1) + Go.col(0));
        gn.col(0) = (tG - Go.col(0) - (dks(0, k - 1))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                 A2 * (dks(i, k - 1) + Go.col(i));
            gn.col(i) = (tG - Gn.col(i - 1) - Go.col(i)
                         - (dks(i - 1, k) + dks(i, k - 1))) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.col(i) = tG;
            dks(i, k) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
        }
// #ifdef _OPENMP
// #pragma omp parallel private(tG)
// {
// #pragma omp for
// #endif
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (m + 1))) +
                 A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (m + 1)));
            gn.col(j * (m + 1)) =
                (tG - Go.col(j * (m + 1)) - Go.col((j - 1) * (m + 1))
                 - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu +
                A2 * go.col(j * (m + 1)) + A3 * go.col((j - 1) * (m + 1));
            Gn.col(j * (m + 1)) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.col(j * (m + 1)).sum() + (mu * gn.col(j * (m + 1))).sum()) / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gn.col(j * (m + 1) + i - 1)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) + Go.col(j * (m + 1) + i)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (m + 1) + i));
                gn.col(j * (m + 1) + i) =
                    (tG - Gn.col(j * (m + 1) + i - 1) - Go.col(j * (m + 1) + i) - Go.col((j - 1) * (m + 1) + i)
                     - ((dks(i - 1, (k - j) + j * (m + 1)) + dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))))) * mu +
                    A1 * gn.col(j * (m + 1) + i - 1) + A2 * go.col(j * (m + 1) + i) + A3 * go.col((j - 1) * (m + 1) + i);
                Gn.col(j * (m + 1) + i) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gn.col(j * (m + 1) + i).sum() + (mu * gn.col(j * (m + 1) + i)).sum()) / (2 * (k + i));
            }
        }
// #ifdef _OPENMP
// }
// #endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (m + 1)));
        gn.col(k * (m + 1)) = (tG - Go.col((k - 1) * (m + 1)) - (dks(0, (k - 1) * (m + 1)))) * mu + A3 * go.col((k - 1) * (m + 1));
        Gn.col(k * (m + 1)) = tG;
        dks(0, k * (m + 1)) = (Gn.col(k * (m + 1)).sum() + (mu * gn.col(k * (m + 1))).sum()) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (m + 1) + i - 1)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (m + 1) + i));
            gn.col(k * (m + 1) + i) =
                (tG - Gn.col(k * (m + 1) + i - 1) - Go.col((k - 1) * (m + 1) + i)
                 - (dks(i - 1, k * (m + 1)) + dks(i, (k - 1) * (m + 1)))) * mu +
                A1 * gn.col(k * (m + 1) + i - 1) + A3 * go.col((k - 1) * (m + 1) + i);
            Gn.col(k * (m + 1) + i) = tG;
            dks(i, k * (m + 1)) = (Gn.col(k * (m + 1) + i).sum() + (mu * gn.col(k * (m + 1) + i)).sum()) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), m + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd Go = MatrixXd::Zero(n, n * (p + 1) * m);
    MatrixXd Gn = MatrixXd::Zero(n, n * (p + 1) * (m + 1));
    MatrixXd go = MatrixXd::Zero(n, (p + 1) * m);
    MatrixXd gn = MatrixXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        gn.col(i) = Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        gn.col(0) = (tG - Go.block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            gn.col(i) = (tG - Go.block(0, i * n, n, n)
                         - (dks(i, (k - 1)) * In)) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) +
                 A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
            gn.col(j * (p + 1)) =
                (tG - Go.block(0, j * n * (p + 1), n, n) - Go.block(0, (j - 1) * n * (p + 1), n, n) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.block(0, j * n * (p + 1), n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                gn.col(j * (p + 1) + i) =
                    (tG - Go.block(0, j * n * (p + 1) + i * n, n, n) - Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n)
                     - ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                    A1 * gn.col(j * (p + 1) + i - 1) + A2 * go.col(j * (p + 1) + i) + A3 * go.col((j - 1) * (p + 1) + i);
                Gn.block(0, j * n * (p + 1) + i * n, n, n) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() + gn.col(j * (p + 1) + i).dot(mu)) / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
        gn.col(k * (p + 1)) = (tG - Go.block(0, (k - 1) * n * (p + 1), n, n) - (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * go.col((k - 1) * (p + 1));
        Gn.block(0, k * n * (p + 1), n, n) = tG;
        dks(0, k * (m + 1)) = (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
            gn.col(k * (p + 1) + i) = (tG - Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n)
                             - (dks(i, (k - 1) * (m + 1)) * In)) * mu +
                            A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.block(0, k * n * (p + 1) + i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() + gn.col(k * (p + 1) + i).dot(mu)) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf) { // , int nthreads) {
// #ifdef _OPENMP
//     if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
//     omp_set_num_threads(nthreads);
// #endif
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd tG(n);
    ArrayXXd Go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd Gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    ArrayXXd go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks(0, k - 1) + Go.col(0));
        gn.col(0) = (tG - Go.col(0) - (dks(0, k - 1))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                 A2 * (dks(i, k - 1) + Go.col(i));
            gn.col(i) = (tG - Go.col(i)
                         - (dks(i, (k - 1)))) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.col(i) = tG;
            dks(i, k) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
        }
// #ifdef _OPENMP
// #pragma omp parallel private(tG)
// {
// #pragma omp for
// #endif
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1))) +
                 A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1)));
            gn.col(j * (p + 1)) =
                (tG - Go.col(j * (p + 1)) - Go.col((j - 1) * (p + 1)) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.col(j * (p + 1)) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gn.col(j * (p + 1) + i - 1)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1) + i)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1) + i));
                gn.col(j * (p + 1) + i) = (tG - Go.col(j * (p + 1) + i) - Go.col((j - 1) * (p + 1) + i)
                             - ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))))) * mu +
                            A1 * gn.col(j * (p + 1) + i - 1) + A2 * go.col(j * (p + 1) + i) + A3 * go.col((j - 1) * (p + 1) + i);
                Gn.col(j * (p + 1) + i) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gn.col(j * (p + 1) + i).sum() + (mu * gn.col(j * (p + 1) + i)).sum()) / (2 * (k + i));
            }
        }
// #ifdef _OPENMP
// }
// #endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1)));
        gn.col(k * (p + 1)) =
            (tG - Go.col((k - 1) * (p + 1)) - (dks(0, (k - 1) * (m + 1)))) * mu +
            A3 * go.col((k - 1) * (p + 1));
        Gn.col(k * (p + 1)) = tG;
        dks(0, k * (m + 1)) = (Gn.col(k * (p + 1)).sum() + (mu * gn.col(k * (p + 1))).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (p + 1) + i - 1)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1) + i));
            gn.col(k * (p + 1) + i) = (tG - Go.col((k - 1) * (p + 1) + i)
                         - (dks(i, (k - 1) * (m + 1)))) * mu +
                        A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.col(k * (p + 1) + i) = tG;
            dks(i, k * (m + 1)) = (Gn.col(k * (p + 1) + i).sum() + (mu * gn.col(k * (p + 1) + i)).sum()) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    MatrixXd tG(n, n);
    MatrixXd Go = MatrixXd::Zero(n, n * (p + 1) * m);
    MatrixXd Gn = MatrixXd::Zero(n, n * (p + 1) * (m + 1));
    MatrixXd go = MatrixXd::Zero(n, (p + 1) * m);
    MatrixXd gn = MatrixXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        gn.col(i) = Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        gn.col(0) = (tG + Go.block(0, 0, n, n) + (dks(0, k - 1) * In)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            gn.col(i) = (tG + Go.block(0, i * n, n, n)
                         + (dks(i, (k - 1)) * In)) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) +
                 A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
            gn.col(j * (p + 1)) =
                (tG + Go.block(0, j * n * (p + 1), n, n) + Go.block(0, (j - 1) * n * (p + 1), n, n) + ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.block(0, j * n * (p + 1), n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                gn.col(j * (p + 1) + i) =
                    (tG + Go.block(0, j * n * (p + 1) + i * n, n, n) + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n)
                     + ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                    A1 * gn.col(j * (p + 1) + i - 1) + A2 * go.col(j * (p + 1) + i) + A3 * go.col((j - 1) * (p + 1) + i);
                Gn.block(0, j * n * (p + 1) + i * n, n, n) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() + gn.col(j * (p + 1) + i).dot(mu)) / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
        gn.col(k * (p + 1)) = (tG + Go.block(0, (k - 1) * n * (p + 1), n, n) + (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * go.col((k - 1) * (p + 1));
        Gn.block(0, k * n * (p + 1), n, n) = tG;
        dks(0, k * (m + 1)) = (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
            gn.col(k * (p + 1) + i) = (tG + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n)
                             + (dks(i, (k - 1) * (m + 1)) * In)) * mu +
                            A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.block(0, k * n * (p + 1) + i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() + gn.col(k * (p + 1) + i).dot(mu)) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf) { // , int nthreads) {
// #ifdef _OPENMP
//     if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
//     omp_set_num_threads(nthreads);
// #endif
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd tG(n);
    ArrayXXd Go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd Gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    ArrayXXd go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks(0, k - 1) + Go.col(0));
        gn.col(0) = (tG + Go.col(0) + (dks(0, k - 1))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                 A2 * (dks(i, k - 1) + Go.col(i));
            gn.col(i) = (tG + Go.col(i)
                         + (dks(i, (k - 1)))) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.col(i) = tG;
            dks(i, k) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
        }
// #ifdef _OPENMP
// #pragma omp parallel private(tG)
// {
// #pragma omp for
// #endif
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1))) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1)));
            gn.col(j * (p + 1)) =
                (tG + Go.col(j * (p + 1)) + Go.col((j - 1) * (p + 1)) +
                 ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.col(j * (p + 1)) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gn.col(j * (p + 1) + i - 1)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1) + i)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1) + i));
                gn.col(j * (p + 1) + i) =
                    (tG + Go.col(j * (p + 1) + i) + Go.col((j - 1) * (p + 1) + i)
                     + ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))))) * mu +
                    A1 * gn.col(j * (p + 1) + i - 1) + A2 * go.col(j * (p + 1) + i) + A3 * go.col((j - 1) * (p + 1) + i);
                Gn.col(j * (p + 1) + i) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gn.col(j * (p + 1) + i).sum() + (mu * gn.col(j * (p + 1) + i)).sum()) / (2 * (k + i));
            }
        }
// #ifdef _OPENMP
// }
// #endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1)));
        gn.col(k * (p + 1)) =
            (tG + Go.col((k - 1) * (p + 1)) + (dks(0, (k - 1) * (m + 1)))) * mu +
            A3 * go.col((k - 1) * (p + 1));
        Gn.col(k * (p + 1)) = tG;
        dks(0, k * (m + 1)) = (Gn.col(k * (p + 1)).sum() + (mu * gn.col(k * (p + 1))).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (p + 1) + i - 1)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1) + i));
            gn.col(k * (p + 1) + i) =
                (tG + Go.col((k - 1) * (p + 1) + i)
                 + (dks(i, (k - 1) * (m + 1)))) * mu +
                A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.col(k * (p + 1) + i) = tG;
            dks(i, k * (m + 1)) = (Gn.col(k * (p + 1) + i).sum() + (mu * gn.col(k * (p + 1) + i)).sum()) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(int j = 0; j <=k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            for(int j = 0; j <=k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil3_pqr_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd mu, const int p, const int q, const int r) { //, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    const int m = q + r;
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (q + 1) * (r + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_np = MatrixXd::Zero(n, n * (p + 1));
    const MatrixXd zeromat_n_p = MatrixXd::Zero(n, p + 1);
    MatrixXd Go = MatrixXd::Zero(n, n * (p + 1) * m);
    MatrixXd Gn = MatrixXd::Zero(n, n * (p + 1) * (m + 1));
    MatrixXd go = MatrixXd::Zero(n, (p + 1) * m);
    MatrixXd gn = MatrixXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        gn.col(i) = Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        if(k <= q) {
            Gn.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
            gn.col(0) = Gn.block(0, 0, n, n) * mu + A2 * go.col(0);
            dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gn.block(0, i * n, n, n) =
                    A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                    A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
                gn.col(i) = Gn.block(0, i * n, n, n) * mu +
                            A1 * gn.col(i - 1) + A2 * go.col(i);
                dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
            }
        } else {
            Gn.block(0, 0, n, n * (p + 1)) = zeromat_n_p;
            gn.block(0, 0, n, p + 1) = zeromat_n_p;
        }
        for(int j = 1; j < k; j++) {
            if(((k - j) <= q) && (j <= r)) {
                Gn.block(0, j * n * (p + 1), n, n) =
                    A2 * (dks(0, (k - j - 1) + j * (q + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) +
                    A3 * (dks(0, (k - j) + (j - 1) * (q + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
                gn.col(j * (p + 1)) =
                    Gn.block(0, j * n * (p + 1), n, n) * mu +
                    A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
                dks(0, (k - j) + j * (q + 1)) = (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
                for(int i = 1; i <= p; i++) {
                    Gn.block(0, j * n * (p + 1) + i * n, n, n) =
                        A1 * (dks(i - 1, (k - j) + j * (q + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                        A2 * (dks(i, (k - j - 1) + j * (q + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                        A3 * (dks(i, (k - j) + (j - 1) * (q + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                    gn.col(j * (p + 1) + i) =
                        Gn.block(0, j * n * (p + 1) + i * n, n, n) * mu +
                        A1 * gn.col(j * (p + 1) + i - 1) + A2 * go.col(j * (p + 1) + i) + A3 * go.col((j - 1) * (p + 1) + i);
                    dks(i, (k - j) + j * (q + 1)) = (Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() + gn.col(j * (p + 1) + i).dot(mu)) / (2 * (k + i));
                }
            } else {
                Gn.block(0, j * n * (p + 1), n, n * (p + 1)) = zeromat_n_p;
                gn.block(0, j * (p + 1), n, p + 1) = zeromat_n_p;
            }
        }
        if(k <= r) {
            Gn.block(0, k * n * (p + 1), n, n) = A3 * (dks(0, (k - 1) * (q + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
            gn.col(k * (p + 1)) = Gn.block(0, k * n * (p + 1), n, n) * mu + A3 * go.col((k - 1) * (p + 1));
            dks(0, k * (q + 1)) = (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gn.block(0, k * n * (p + 1) + i * n, n, n) =
                    A1 * (dks(i - 1, k * (q + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                    A3 * (dks(i, (k - 1) * (q + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
                gn.col(k * (p + 1) + i) =
                    Gn.block(0, k * n * (p + 1) + i * n, n, n) * mu +
                    A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
                dks(i, k * (q + 1)) = (Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() + gn.col(k * (p + 1) + i).dot(mu)) / (2 * (k + i));
            }
        } else {
            Gn.block(0, k * n * (p + 1), n, n * (p + 1)) = zeromat_n_p;
            gn.block(0, k * (p + 1), n, p + 1) = zeromat_n_p;
        }
        // if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
        //     dks /= 1e10;
        //     Gn /= 1e10;
        //     gn /= 1e10;
        //     lscf -= log(1e10);
        // }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil3_pqr_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int p, const int q, const int r) { //, Eigen::ArrayXXd& lscf) {
    const int n = A1.rows();
    const int m = q + r;
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (q + 1) * (r + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_p = ArrayXXd::Zero(n, p + 1);
    ArrayXXd Go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd Gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    ArrayXXd go = ArrayXXd::Zero(n, (p + 1) * m);
    ArrayXXd gn = ArrayXXd::Zero(n, (p + 1) * (m + 1));
    for(int i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(int k = 1; k <= m; k++) {
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        if(k <= q) {
            Gn.col(0) = A2 * (dks(0, k - 1) + Go.col(0));
            gn.col(0) = Gn.col(0) * mu + A2 * go.col(0);
            dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gn.col(i) = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                            A2 * (dks(i, k - 1) + Go.col(i));
                gn.col(i) = Gn.col(i) * mu +
                            A1 * gn.col(i - 1) + A2 * go.col(i);
                dks(i, k) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
            }
        } else {
            Gn.block(0, 0, n, p + 1) = zeromat_n_p;
            gn.block(0, 0, n, p + 1) = zeromat_n_p;
        }
        for(int j = 1; j < k; j++) {
            if(((k - j) <= q) && (j <= r)) {
                Gn.col(j * (p + 1)) =
                    A2 * (dks(0, (k - j - 1) + j * (q + 1)) + Go.col(j * (p + 1))) +
                    A3 * (dks(0, (k - j) + (j - 1) * (q + 1)) + Go.col((j - 1) * (p + 1)));
                gn.col(j * (p + 1)) =
                    Gn.col(j * (p + 1)) * mu +
                    A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
                dks(0, (k - j) + j * (q + 1)) = (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
                for(int i = 1; i <= p; i++) {
                    Gn.col(j * (p + 1) + i) =
                        A1 * (dks(i - 1, (k - j) + j * (q + 1)) + Gn.col(j * (p + 1) + i - 1)) +
                        A2 * (dks(i, (k - j - 1) + j * (q + 1)) + Go.col(j * (p + 1) + i)) +
                        A3 * (dks(i, (k - j) + (j - 1) * (q + 1)) + Go.col((j - 1) * (p + 1) + i));
                    gn.col(j * (p + 1) + i) =
                        Gn.col(j * (p + 1) + i) * mu +
                        A1 * gn.col(j * (p + 1) + i - 1) +
                        A2 * go.col(j * (p + 1) + i) + A3 *
                        go.col((j - 1) * (p + 1) + i);
                    dks(i, (k - j) + j * (q + 1)) = (Gn.col(j * (p + 1) + i).sum() + (mu * gn.col(j * (p + 1) + i)).sum()) / (2 * (k + i));
                }
            } else {
                Gn.block(0, j * (p + 1), n, p + 1) = zeromat_n_p;
                gn.block(0, j * (p + 1), n, p + 1) = zeromat_n_p;
            }
        }
        if(k <= r) {
            Gn.col(k * (p + 1)) = A3 * (dks(0, (k - 1) * (q + 1)) + Go.col((k - 1) * (p + 1)));
            gn.col(k * (p + 1)) =
                Gn.col(k * (p + 1)) * mu +
                A3 * go.col((k - 1) * (p + 1));
            dks(0, k * (q + 1)) = (Gn.col(k * (p + 1)).sum() + (mu * gn.col(k * (p + 1))).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gn.col(k * (p + 1) + i) =
                    A1 * (dks(i - 1, k * (q + 1)) + Gn.col(k * (p + 1) + i - 1)) +
                    A3 * (dks(i, (k - 1) * (q + 1)) + Go.col((k - 1) * (p + 1) + i));
                gn.col(k * (p + 1) + i) = Gn.col(k * (p + 1) + i) * mu +
                            A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
                dks(i, k * (q + 1)) = (Gn.col(k * (p + 1) + i).sum() + (mu * gn.col(k * (p + 1) + i)).sum()) / (2 * (k + i));
            }
        } else {
            Gn.block(0, k * (p + 1), n, p + 1) = zeromat_n_p;
            gn.block(0, k * (p + 1), n, p + 1) = zeromat_n_p;
        }
        // if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
        //     dks /= 1e10;
        //     Gn /= 1e10;
        //     gn /= 1e10;
        //     lscf -= log(1e10);
        // }
    }
    return dks;
}
