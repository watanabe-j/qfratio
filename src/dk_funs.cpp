#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <limits>
#include "dk_funs.h"

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// // [[Rcpp::export]]
Eigen::ArrayXd d1_i_vE(const Eigen::ArrayXd& L, const int m, double& lscf) {
    int n = L.size();
    ArrayXd dks = ArrayXd::Zero(m + 1);
    dks(0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    ArrayXd uk = ArrayXd::Zero(n);
    for(int k = 1; k <= m; k++) {
        uk = L * (dks(k - 1) + uk);
        dks(k) = uk.sum() / (2 * k);
        if(uk.maxCoeff() > thr) {
            dks /= 1e10;
            uk /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXd d1_i_mE(const Eigen::MatrixXd& A, const int m, double& lscf) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A, Eigen::EigenvaluesOnly);
    ArrayXd L = eigA.eigenvalues();
    return d1_i_vE(L, m, lscf);
}

// // [[Rcpp::export]]
Eigen::ArrayXd dtil1_i_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& mu,
    const int m, double& lscf) {
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
            dks /= 1e10;
            uk /= 1e10;
            vk /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXd dtil1_i_mE(const Eigen::MatrixXd& A, const Eigen::VectorXd& mu,
        const int m, double& lscf) {
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
    ArrayXd L = eigA.eigenvalues();
    ArrayXd mud = eigA.eigenvectors().transpose() * mu;
    return dtil1_i_vE(L, mud, m, lscf);
}

// // [[Rcpp::export]]
Eigen::ArrayXXd arl_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& D, const int m) {
    const int n = L.size();
    double lscf = 0;
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
         const int m, const int p, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const int m, const int p, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}


// // [[Rcpp::export]]
Eigen::ArrayXXd d2_ij_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const int m, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d2_ij_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const int m, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h2_ij_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h2_ij_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int m, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, const int p, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int m, const int p, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, const int p, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
         const Eigen::ArrayXd& mu, const int m, const int p, double& lscf) {
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
            dks /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil2_pq_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int p, const int q) { //, double& lscf) {
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
         const Eigen::ArrayXd& mu, const int p, const int q) { //, double& lscf) {
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
                            const int m, double& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_nm = MatrixXd::Zero(n, n * (m + 1));
    MatrixXd Gc = MatrixXd::Zero(n, n * (m + 1));
    MatrixXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_nm;
    MatrixXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_nm;
    for(int i = 1; i <= m; i++) {
        Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gc.block(0, (i - 1) * n, n, n));
        dks(i, 0) = Gc.block(0, i * n, n, n).trace() / (2 * i);
    }
    Gn[0] = Gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        Gc.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + (Go[0]).block(0, 0, n, n));
        dks(0, k) = Gc.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + Gc.block(0, (i - 1) * n, n, n)) +
                                       A2 * (dks(i, k - 1) * In + (Go[0]).block(0, i * n, n, n));
            dks(i, k) = Gc.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        Gn[0] = Gc;
        for(int j = 1; j < k; j++) {
            Gc.block(0, 0, n, n) = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, 0, n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, 0, n, n));
            dks(0, (k - j) + j * (m + 1)) = Gc.block(0, 0, n, n).trace() / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                                           A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, i * n, n, n)) +
                                           A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, i * n, n, n));
                dks(i, (k - j) + j * (m + 1)) = Gc.block(0, i * n, n, n).trace() / (2 * (k + i));
            }
            Gn[j] = Gc;
        }
        Gc.block(0, 0, n, n) = A3 * (dks(0, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, 0, n, n));
        dks(0, k * (m + 1)) = Gc.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, k * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                                       A3 * (dks(i, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, i * n, n, n));
            dks(i, k * (m + 1)) = Gc.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        Gn[k] = Gc;
        if(Gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d3_ijk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const int m, double& lscf) {
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_m = ArrayXXd::Zero(n, m + 1);
    ArrayXXd Gc = ArrayXXd::Zero(n, (m + 1));
    ArrayXXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_m;
    ArrayXXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_m;
    for(int i = 1; i <= m; i++) {
        Gc.col(i) = A1 * (dks(i - 1, 0) + Gc.col(i - 1));
        dks(i, 0) = Gc.col(i).sum() / (2 * i);
    }
    Gn[0] = Gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        Gc.col(0) = A2 * (dks(0, k - 1) + (Go[0]).col(0));
        dks(0, k) = Gc.col(0).sum() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gc.col(i) = A1 * (dks(i - 1, k) + Gc.col((i - 1))) +
                        A2 * (dks(i, k - 1) + (Go[0]).col(i));
            dks(i, k) = Gc.col(i).sum() / (2 * (k + i));
        }
        Gn[0] = Gc;
        for(int j = 1; j < k; j++) {
            Gc.col(0) = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + (Go[j]).col(0)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(0));
            dks(0, (k - j) + j * (m + 1)) = Gc.col(0).sum() / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                Gc.col(i) = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gc.col((i - 1))) +
                            A2 * (dks(i, (k - j - 1) + j * (m + 1)) + (Go[j]).col(i)) +
                            A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(i));
                dks(i, (k - j) + j * (m + 1)) = Gc.col(i).sum() / (2 * (k + i));
            }
            Gn[j] = Gc;
        }
        Gc.col(0) = A3 * (dks(0, (k - 1) * (m + 1)) + (Go[k - 1]).col(0));
        dks(0, k * (m + 1)) = Gc.col(0).sum() / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            Gc.col(i) = A1 * (dks(i - 1, k * (m + 1)) + Gc.col((i - 1))) +
                        A3 * (dks(i, (k - 1) * (m + 1)) + (Go[k - 1]).col(i));
            dks(i, k * (m + 1)) = Gc.col(i).sum() / (2 * (k + i));
        }
        Gn[k] = Gc;
        if(Gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                             const int m, const int p, double& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_np = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd Gc = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_np;
    MatrixXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_np;
    for(int i = 1; i <= p; i++) {
        Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gc.block(0, (i - 1) * n, n, n));
        dks(i, 0) = Gc.block(0, i * n, n, n).trace() / (2 * i);
    }
    Gn[0] = Gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        Gc.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + (Go[0]).block(0, 0, n, n));
        dks(0, k) = Gc.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + Gc.block(0, (i - 1) * n, n, n)) +
                                       A2 * (dks(i, k - 1) * In + (Go[0]).block(0, i * n, n, n));
            dks(i, k) = Gc.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        Gn[0] = Gc;
        for(int j = 1; j < k; j++) {
            Gc.block(0, 0, n, n) = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, 0, n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, 0, n, n));
            dks(0, (k - j) + j * (m + 1)) = Gc.block(0, 0, n, n).trace() / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                                           A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, i * n, n, n)) +
                                           A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, i * n, n, n));
                dks(i, (k - j) + j * (m + 1)) = Gc.block(0, i * n, n, n).trace() / (2 * (k + i));
            }
            Gn[j] = Gc;
        }
        Gc.block(0, 0, n, n) = A3 * (dks(0, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, 0, n, n));
        dks(0, k * (m + 1)) = Gc.block(0, 0, n, n).trace() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, k * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                                       A3 * (dks(i, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, i * n, n, n));
            dks(i, k * (m + 1)) = Gc.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        Gn[k] = Gc;
        if(Gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd d3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const int m, const int p, double& lscf) {
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_p = ArrayXXd::Zero(n, p + 1);
    ArrayXXd Gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_p;
    ArrayXXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.col(i) = A1 * (dks(i - 1, 0) + Gc.col(i - 1));
        dks(i, 0) = Gc.col(i).sum() / (2 * i);
    }
    Gn[0] = Gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        Gc.col(0) = A2 * (dks(0, k - 1) + (Go[0]).col(0));
        dks(0, k) = Gc.col(0).sum() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gc.col(i) = A1 * (dks(i - 1, k) + Gc.col((i - 1))) +
                        A2 * (dks(i, k - 1) + (Go[0]).col(i));
            dks(i, k) = Gc.col(i).sum() / (2 * (k + i));
        }
        Gn[0] = Gc;
        for(int j = 1; j < k; j++) {
            Gc.col(0) = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + (Go[j]).col(0)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(0));
            dks(0, (k - j) + j * (m + 1)) = Gc.col(0).sum() / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gc.col(i) = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gc.col((i - 1))) +
                            A2 * (dks(i, (k - j - 1) + j * (m + 1)) + (Go[j]).col(i)) +
                            A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(i));
                dks(i, (k - j) + j * (m + 1)) = Gc.col(i).sum() / (2 * (k + i));
            }
            Gn[j] = Gc;
        }
        Gc.col(0) = A3 * (dks(0, (k - 1) * (m + 1)) + (Go[k - 1]).col(0));
        dks(0, k * (m + 1)) = Gc.col(0).sum() / (2 * k);
        for(int i = 1; i <= p; i++) {
            Gc.col(i) = A1 * (dks(i - 1, k * (m + 1)) + Gc.col((i - 1))) +
                        A3 * (dks(i, (k - 1) * (m + 1)) + (Go[k - 1]).col(i));
            dks(i, k * (m + 1)) = Gc.col(i).sum() / (2 * (k + i));
        }
        Gn[k] = Gc;
        if(Gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h3_ijk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd& mu, const int m, double& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_nm = MatrixXd::Zero(n, n * (m + 1));
    const MatrixXd zeromat_n_m = MatrixXd::Zero(n, m + 1);
    MatrixXd Gc = MatrixXd::Zero(n, n * (m + 1));
    MatrixXd gc = MatrixXd::Zero(n, (m + 1));
    MatrixXd tG(n, n);
    MatrixXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_nm;
    MatrixXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_m;
    MatrixXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_nm;
    MatrixXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_m;
    for(int i = 1; i <= m; i++) {
        tG = A1 * (dks(i - 1, 0) * In + Gc.block(0, (i - 1) * n, n, n));
        gc.col(i) = (tG - Gc.block(0, (i - 1) * n, n, n) - (dks(i - 1, 0) * In)) * mu + A1 * gc.col(i - 1);
        Gc.block(0, i * n, n, n) = tG;
        dks(i, 0) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        tG = A2 * (dks(0, k - 1) * In + (Go[0]).block(0, 0, n, n));
        gc.col(0) = (tG - (Go[0]).block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * (go[0]).col(0);
        Gc.block(0, 0, n, n) = tG;
        dks(0, k) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k) * In + Gc.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + (Go[0]).block(0, i * n, n, n));
            gc.col(i) = (tG - Gc.block(0, (i - 1) * n, n, n) - (Go[0]).block(0, i * n, n, n)
                         - ((dks(i - 1, k) + dks(i, (k - 1))) * In)) * mu +
                        A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
            Gc.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, 0, n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, 0, n, n));
            gc.col(0) = (tG - (Go[j]).block(0, 0, n, n) - (Go[j - 1]).block(0, 0, n, n) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
            Gc.block(0, 0, n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, i * n, n, n)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, i * n, n, n));
                gc.col(i) = (tG - Gc.block(0, (i - 1) * n, n, n) - (Go[j]).block(0, i * n, n, n) - (Go[j - 1]).block(0, i * n, n, n)
                             - ((dks(i - 1, (k - j) + j * (m + 1)) + dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                            A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                Gc.block(0, i * n, n, n) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        tG = A3 * (dks(0, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, 0, n, n));
        gc.col(0) = (tG - (Go[k - 1]).col(0) - (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * (go[k - 1]).col(0);
        Gc.block(0, 0, n, n) = tG;
        dks(0, k * (m + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, i * n, n, n));
            gc.col(i) = (tG - Gc.block(0, (i - 1) * n, n, n) - (Go[k - 1]).block(0, i * n, n, n)
                         - ((dks(i - 1, k * (m + 1)) + dks(i, (k - 1) * (m + 1))) * In)) * mu +
                        A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
            Gc.block(0, i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
        }
        Gn[k] = Gc;
        gn[k] = gc;
        if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd h3_ijk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int m, double& lscf) {
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_m = ArrayXXd::Zero(n, m + 1);
    ArrayXXd Gc = ArrayXXd::Zero(n, (m + 1));
    ArrayXXd gc = ArrayXXd::Zero(n, (m + 1));
    ArrayXXd tG(n, n);
    ArrayXXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_m;
    ArrayXXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_m;
    ArrayXXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_m;
    ArrayXXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_m;
    for(int i = 1; i <= m; i++) {
        tG = A1 * (dks(i - 1, 0) + Gc.col(i - 1));
        gc.col(i) = (tG - Gc.col(i - 1) - dks(i - 1, 0)) * mu + A1 * gc.col(i - 1);
        Gc.col(i) = tG;
        dks(i, 0) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        tG = A2 * (dks(0, k - 1) + (Go[0]).col(0));
        gc.col(0) = (tG - (Go[0]).col(0) - (dks(0, k - 1))) * mu + A2 * (go[0]).col(0);
        Gc.col(0) = tG;
        dks(0, k) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k) + Gc.col((i - 1))) +
                 A2 * (dks(i, k - 1) + (Go[0]).col(i));
            gc.col(i) = (tG - Gc.col(i - 1) - (Go[0]).col(i)
                         - (dks(i - 1, k) + dks(i, k - 1))) * mu +
                        A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
            Gc.col(i) = tG;
            dks(i, k) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + (Go[j]).col(0)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(0));
            gc.col(0) = (tG - (Go[j]).col(0) - (Go[j - 1]).col(0) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
            Gc.col(0) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
            for(int i = 1; i <= (m - k); i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gc.col((i - 1))) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) + (Go[j]).col(i)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(i));
                gc.col(i) = (tG - Gc.col(i - 1) - (Go[j]).col(i) - (Go[j - 1]).col(i)
                             - ((dks(i - 1, (k - j) + j * (m + 1)) + dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))))) * mu +
                            A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                Gc.col(i) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + (Go[k - 1]).col(0));
        gc.col(0) = (tG - (Go[k - 1]).col(0) - (dks(0, (k - 1) * (m + 1)))) * mu + A3 * (go[k - 1]).col(0);
        Gc.col(0) = tG;
        dks(0, k * (m + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= (m - k); i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gc.col((i - 1))) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + (Go[k - 1]).col(i));
            gc.col(i) = (tG - Gc.col(i - 1) - (Go[k - 1]).col(i)
                         - (dks(i - 1, k * (m + 1)) + dks(i, (k - 1) * (m + 1)))) * mu +
                        A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
            Gc.col(i) = tG;
            dks(i, k * (m + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
        }
        Gn[k] = Gc;
        gn[k] = gc;
        if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd& mu, const int m, const int p, double& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_np = MatrixXd::Zero(n, n * (p + 1));
    const MatrixXd zeromat_n_p = MatrixXd::Zero(n, p + 1);
    MatrixXd Gc = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd gc = MatrixXd::Zero(n, (p + 1));
    MatrixXd tG(n, n);
    MatrixXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_np;
    MatrixXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_p;
    MatrixXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_np;
    MatrixXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gc.block(0, (i - 1) * n, n, n));
        gc.col(i) = Gc.block(0, i * n, n, n) * mu + A1 * gc.col(i - 1);
        dks(i, 0) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        tG = A2 * (dks(0, k - 1) * In + (Go[0]).block(0, 0, n, n));
        gc.col(0) = (tG - (Go[0]).block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * (go[0]).col(0);
        Gc.block(0, 0, n, n) = tG;
        dks(0, k) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) * In + Gc.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + (Go[0]).block(0, i * n, n, n));
            gc.col(i) = (tG - (Go[0]).block(0, i * n, n, n)
                         - (dks(i, (k - 1)) * In)) * mu +
                        A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
            Gc.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, 0, n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, 0, n, n));
            gc.col(0) = (tG - (Go[j]).block(0, 0, n, n) - (Go[j - 1]).block(0, 0, n, n) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
            Gc.block(0, 0, n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, i * n, n, n)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, i * n, n, n));
                gc.col(i) = (tG - (Go[j]).block(0, i * n, n, n) - (Go[j - 1]).block(0, i * n, n, n)
                             - ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                            A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                Gc.block(0, i * n, n, n) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        tG = A3 * (dks(0, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, 0, n, n));
        gc.col(0) = (tG - (Go[k - 1]).col(0) - (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * (go[k - 1]).col(0);
        Gc.block(0, 0, n, n) = tG;
        dks(0, k * (m + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, i * n, n, n));
            gc.col(i) = (tG - (Go[k - 1]).block(0, i * n, n, n)
                             - (dks(i, (k - 1) * (m + 1)) * In)) * mu +
                            A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
            Gc.block(0, i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
        }
        Gn[k] = Gc;
        gn[k] = gc;
        if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd htil3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int m, const int p, double& lscf) {
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_p = ArrayXXd::Zero(n, p + 1);
    ArrayXXd Gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXXd gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXd tG(n);
    ArrayXXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_p;
    ArrayXXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_p;
    ArrayXXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_p;
    ArrayXXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.col(i) = A1 * (dks(i - 1, 0) + Gc.col(i - 1));
        gc.col(i) = Gc.col(i) * mu + A1 * gc.col(i - 1);
        dks(i, 0) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        tG = A2 * (dks(0, k - 1) + (Go[0]).col(0));
        gc.col(0) = (tG - (Go[0]).col(0) - (dks(0, k - 1))) * mu + A2 * (go[0]).col(0);
        Gc.col(0) = tG;
        dks(0, k) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + Gc.col((i - 1))) +
                 A2 * (dks(i, k - 1) + (Go[0]).col(i));
            gc.col(i) = (tG - (Go[0]).col(i)
                         - (dks(i, (k - 1)))) * mu +
                        A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
            Gc.col(i) = tG;
            dks(i, k) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + (Go[j]).col(0)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(0));
            gc.col(0) = (tG - (Go[j]).col(0) - (Go[j - 1]).col(0) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
            Gc.col(0) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gc.col((i - 1))) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) + (Go[j]).col(i)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(i));
                gc.col(i) = (tG - (Go[j]).col(i) - (Go[j - 1]).col(i)
                             - ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))))) * mu +
                            A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                Gc.col(i) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + (Go[k - 1]).col(0));
        gc.col(0) = (tG - (Go[k - 1]).col(0) - (dks(0, (k - 1) * (m + 1)))) * mu + A3 * (go[k - 1]).col(0);
        Gc.col(0) = tG;
        dks(0, k * (m + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gc.col((i - 1))) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + (Go[k - 1]).col(i));
            gc.col(i) = (tG - (Go[k - 1]).col(i)
                         - (dks(i, (k - 1) * (m + 1)))) * mu +
                        A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
            Gc.col(i) = tG;
            dks(i, k * (m + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
        }
        Gn[k] = Gc;
        gn[k] = gc;
        if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd& mu, const int m, const int p, double& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_np = MatrixXd::Zero(n, n * (p + 1));
    const MatrixXd zeromat_n_p = MatrixXd::Zero(n, p + 1);
    MatrixXd Gc = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd gc = MatrixXd::Zero(n, (p + 1));
    MatrixXd tG(n, n);
    MatrixXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_np;
    MatrixXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_p;
    MatrixXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_np;
    MatrixXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gc.block(0, (i - 1) * n, n, n));
        gc.col(i) = Gc.block(0, i * n, n, n) * mu + A1 * gc.col(i - 1);
        dks(i, 0) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        tG = A2 * (dks(0, k - 1) * In + (Go[0]).block(0, 0, n, n));
        gc.col(0) = (tG + (Go[0]).block(0, 0, n, n) + (dks(0, k - 1) * In)) * mu + A2 * (go[0]).col(0);
        Gc.block(0, 0, n, n) = tG;
        dks(0, k) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) * In + Gc.block(0, (i - 1) * n, n, n)) +
                 A2 * (dks(i, k - 1) * In + (Go[0]).block(0, i * n, n, n));
            gc.col(i) = (tG + (Go[0]).block(0, i * n, n, n)
                         + (dks(i, (k - 1)) * In)) * mu +
                        A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
            Gc.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, 0, n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, 0, n, n));
            gc.col(0) = (tG + (Go[j]).block(0, 0, n, n) + (Go[j - 1]).block(0, 0, n, n) + ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
            Gc.block(0, 0, n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + (Go[j]).block(0, i * n, n, n)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + (Go[j - 1]).block(0, i * n, n, n));
                gc.col(i) = (tG + (Go[j]).block(0, i * n, n, n) + (Go[j - 1]).block(0, i * n, n, n)
                             + ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                            A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                Gc.block(0, i * n, n, n) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        tG = A3 * (dks(0, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, 0, n, n));
        gc.col(0) = (tG + (Go[k - 1]).col(0) + (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * (go[k - 1]).col(0);
        Gc.block(0, 0, n, n) = tG;
        dks(0, k * (m + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) * In + (Go[k - 1]).block(0, i * n, n, n));
            gc.col(i) = (tG + (Go[k - 1]).block(0, i * n, n, n)
                         + (dks(i, (k - 1) * (m + 1)) * In)) * mu +
                        A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
            Gc.block(0, i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
        }
        Gn[k] = Gc;
        gn[k] = gc;
        if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd hhat3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int m, const int p, double& lscf) {
    const int n = A1.rows();
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_p = ArrayXXd::Zero(n, p + 1);
    ArrayXXd Gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXXd gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXd tG(n);
    ArrayXXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_p;
    ArrayXXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_p;
    ArrayXXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_p;
    ArrayXXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.col(i) = A1 * (dks(i - 1, 0) + Gc.col(i - 1));
        gc.col(i) = Gc.col(i) * mu + A1 * gc.col(i - 1);
        dks(i, 0) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        tG = A2 * (dks(0, k - 1) + (Go[0]).col(0));
        gc.col(0) = (tG + (Go[0]).col(0) + (dks(0, k - 1))) * mu + A2 * (go[0]).col(0);
        Gc.col(0) = tG;
        dks(0, k) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + Gc.col((i - 1))) +
                 A2 * (dks(i, k - 1) + (Go[0]).col(i));
            gc.col(i) = (tG + (Go[0]).col(i) + (dks(i, (k - 1)))) * mu +
                        A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
            Gc.col(i) = tG;
            dks(i, k) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + (Go[j]).col(0)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(0));
            gc.col(0) = (tG + (Go[j]).col(0) + (Go[j - 1]).col(0) + ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
            Gc.col(0) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                tG = A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gc.col((i - 1))) +
                     A2 * (dks(i, (k - j - 1) + j * (m + 1)) + (Go[j]).col(i)) +
                     A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + (Go[j - 1]).col(i));
                gc.col(i) = (tG + (Go[j]).col(i) + (Go[j - 1]).col(i)
                             + ((dks(i, (k - j - 1) + j * (m + 1)) + dks(i, (k - j) + (j - 1) * (m + 1))))) * mu +
                            A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                Gc.col(i) = tG;
                dks(i, (k - j) + j * (m + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + (Go[k - 1]).col(0));
        gc.col(0) = (tG + (Go[k - 1]).col(0) + (dks(0, (k - 1) * (m + 1)))) * mu + A3 * (go[k - 1]).col(0);
        Gc.col(0) = tG;
        dks(0, k * (m + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
        for(int i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gc.col((i - 1))) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + (Go[k - 1]).col(i));
            gc.col(i) = (tG + (Go[k - 1]).col(i)
                         + (dks(i, (k - 1) * (m + 1)))) * mu +
                        A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
            Gc.col(i) = tG;
            dks(i, k * (m + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
        }
        Gn[k] = Gc;
        gn[k] = gc;
        if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
            dks /= 1e10;
            for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
            for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
            lscf -= log(1e10);
        }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil3_pqr_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2, const Eigen::MatrixXd& A3,
                            const Eigen::VectorXd mu, const int p, const int q, const int r) { //, double& lscf) {
    const int n = A1.rows();
    const MatrixXd In = MatrixXd::Identity(n, n);
    const int m = q + r;
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (q + 1) * (r + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const MatrixXd zeromat_n_np = MatrixXd::Zero(n, n * (p + 1));
    const MatrixXd zeromat_n_p = MatrixXd::Zero(n, p + 1);
    MatrixXd Gc = MatrixXd::Zero(n, n * (p + 1));
    MatrixXd gc = MatrixXd::Zero(n, (p + 1));
    MatrixXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_np;
    MatrixXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_p;
    MatrixXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_np;
    MatrixXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + Gc.block(0, (i - 1) * n, n, n));
        gc.col(i) = Gc.block(0, i * n, n, n) * mu + A1 * gc.col(i - 1);
        dks(i, 0) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        if(k <= q) {
            Gc.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + (Go[0]).block(0, 0, n, n));
            gc.col(0) = Gc.block(0, 0, n, n) * mu + A2 * (go[0]).col(0);
            dks(0, k) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + Gc.block(0, (i - 1) * n, n, n)) +
                A2 * (dks(i, k - 1) * In + (Go[0]).block(0, i * n, n, n));
                gc.col(i) = Gc.block(0, i * n, n, n) * mu +
                A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
                dks(i, k) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
            }
        } else {
            Gc = zeromat_n_np;
            gc = zeromat_n_p;
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            if(((k - j) <= q) && (j <= r)) {
                Gc.block(0, 0, n, n) = A2 * (dks(0, (k - j - 1) + j * (q + 1)) * In + (Go[j]).block(0, 0, n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (q + 1)) * In + (Go[j - 1]).block(0, 0, n, n));
                gc.col(0) = Gc.block(0, 0, n, n) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
                dks(0, (k - j) + j * (q + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
                for(int i = 1; i <= p; i++) {
                    Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, (k - j) + j * (q + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                    A2 * (dks(i, (k - j - 1) + j * (q + 1)) * In + (Go[j]).block(0, i * n, n, n)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (q + 1)) * In + (Go[j - 1]).block(0, i * n, n, n));
                    gc.col(i) = Gc.block(0, i * n, n, n) * mu +
                    A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                    dks(i, (k - j) + j * (q + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
                }
            } else {
                Gc = zeromat_n_np;
                gc = zeromat_n_p;
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        if(k <= r) {
            Gc.block(0, 0, n, n) = A3 * (dks(0, (k - 1) * (q + 1)) * In + (Go[k - 1]).block(0, 0, n, n));
            gc.col(0) = Gc.block(0, 0, n, n) * mu + A3 * (go[k - 1]).col(0);
            dks(0, k * (q + 1)) = (Gc.block(0, 0, n, n).trace() + gc.col(0).dot(mu)) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gc.block(0, i * n, n, n) = A1 * (dks(i - 1, k * (q + 1)) * In + Gc.block(0, (i - 1) * n, n, n)) +
                A3 * (dks(i, (k - 1) * (q + 1)) * In + (Go[k - 1]).block(0, i * n, n, n));
                gc.col(i) = Gc.block(0, i * n, n, n) * mu +
                A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
                dks(i, k * (q + 1)) = (Gc.block(0, i * n, n, n).trace() + gc.col(i).dot(mu)) / (2 * (k + i));
            }
        } else {
            Gc = zeromat_n_np;
            gc = zeromat_n_p;
        }
        Gn[k] = Gc;
        gn[k] = gc;
        // if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
        //     dks /= 1e10;
        //     for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
        //     for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
        //     lscf -= log(1e10);
        // }
    }
    return dks;
}

// // [[Rcpp::export]]
Eigen::ArrayXXd dtil3_pqr_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2, const Eigen::ArrayXd& A3,
                            const Eigen::ArrayXd& mu, const int p, const int q, const int r) { //, double& lscf) {
    const int n = A1.rows();
    const int m = q + r;
    ArrayXXd dks = ArrayXXd::Zero(p + 1, (q + 1) * (r + 1));
    dks(0, 0) = 1;
    double thr = std::numeric_limits<double>::max() / 100 / double(n);
    const ArrayXXd zeromat_n_p = ArrayXXd::Zero(n, p + 1);
    ArrayXXd Gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXXd gc = ArrayXXd::Zero(n, (p + 1));
    ArrayXXd Go[m];
    for(int ii = 0; ii < m; ii++) Go[ii] = zeromat_n_p;
    ArrayXXd go[m];
    for(int ii = 0; ii < m; ii++) go[ii] = zeromat_n_p;
    ArrayXXd Gn[m + 1];
    for(int ii = 0; ii <= m; ii++) Gn[ii] = zeromat_n_p;
    ArrayXXd gn[m + 1];
    for(int ii = 0; ii <= m; ii++) gn[ii] = zeromat_n_p;
    for(int i = 1; i <= p; i++) {
        Gc.col(i) = A1 * (dks(i - 1, 0) + Gc.col(i - 1));
        gc.col(i) = Gc.col(i) * mu + A1 * gc.col(i - 1);
        dks(i, 0) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * i);
    }
    Gn[0] = Gc;
    gn[0] = gc;
    for(int k = 1; k <= m; k++) {
        for(int ii = 0; ii < k; ii++) Go[ii] = Gn[ii];
        for(int ii = 0; ii < k; ii++) go[ii] = gn[ii];
        if(k <= q) {
            Gc.col(0) = A2 * (dks(0, k - 1) + (Go[0]).col(0));
            gc.col(0) = Gc.col(0) * mu + A2 * (go[0]).col(0);
            dks(0, k) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gc.col(i) = A1 * (dks(i - 1, k) + Gc.col((i - 1))) +
                A2 * (dks(i, k - 1) + (Go[0]).col(i));
                gc.col(i) = Gc.col(i) * mu +
                A1 * gc.col(i - 1) + A2 * (go[0]).col(i);
                dks(i, k) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
            }
        } else {
            Gc = zeromat_n_p;
            gc = zeromat_n_p;
        }
        Gn[0] = Gc;
        gn[0] = gc;
        for(int j = 1; j < k; j++) {
            if(((k - j) <= q) && (j <= r)) {
                Gc.col(0) = A2 * (dks(0, (k - j - 1) + j * (q + 1)) + (Go[j]).col(0)) + A3 * (dks(0, (k - j) + (j - 1) * (q + 1)) + (Go[j - 1]).col(0));
                gc.col(0) = Gc.col(0) * mu + A2 * (go[j]).col(0) + A3 * (go[j - 1]).col(0);
                dks(0, (k - j) + j * (q + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
                for(int i = 1; i <= p; i++) {
                    Gc.col(i) = A1 * (dks(i - 1, (k - j) + j * (q + 1)) + Gc.col((i - 1))) +
                    A2 * (dks(i, (k - j - 1) + j * (q + 1)) + (Go[j]).col(i)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (q + 1)) + (Go[j - 1]).col(i));
                    gc.col(i) = Gc.col(i) * mu +
                    A1 * gc.col(i - 1) + A2 * (go[j]).col(i) + A3 * (go[j - 1]).col(i);
                    dks(i, (k - j) + j * (q + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
                }
            } else {
                Gc = zeromat_n_p;
                gc = zeromat_n_p;
            }
            Gn[j] = Gc;
            gn[j] = gc;
        }
        if(k <= r) {
            Gc.col(0) = A3 * (dks(0, (k - 1) * (q + 1)) + (Go[k - 1]).col(0));
            gc.col(0) = Gc.col(0) * mu + A3 * (go[k - 1]).col(0);
            dks(0, k * (q + 1)) = (Gc.col(0).sum() + (mu * gc.col(0)).sum()) / (2 * k);
            for(int i = 1; i <= p; i++) {
                Gc.col(i) = A1 * (dks(i - 1, k * (q + 1)) + Gc.col((i - 1))) +
                A3 * (dks(i, (k - 1) * (q + 1)) + (Go[k - 1]).col(i));
                gc.col(i) = Gc.col(i) * mu +
                A1 * gc.col(i - 1) + A3 * (go[k - 1]).col(i);
                dks(i, k * (q + 1)) = (Gc.col(i).sum() + (mu * gc.col(i)).sum()) / (2 * (k + i));
            }
        } else {
            Gc = zeromat_n_p;
            gc = zeromat_n_p;
        }
        Gn[k] = Gc;
        gn[k] = gc;
        // if(Gc.maxCoeff() > thr || gc.maxCoeff() > thr) {
        //     dks /= 1e10;
        //     for(int ii = 0; ii <= k; ii++) Gn[ii] /= 1e10;
        //     for(int ii = 0; ii <= k; ii++) gn[ii] /= 1e10;
        //     lscf -= log(1e10);
        // }
    }
    return dks;
}
