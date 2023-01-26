#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <limits>
#include "dk_funs.h"

#ifdef _OPENMP
    #include <omp.h>
    // [[Rcpp::plugins(openmp)]]
#endif

#define LN1E10 M_LN10 * 10

using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix;
using Eigen::Array;
using Eigen::MatrixBase;
using Eigen::ArrayBase;
using Eigen::Dynamic;
using Eigen::Index;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXl;
typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> VectorXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;


// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d1_i_vE(const Eigen::ArrayBase<Derived>& L, const Eigen::Index m, 
        Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
        const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    Index n = L.size();
    ArrayXx dks = ArrayXx::Zero(m + 1);
    dks(0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx uk = ArrayXx::Zero(n);
    for(Index k = 1; k <= m; k++) {
        uk = L * (dks(k - 1) + uk);
        dks(k) = uk.sum() / (2 * k);
        if(uk.maxCoeff() > thr) {
            dks(k) /= 1e10;
            uk /= 1e10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXd d1_i_vE(const ArrayBase<ArrayXd>& L, const Index m,
                         ArrayXd &lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d1_i_mE(const Eigen::MatrixBase<Derived>& A, const Eigen::Index m,
        Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
        const typename Derived::Scalar thr_margin) {
    Eigen::SelfAdjointEigenSolver<Derived> eigA(A, Eigen::EigenvaluesOnly);
    Array<typename Derived::Scalar, Dynamic, 1> L = eigA.eigenvalues();
    return d1_i_vE(L, m, lscf, thr_margin);
}
template ArrayXd d1_i_mE(const MatrixBase<MatrixXd>& A, const Index m,
                         ArrayXd& lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
dtil1_i_vE(const Eigen::ArrayBase<Derived>& L,
           const Eigen::ArrayBase<Derived>& mu, const Eigen::Index m,
           Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
           const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = L.size();
    ArrayXx D = square(mu);
    ArrayXx dks = ArrayXx::Zero(m + 1);
    dks(0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx uk = ArrayXx::Zero(n);
    ArrayXx vk = ArrayXx::Zero(n);
    for(Index k = 1; k <= m; k++) {
        uk = L * (dks(k - 1) + uk);
        vk = D * uk + L * vk;
        dks(k) = (uk + vk).sum() / (2 * k);
        if(uk.maxCoeff() > thr || vk.maxCoeff() > thr) {
            dks(k) /= 1e10;
            uk /= 1e10;
            vk /= 1e10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXd dtil1_i_vE(const ArrayBase<ArrayXd>& L,
                            const ArrayBase<ArrayXd>& mu, const Index m,
                            ArrayXd& lscf, const double thr_margin);
template ArrayXl dtil1_i_vE(const ArrayBase<ArrayXl>& L,
                            const ArrayBase<ArrayXl>& mu, const Index m,
                            ArrayXl& lscf, const long double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
dtil1_i_mE(const Eigen::MatrixBase<Derived>& A,
           const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
           const Eigen::Index m,
           Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
           const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    Eigen::SelfAdjointEigenSolver<MatrixXx> eigA(A);
    ArrayXx L = eigA.eigenvalues();
    ArrayXx mud = eigA.eigenvectors().transpose() * mu;
    return dtil1_i_vE(L, mud, m, lscf, thr_margin);
}
template ArrayXd dtil1_i_mE(const MatrixBase<MatrixXd>& A,
                            const VectorXd& mu, const Index m,
                            ArrayXd& lscf, const double thr_margin);
template ArrayXl dtil1_i_mE(const MatrixBase<MatrixXl>& A,
                            const VectorXl& mu, const Index m,
                            ArrayXl& lscf, const long double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
arl_vE(const Eigen::ArrayBase<Derived>& L, const Eigen::ArrayBase<Derived>& D,
       const Eigen::Index m, const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = L.size();
    ArrayXx lscf = ArrayXx::Zero(m + 1);
    ArrayXXx arls = ArrayXXx::Zero(m + 1, m + 1);
    arls.col(0) = d1_i_vE(L, m, lscf, thr_margin);
    ArrayXXx wrls = ArrayXXx::Zero(n, m);
    for(Index k = 0; k < m; k++) {
        for(Index l = 0; l <= k; l++) {
            wrls.col(l) = L * (arls(k, l) + wrls.col(l));
            arls(k + 1, l + 1) = (D * wrls.col(l)).sum();
        }
    }
    return arls;
}
template ArrayXXd arl_vE(const ArrayBase<ArrayXd>& L,
                         const ArrayBase<ArrayXd>& D, const Index m,
                         const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
arl_mE(const Eigen::MatrixBase<Derived>& A,
       const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
       const Eigen::Index m, const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    Eigen::SelfAdjointEigenSolver<Derived> eigA(A);
    ArrayXx L = eigA.eigenvalues();
    ArrayXx mud = eigA.eigenvectors().transpose() * mu;
    ArrayXx D = square(mud);
    return arl_vE(L, D, m, thr_margin);
}
template ArrayXXd arl_mE(const MatrixBase<MatrixXd>& A, const VectorXd& mu,
                         const Index m, const double thr_margin);


// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_pj_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const Eigen::Index m, const Eigen::Index p,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx G_k_i = MatrixXx::Zero(n, n * (p + 1));
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(G_k_i.block(0, i * n, n, n).data(), n, n).noalias() =
            A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        dks(i, 0) = G_k_i.block(0, i * n, n, n).trace() / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        G_k_i.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        dks(0, k) = G_k_i.block(0, 0, n, n).trace() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                                          A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            dks(i, k) = G_k_i.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            //lscf.rightCols(m + 1 - k) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d2_pj_mE(const MatrixBase<MatrixXd>& A1,
                           const MatrixBase<MatrixXd>& A2,
                           const Index m, const Index p,
                           ArrayXd& lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_pj_vE(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::Index m, const Eigen::Index p,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.size();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx G_k_i = ArrayXXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        dks(i, 0) = G_k_i.col(i).sum() / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        G_k_i.col(0) = A2 * (dks(0, k - 1) + G_k_i.col(0));
        dks(0, k) = G_k_i.col(0).sum() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            G_k_i.col(i) = A1 * (dks(i - 1, k) + G_k_i.col(i - 1)) +
                           A2 * (dks(i, k - 1) + G_k_i.col(i));
            dks(i, k) = G_k_i.col(i).sum() / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            //lscf.rightCols(m + 1 - k) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d2_pj_vE(const ArrayBase<ArrayXd>& A1,
                           const ArrayBase<ArrayXd>& A2,
                           const Index m, const Index p,
                           ArrayXd& lscf, const double thr_margin);


// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx Go = MatrixXx::Zero(n, n * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 1));
    for(Index k = 1; k <= m; k++) {
        if(k % 200 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k) = Gn.block(0, 0, n, n * k);

        MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =
            A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        dks(0, k) = Gn.block(0, 0, n, n).trace() / (2 * k);
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            MatrixXx::Map(Gn.block(0, n * i1, n, n).data(), n, n).noalias() =
                A1 * (dks(i1 - 1, k - i1) * In + Go.block(0, n * (i1 - 1), n, n)) +
                A2 * (dks(i1, k - i1 - 1) * In + Go.block(0, n * i1, n, n));
            dks(i1, k - i1) = Gn.block(0, n * i1, n, n).trace() / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        MatrixXx::Map(Gn.block(0, n * k, n, n).data(), n, n).noalias() =
            A1 * (dks(k - 1, 0) * In + Go.block(0, n * (k - 1), n, n));
        dks(k, 0) = Gn.block(0, n * k, n, n).trace() / (2 * k);

        if(Gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) dks(i1, k - i1) /= 1e10;
            Gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = std::max(k - i1, 0); i2 <= m - i1; i2++) {
            //         lscf(i1, i2) -= LN1E10;
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d2_ij_mE(const MatrixBase<MatrixXd>& A1,
                           const MatrixBase<MatrixXd>& A2,
                           const Index m, ArrayXd& lscf,
                           const double thr_margin, int nthreads);
template ArrayXXl d2_ij_mE(const MatrixBase<MatrixXl>& A1,
                           const MatrixBase<MatrixXl>& A2,
                           const Index m, ArrayXl& lscf,
                           const long double thr_margin, int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2,
         const Eigen::Index m, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx Go = ArrayXXx::Zero(n, m);
    ArrayXXx Gn = ArrayXXx::Zero(n, m + 1);
    for(Index k = 1; k <= m; k++) {
        if(k % 2000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k) = Gn.block(0, 0, n, k);

        Gn.col(0) = A2 * (dks(0, k - 1) + Go.col(0));
        dks(0, k) = Gn.col(0).sum() / (2 * k);
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            Gn.col(i1) = A1 * (dks(i1 - 1, k - i1) + Go.col(i1 - 1)) +
                         A2 * (dks(i1, k - i1 - 1) + Go.col(i1));
            dks(i1, k - i1) = Gn.col(i1).sum() / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        Gn.col(k) = A1 * (dks(k - 1, 0) + Go.col(k - 1));
        dks(k, 0) = Gn.col(k).sum() / (2 * k);

        if(Gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) dks(i1, k - i1) /= 1e10;
            Gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = std::max(k - i1, 0); i2 <= m - i1; i2++) {
            //         lscf(i1, i2) -= LN1E10;
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d2_ij_vE(const ArrayBase<ArrayXd>& A1,
                           const ArrayBase<ArrayXd>& A2,
                           const Index m, ArrayXd& lscf,
                           const double thr_margin, int nthreads);
template ArrayXXl d2_ij_vE(const ArrayBase<ArrayXl>& A1,
                           const ArrayBase<ArrayXl>& A2,
                           const Index m, ArrayXl& lscf,
                           const long double thr_margin, int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 1));
    MatrixXx go = MatrixXx::Zero(n, m);
    MatrixXx gn = MatrixXx::Zero(n, m + 1);
    for(Index k = 1; k <= m; k++) {
        if(k % 200 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k) = Gn.block(0, 0, n, n * k);
        go.block(0, 0, n, k) = gn.block(0, 0, n, k);

        tG.noalias() = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG - Go.block(0, 0, n, n)
             - ((dks(0, k - 1)) * In)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            tG.noalias() = A1 * (dks(i1 - 1, k - i1) * In + Go.block(0, n * (i1 - 1), n, n)) +
                           A2 * (dks(i1, k - i1 - 1) * In + Go.block(0, n * i1, n, n));
            MatrixXx::Map(gn.col(i1).data(), n, 1).noalias() =
                (tG - Go.block(0, n * (i1 - 1), n, n) -
                 Go.block(0, n * i1, n, n)
                 - ((dks(i1 - 1, k - i1) + dks(i1, k - i1 - 1)) * In)) * mu +
                A1 * go.col(i1 - 1) + A2 * go.col(i1);
            Gn.block(0, n * i1, n, n) = tG;
            dks(i1, k - i1) = (Gn.block(0, n * i1, n, n).trace() + gn.col(i1).dot(mu)) / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        tG.noalias() = A1 * (dks(k - 1, 0) * In + Go.block(0, n * (k - 1), n, n));
        MatrixXx::Map(gn.col(k).data(), n, 1).noalias() =
            (tG - Go.block(0, n * (k - 1), n, n)
             - ((dks(k - 1, 0)) * In)) * mu + A1 * go.col(k - 1);
        Gn.block(0, n * k, n, n) = tG;
        dks(k, 0) = (Gn.block(0, n * k, n, n).trace() + gn.col(k).dot(mu)) / (2 * k);

        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) dks(i1, k - i1) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = std::max(k - i1, 0); i2 <= m - i1; i2++) {
            //         lscf(i1, i2) -= LN1E10;
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd h2_ij_mE(const MatrixBase<MatrixXd>& A1,
                           const MatrixBase<MatrixXd>& A2,
                           const VectorXd& mu, const Index m,
                           ArrayXd& lscf, const double thr_margin,
                           int nthreads);
template ArrayXXl h2_ij_mE(const MatrixBase<MatrixXl>& A1,
                           const MatrixBase<MatrixXl>& A2,
                           const VectorXl& mu, const Index m,
                           ArrayXl& lscf, const long double thr_margin,
                           int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_vE(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::ArrayBase<Derived>& mu,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(m + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, m);
    ArrayXXx Gn = ArrayXXx::Zero(n, m + 1);
    ArrayXXx go = ArrayXXx::Zero(n, m);
    ArrayXXx gn = ArrayXXx::Zero(n, m + 1);
    for(Index k = 1; k <= m; k++) {
        if(k % 2000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k) = Gn.block(0, 0, n, k);
        go.block(0, 0, n, k) = gn.block(0, 0, n, k);

        tG = A2 * (dks(0, k - 1) + Go.col(0));
        gn.col(0) = (tG - Go.col(0)
             - ((dks(0, k - 1)))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            tG = A1 * (dks(i1 - 1, k - i1) + Go.col(i1 - 1)) +
                 A2 * (dks(i1, k - i1 - 1) + Go.col(i1));
            gn.col(i1) =
                (tG - Go.col(i1 - 1) -
                 Go.col(i1)
                 - ((dks(i1 - 1, k - i1) + dks(i1, k - i1 - 1)))) * mu +
                A1 * go.col(i1 - 1) + A2 * go.col(i1);
            Gn.col(i1) = tG;
            dks(i1, k - i1) = (Gn.col(i1).sum() + (mu * gn.col(i1)).sum()) / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        tG = A1 * (dks(k - 1, 0) + Go.col(k - 1));
        gn.col(k) = (tG - Go.col(k - 1)
             - ((dks(k - 1, 0)))) * mu + A1 * go.col(k - 1);
        Gn.col(k) = tG;
        dks(k, 0) = (Gn.col(k).sum() + (mu * gn.col(k)).sum()) / (2 * k);

        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) dks(i1, k - i1) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = std::max(k - i1, 0); i2 <= m - i1; i2++) {
            //         lscf(i1, i2) -= LN1E10;
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd h2_ij_vE(const ArrayBase<ArrayXd>& A1,
                           const ArrayBase<ArrayXd>& A2,
                           const ArrayBase<ArrayXd>& mu,
                           const Index m, ArrayXd& lscf,
                           const double thr_margin, int nthreads);
template ArrayXXl h2_ij_vE(const ArrayBase<ArrayXl>& A1,
                           const ArrayBase<ArrayXl>& A2,
                           const ArrayBase<ArrayXl>& mu,
                           const Index m, ArrayXl& lscf,
                           const long double thr_margin, int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil2_pj_mE(const Eigen::MatrixBase<Derived>& A1,
            const Eigen::MatrixBase<Derived>& A2,
            const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
            const Eigen::Index m, const Eigen::Index p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx G_k_i = MatrixXx::Zero(n, n * (p + 1));
    MatrixXx g_k_i = MatrixXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        g_k_i.col(i) = G_k_i.block(0, i * n, n, n) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        tG.noalias() = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = (tG - G_k_i.block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * g_k_i.col(0);
        G_k_i.block(0, 0, n, n) = tG;
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                           A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            g_k_i.col(i) = (tG - G_k_i.block(0, i * n, n, n) // require aliasing
                            - (dks(i, k - 1) * In)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.block(0, i * n, n, n) = tG;
            dks(i, k) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            //lscf.rightCols(m + 1 - k) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd htil2_pj_mE(const MatrixBase<MatrixXd>& A1, 
                              const MatrixBase<MatrixXd>& A2,
                              const VectorXd& mu, const Index m, const Index p,
                              ArrayXd& lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil2_pj_vE(const Eigen::ArrayBase<Derived>& A1,
            const Eigen::ArrayBase<Derived>& A2,
            const Eigen::ArrayBase<Derived>& mu, const Eigen::Index m, const Eigen::Index p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.size();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx G_k_i = ArrayXXx::Zero(n, p + 1);
    ArrayXXx g_k_i = ArrayXXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = G_k_i.col(i) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        tG = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = (tG - G_k_i.col(0) - dks(0, k - 1)) * mu + A2 * g_k_i.col(0);
        G_k_i.col(0) = tG;
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
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
            //lscf.rightCols(m + 1 - k) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd htil2_pj_vE(const ArrayBase<ArrayXd>& A1,
                              const ArrayBase<ArrayXd>& A2,
                              const ArrayBase<ArrayXd>& mu,
                              const Index m, const Index p,
                              ArrayXd& lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat2_pj_mE(const Eigen::MatrixBase<Derived>& A1,
            const Eigen::MatrixBase<Derived>& A2,
            const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
            const Eigen::Index m, const Eigen::Index p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx G_k_i = MatrixXx::Zero(n, n * (p + 1));
    MatrixXx g_k_i = MatrixXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(G_k_i.block(0, i * n, n, n).data(), n, n).noalias() =
            A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        MatrixXx::Map(g_k_i.col(i).data(), n, 1).noalias() =
            G_k_i.block(0, i * n, n, n) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        tG.noalias() = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = (tG + G_k_i.block(0, 0, n, n) + (dks(0, k - 1) * In)) * mu + A2 * g_k_i.col(0);
        G_k_i.block(0, 0, n, n) = tG;
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * (dks(i - 1, k) * In + G_k_i.block(0, (i - 1) * n, n, n)) +
                           A2 * (dks(i, k - 1) * In + G_k_i.block(0, i * n, n, n));
            g_k_i.col(i) = (tG + G_k_i.block(0, i * n, n, n) // require aliasing
                            + (dks(i, k - 1) * In)) * mu +
                           A1 * g_k_i.col(i - 1) + A2 * g_k_i.col(i);
            G_k_i.block(0, i * n, n, n) = tG;
            dks(i, k) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * (k + i));
        }
        if(G_k_i.maxCoeff() > thr || g_k_i.maxCoeff() > thr) {
            dks.col(k) /= 1e10;
            G_k_i /= 1e10;
            g_k_i /= 1e10;
            //lscf.rightCols(m + 1 - k) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd hhat2_pj_mE(const MatrixBase<MatrixXd>& A1, 
                              const MatrixBase<MatrixXd>& A2,
                              const VectorXd& mu,
                              const Index m, const Index p,
                              ArrayXd& lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat2_pj_vE(const Eigen::ArrayBase<Derived>& A1,
            const Eigen::ArrayBase<Derived>& A2,
            const Eigen::ArrayBase<Derived>& mu,
            const Eigen::Index m, const Eigen::Index p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.size();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, m + 1);
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx G_k_i = ArrayXXx::Zero(n, p + 1);
    ArrayXXx g_k_i = ArrayXXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = G_k_i.col(i) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        tG = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = (tG + G_k_i.col(0) + dks(0, k - 1)) * mu + A2 * g_k_i.col(0);
        G_k_i.col(0) = tG;
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
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
            //lscf.rightCols(m + 1 - k) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd hhat2_pj_vE(const ArrayBase<ArrayXd>& A1,
                              const ArrayBase<ArrayXd>& A2,
                              const ArrayBase<ArrayXd>& mu,
                              const Index m, const Index p,
                              ArrayXd& lscf, const double thr_margin);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil2_pq_mE(const Eigen::MatrixBase<Derived>& A1,
            const Eigen::MatrixBase<Derived>& A2,
            const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
            const Eigen::Index p, const Eigen::Index q) { //, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf, const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, q + 1);
    dks(0, 0) = 1;
    // Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx G_k_i = MatrixXx::Zero(n, n * (p + 1));
    MatrixXx g_k_i = MatrixXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        G_k_i.block(0, i * n, n, n) = A1 * (dks(i - 1, 0) * In + G_k_i.block(0, (i - 1) * n, n, n));
        g_k_i.col(i) = G_k_i.block(0, i * n, n, n) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.block(0, i * n, n, n).trace() + g_k_i.col(i).dot(mu)) / (2 * i);
    }
    for(Index k = 1; k <= q; k++) {
        if(k % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        G_k_i.block(0, 0, n, n) = A2 * (dks(0, k - 1) * In + G_k_i.block(0, 0, n, n));
        g_k_i.col(0) = G_k_i.block(0, 0, n, n) * mu + A2 * g_k_i.col(0);
        dks(0, k) = (G_k_i.block(0, 0, n, n).trace() + g_k_i.col(0).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
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
        //     lscf -= LN1E10;
        // }
    }
    return dks;
}
template ArrayXXd dtil2_pq_mE(const MatrixBase<MatrixXd>& A1,
                              const MatrixBase<MatrixXd>& A2,
                              const VectorXd& mu, const Index p, const Index q);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil2_pq_vE(const Eigen::ArrayBase<Derived>& A1,
            const Eigen::ArrayBase<Derived>& A2,
            const Eigen::ArrayBase<Derived>& mu,
            const Eigen::Index p, const Eigen::Index q) { //, EEigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf, const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.size();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, q + 1);
    dks(0, 0) = 1;
    // Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx G_k_i = ArrayXXx::Zero(n, p + 1);
    ArrayXXx g_k_i = ArrayXXx::Zero(n, p + 1);
    for(Index i = 1; i <= p; i++) {
        G_k_i.col(i) = A1 * (dks(i - 1, 0) + G_k_i.col(i - 1));
        g_k_i.col(i) = G_k_i.col(i) * mu + A1 * g_k_i.col(i - 1);
        dks(i, 0) = (G_k_i.col(i).sum() + (g_k_i.col(i) * mu).sum()) / (2 * i);
    }
    for(Index k = 1; k <= q; k++) {
        if(k % 1000000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        G_k_i.col(0) = A2 * (dks(0, k - 1) + G_k_i.col(0));
        g_k_i.col(0) = G_k_i.col(0) * mu + A2 * g_k_i.col(0);
        dks(0, k) = (G_k_i.col(0).sum() + (g_k_i.col(0) * mu).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
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
        //     lscf -= LN1E10;
        // }
    }
    return dks;
}
template ArrayXXd dtil2_pq_vE(const ArrayBase<ArrayXd>& A1,
                              const ArrayBase<ArrayXd>& A2,
                              const ArrayBase<ArrayXd>& mu,
                              const Index p, const Index q);


// This is a utility function used in d/h3_ijk_*E() for indexing of
// Gn and the like, which represent triangular grids of order k:
//    0 ... k (i2)
//  0 x x x x
//  . x x x o
//  . x x o o
//  k x o o o
// (i1)
// This function returns the index of (i1, i2)-th x.
inline Eigen::Index id3(Eigen::Index i1, Eigen::Index i2, Eigen::Index k) {
    return i1 + (k + 1) * i2 - (i2 - 1) * i2 / 2;
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx Go = MatrixXx::Zero(n, n * (m + 1) * m / 2);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 2) * (m + 1) / 2);
    for(Index k = 1; k <= m; k++) {
        if(k % 50 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k * (k + 1) / 2) = Gn.block(0, 0, n, n * k * (k + 1) / 2);

        MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =
            A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, 0, n, n));
        dks(0, k * (m + 1)) = Gn.block(0, 0, n, n).trace() / (2 * k);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            MatrixXx::Map(Gn.block(0, n * id3(0, i2, k), n, n).data(), n, n).noalias() =
                A2 * (dks(0, i2 - 1 + i3 * (m + 1)) * In + Go.block(0, n * id3(0, i2 - 1, k - 1), n, n)) +
                A3 * (dks(0, i2 + (i3 - 1) * (m + 1)) * In + Go.block(0, n * id3(0, i2, k - 1), n, n));
            dks(0, i2 + i3 * (m + 1)) = Gn.block(0, n * id3(0, i2, k), n, n).trace() / (2 * k);
        }
        MatrixXx::Map(Gn.block(0, n * id3(0, k, k), n, n).data(), n, n).noalias() =
            A2 * (dks(0, k - 1) * In + Go.block(0, n * id3(0, k - 1, k - 1), n, n));
        dks(0, k) = Gn.block(0, n * id3(0, k, k), n, n).trace() / (2 * k);
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            MatrixXx::Map(Gn.block(0, n * i1, n, n).data(), n, n).noalias() =
                A1 * (dks(i1 - 1, (k - i1) * (m + 1)) * In + Go.block(0, n * (i1 - 1), n, n)) +
                A3 * (dks(i1, (k - i1 - 1) * (m + 1)) * In + Go.block(0, n * i1, n, n));
            dks(i1, (k - i1) * (m + 1)) = Gn.block(0, n * i1, n, n).trace() / (2 * k);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                MatrixXx::Map(Gn.block(0, n * id3(i1, i2, k), n, n).data(), n, n).noalias() =
                    A1 * (dks(i1 - 1, i2 + i3 * (m + 1)) * In + Go.block(0, n * id3(i1 - 1, i2, k - 1), n, n)) +
                    A2 * (dks(i1, i2 - 1 + i3 * (m + 1)) * In + Go.block(0, n * id3(i1, i2 - 1, k - 1), n, n)) +
                    A3 * (dks(i1, i2 + (i3 - 1) * (m + 1)) * In + Go.block(0, n * id3(i1, i2, k - 1), n, n));
                dks(i1, i2 + i3 * (m + 1)) = Gn.block(0, n * id3(i1, i2, k), n, n).trace() / (2 * k);
            }
            MatrixXx::Map(Gn.block(0, n * id3(i1, k - i1, k), n, n).data(), n, n).noalias() =
                A1 * (dks(i1 - 1, k - i1) * In + Go.block(0, n * id3(i1 - 1, k - i1, k - 1), n, n)) +
                A2 * (dks(i1, k - i1 - 1) * In + Go.block(0, n * id3(i1, k - i1 - 1, k - 1), n, n));
            dks(i1, k - i1) = Gn.block(0, n * id3(i1, k - i1, k), n, n).trace() / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        //
        MatrixXx::Map(Gn.block(0, n * k, n, n).data(), n, n).noalias() =
            A1 * (dks(k - 1, 0) * In + Go.block(0, n * (k - 1), n, n));
        dks(k, 0) = Gn.block(0, n * k, n, n).trace() / (2 * k);

        if(Gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) {
                for(Index i2 = 0; i2 <= k - i1; i2++) {
                    dks(i1, i2 + (k - i1 - i2) * (m + 1)) /= 1e10;
                }
            }
            Gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = 0; i2 <= m - i1; i2++) {
            //         for(Index i3 = std::max(k - i1 - i2, 0); i3 <= m - i1 - i2; i3++) {
            //             lscf(i1, i2 + i3 * (m + 1)) -= LN1E10;
            //         }
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d3_ijk_mE(const MatrixBase<MatrixXd>& A1,
                            const MatrixBase<MatrixXd>& A2,
                            const MatrixBase<MatrixXd>& A3,
                            const Index m, ArrayXd& lscf, 
                            const double thr_margin, int nthreads);
template ArrayXXl d3_ijk_mE(const MatrixBase<MatrixXl>& A1,
                            const MatrixBase<MatrixXl>& A2,
                            const MatrixBase<MatrixXl>& A3,
                            const Index m, ArrayXl& lscf, 
                            const long double thr_margin, int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (m + 1) * m / 2);
    ArrayXXx Gn = ArrayXXx::Zero(n, (m + 2) * (m + 1) / 2);
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (k + 1) / 2) = Gn.block(0, 0, n, k * (k + 1) / 2);

        Gn.col(0) = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col(0));
        dks(0, k * (m + 1)) = Gn.col(0).sum() / (2 * k);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            Gn.col(id3(0, i2, k)) =
                A2 * (dks(0, i2 - 1 + i3 * (m + 1)) + Go.col(id3(0, i2 - 1, k - 1))) +
                A3 * (dks(0, i2 + (i3 - 1) * (m + 1)) + Go.col(id3(0, i2, k - 1)));
            dks(0, i2 + i3 * (m + 1)) = Gn.col(id3(0, i2, k)).sum() / (2 * k);
        }
        Gn.col(id3(0, k, k)) = A2 * (dks(0, k - 1) + Go.col(id3(0, k - 1, k - 1)));
        dks(0, k) = Gn.col(id3(0, k, k)).sum() / (2 * k);
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            Gn.col(i1) =
                A1 * (dks(i1 - 1, (k - i1) * (m + 1)) + Go.col(i1 - 1)) +
                A3 * (dks(i1, (k - i1 - 1) * (m + 1)) + Go.col(i1));
            dks(i1, (k - i1) * (m + 1)) = Gn.col(i1).sum() / (2 * k);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                Gn.col(id3(i1, i2, k)) =
                    A1 * (dks(i1 - 1, i2 + i3 * (m + 1)) + Go.col(id3(i1 - 1, i2, k - 1))) +
                    A2 * (dks(i1, i2 - 1 + i3 * (m + 1)) + Go.col(id3(i1, i2 - 1, k - 1))) +
                    A3 * (dks(i1, i2 + (i3 - 1) * (m + 1)) + Go.col(id3(i1, i2, k - 1)));
                dks(i1, i2 + i3 * (m + 1)) = Gn.col(id3(i1, i2, k)).sum() / (2 * k);
            }
            Gn.col(id3(i1, k - i1, k)) =
                A1 * (dks(i1 - 1, k - i1) + Go.col(id3(i1 - 1, k - i1, k - 1))) +
                A2 * (dks(i1, k - i1 - 1) + Go.col(id3(i1, k - i1 - 1, k - 1)));
            dks(i1, k - i1) = Gn.col(id3(i1, k - i1, k)).sum() / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        Gn.col(k) = A1 * (dks(k - 1, 0) + Go.col(k - 1));
        dks(k, 0) = Gn.col(k).sum() / (2 * k);

        if(Gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) {
                for(Index i2 = 0; i2 <= k - i1; i2++) {
                    dks(i1, i2 + (k - i1 - i2) * (m + 1)) /= 1e10;
                }
            }
            Gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = 0; i2 <= m - i1; i2++) {
            //         for(Index i3 = std::max(k - i1 - i2, 0); i3 <= m - i1 - i2; i3++) {
            //             lscf(i1, i2 + i3 * (m + 1)) -= LN1E10;
            //         }
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d3_ijk_vE(const ArrayBase<ArrayXd>& A1,
                            const ArrayBase<ArrayXd>& A2,
                            const ArrayBase<ArrayXd>& A3,
                            const Index m, ArrayXd& lscf,
                            const double thr_margin, int nthreads);
template ArrayXXl d3_ijk_vE(const ArrayBase<ArrayXl>& A1,
                            const ArrayBase<ArrayXl>& A2,
                            const ArrayBase<ArrayXl>& A3,
                            const Index m, ArrayXl& lscf,
                            const long double thr_margin, int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Index m, const Eigen::Index p,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf, 
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx Go = MatrixXx::Zero(n, n * (p + 1) * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
            A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        dks(i, 0) = Gn.block(0, i * n, n, n).trace() / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =
            A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        dks(0, k) = Gn.block(0, 0, n, n).trace() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
                A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            dks(i, k) = Gn.block(0, i * n, n, n).trace() / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            MatrixXx::Map(Gn.block(0, j * n * (p + 1), n, n).data(), n, n).noalias() =
                A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
            dks(0, (k - j) + j * (m + 1)) = Gn.block(0, j * n * (p + 1), n, n).trace() / (2 * k);
            for(Index i = 1; i <= p; i++) {
                MatrixXx::Map(Gn.block(0, j * n * (p + 1) + i * n, n, n).data(), n, n).noalias() =
                    A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                    A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                dks(i, (k - j) + j * (m + 1)) = Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        MatrixXx::Map(Gn.block(0, k * n * (p + 1), n, n).data(), n, n).noalias() =
            A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
        dks(0, k * (m + 1)) = Gn.block(0, k * n * (p + 1), n, n).trace() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            MatrixXx::Map(Gn.block(0, k * n * (p + 1) + i * n, n, n).data(), n, n).noalias() =
                A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
            dks(i, k * (m + 1)) = Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr) {
            for(Index j = 0; j <= k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            // for(Index j = 0; j < k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= LN1E10;
            // lscf.block(0, k * (m + 1), p + 1, (m + 1) * (m + 1 - k)) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d3_pjk_mE(const MatrixBase<MatrixXd>& A1,
                            const MatrixBase<MatrixXd>& A2,
                            const MatrixBase<MatrixXd>& A3,
                            const Eigen::Index m, const Eigen::Index p,
                            ArrayXd& lscf, const double thr_margin,
                            int nthreads);
template ArrayXXl d3_pjk_mE(const MatrixBase<MatrixXl>& A1,
                            const MatrixBase<MatrixXl>& A2,
                            const MatrixBase<MatrixXl>& A3,
                            const Eigen::Index m, const Eigen::Index p,
                            ArrayXl& lscf, const long double thr_margin,
                            int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::Index m, const Eigen::Index p,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx Gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        dks(i, 0) = Gn.col(i).sum() / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 500 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        Gn.col(0) = A2 * (dks(0, k - 1) + Go.col(0));
        dks(0, k) = Gn.col(0).sum() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            Gn.col(i) = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                        A2 * (dks(i, k - 1) + Go.col(i));
            dks(i, k) = Gn.col(i).sum() / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            Gn.col(j * (p + 1)) =
                A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1))) +
                A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1)));
            dks(0, (k - j) + j * (m + 1)) = Gn.col(j * (p + 1)).sum() / (2 * k);
            for(Index i = 1; i <= p; i++) {
                Gn.col(j * (p + 1) + i) =
                    A1 * (dks(i - 1, (k - j) + j * (m + 1)) + Gn.col(j * (p + 1) + (i - 1))) +
                    A2 * (dks(i, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1) + i)) +
                    A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1) + i));
                dks(i, (k - j) + j * (m + 1)) = Gn.col(j * (p + 1) + i).sum() / (2 * (k + i));
            }
        }
#ifdef _OPENMP
}
#endif
        Gn.col(k * (p + 1)) = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1)));
        dks(0, k * (m + 1)) = Gn.col(k * (p + 1)).sum() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            Gn.col(k * (p + 1) + i) =
                A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (p + 1) + i - 1)) +
                A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1) + i));
            dks(i, k * (m + 1)) = Gn.col(k * (p + 1) + i).sum() / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr) {
            for(Index j = 0; j <= k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            // for(Index j = 0; j < k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= LN1E10;
            // lscf.block(0, k * (m + 1), p + 1, (m + 1) * (m + 1 - k)) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd d3_pjk_vE(const ArrayBase<ArrayXd>& A1,
                            const ArrayBase<ArrayXd>& A2,
                            const ArrayBase<ArrayXd>& A3,
                            const Index m, const Index p,
                            ArrayXd& lscf, const double thr_margin,
                            int nthreads);
template ArrayXXl d3_pjk_vE(const ArrayBase<ArrayXl>& A1,
                            const ArrayBase<ArrayXl>& A2,
                            const ArrayBase<ArrayXl>& A3,
                            const Index m, const Index p,
                            ArrayXl& lscf, const long double thr_margin,
                            int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * (m + 1) * m / 2);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 2) * (m + 1) / 2);
    MatrixXx go = MatrixXx::Zero(n, (m + 1) * m / 2);
    MatrixXx gn = MatrixXx::Zero(n, (m + 2) * (m + 1) / 2);
    for(Index k = 1; k <= m; k++) {
        if(k % 50 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k * (k + 1) / 2) = Gn.block(0, 0, n, n * k * (k + 1) / 2);
        go.block(0, 0, n, k * (k + 1) / 2) = gn.block(0, 0, n, k * (k + 1) / 2);

        tG.noalias() = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, 0, n, n));
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG - Go.block(0, 0, n, n)
             - (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k * (m + 1)) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            tG.noalias() = A2 * (dks(0, i2 - 1 + i3 * (m + 1)) * In + Go.block(0, n * id3(0, i2 - 1, k - 1), n, n)) +
                           A3 * (dks(0, i2 + (i3 - 1) * (m + 1)) * In + Go.block(0, n * id3(0, i2, k - 1), n, n));
            MatrixXx::Map(gn.col(id3(0, i2, k)).data(), n, 1).noalias() =
                (tG -
                 Go.block(0, n * id3(0, i2 - 1, k - 1), n, n) - Go.block(0, n * id3(0, i2, k - 1), n, n)
                 - ((dks(0, i2 - 1 + i3 * (m + 1)) +
                    dks(0, i2 + (i3 - 1) * (m + 1))) * In)) * mu +
                A2 * go.col(id3(0, i2 - 1, k - 1)) + A3 * go.col(id3(0, i2, k - 1));
            Gn.block(0, n * id3(0, i2, k), n, n) = tG;
            dks(0, i2 + i3 * (m + 1)) = (Gn.block(0, n * id3(0, i2, k), n, n).trace() + gn.col(id3(0, i2, k)).dot(mu)) / (2 * k);
        }
        tG.noalias() = A2 * (dks(0, k - 1) * In + Go.block(0, n * id3(0, k - 1, k - 1), n, n));
        MatrixXx::Map(gn.col(id3(0, k, k)).data(), n, 1).noalias() =
            (tG - Go.block(0, n * id3(0, k - 1, k - 1), n, n)
             - ((dks(0, k - 1)) * In)) * mu + A2 * go.col(id3(0, k - 1, k - 1));
        Gn.block(0, n * id3(0, k, k), n, n) = tG;
        dks(0, k) = (Gn.block(0, n * id3(0, k, k), n, n).trace() + gn.col(id3(0, k, k)).dot(mu)) / (2 * k);
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            tG.noalias() = A1 * (dks(i1 - 1, (k - i1) * (m + 1)) * In + Go.block(0, n * (i1 - 1), n, n)) +
                           A3 * (dks(i1, (k - i1 - 1) * (m + 1)) * In + Go.block(0, n * i1, n, n));
            MatrixXx::Map(gn.col(i1).data(), n, 1).noalias() =
                (tG - Go.block(0, n * (i1 - 1), n, n) - Go.block(0, n * i1, n, n)
                 - ((dks(i1 - 1, (k - i1) * (m + 1)) + dks(i1, (k - i1 - 1) * (m + 1))) * In)) * mu +
                A1 * go.col(i1 - 1) + A3 * go.col(i1);
            Gn.block(0, n * i1, n, n) = tG;
            dks(i1, (k - i1) * (m + 1)) = (Gn.block(0, n * i1, n, n).trace() + gn.col(i1).dot(mu)) / (2 * k);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                tG.noalias() = A1 * (dks(i1 - 1, i2 + i3 * (m + 1)) * In + Go.block(0, n * id3(i1 - 1, i2, k - 1), n, n)) +
                               A2 * (dks(i1, i2 - 1 + i3 * (m + 1)) * In + Go.block(0, n * id3(i1, i2 - 1, k - 1), n, n)) +
                               A3 * (dks(i1, i2 + (i3 - 1) * (m + 1)) * In + Go.block(0, n * id3(i1, i2, k - 1), n, n));
                MatrixXx::Map(gn.col(id3(i1, i2, k)).data(), n, 1).noalias() =
                    (tG - Go.block(0, n * id3(i1 - 1, i2, k - 1), n, n) -
                     Go.block(0, n * id3(i1, i2 - 1, k - 1), n, n) - Go.block(0, n * id3(i1, i2, k - 1), n, n)
                     - ((dks(i1 - 1, i2 + i3 * (m + 1)) + dks(i1, i2 - 1 + i3 * (m + 1)) +
                        dks(i1, i2 + (i3 - 1) * (m + 1))) * In)) * mu +
                    A1 * go.col(id3(i1 - 1, i2, k - 1)) + A2 * go.col(id3(i1, i2 - 1, k - 1)) + A3 * go.col(id3(i1, i2, k - 1));
                Gn.block(0, n * id3(i1, i2, k), n, n) = tG;
                dks(i1, i2 + i3 * (m + 1)) = (Gn.block(0, n * id3(i1, i2, k), n, n).trace() + gn.col(id3(i1, i2, k)).dot(mu)) / (2 * k);
            }
            tG.noalias() = A1 * (dks(i1 - 1, k - i1) * In + Go.block(0, n * id3(i1 - 1, k - i1, k - 1), n, n)) +
                           A2 * (dks(i1, k - i1 - 1) * In + Go.block(0, n * id3(i1, k - i1 - 1, k - 1), n, n));
            MatrixXx::Map(gn.col(id3(i1, k - i1, k)).data(), n, 1).noalias() =
                (tG - Go.block(0, n * id3(i1 - 1, k - i1, k - 1), n, n) -
                 Go.block(0, n * id3(i1, k - i1 - 1, k - 1), n, n)
                 - ((dks(i1 - 1, k - i1) + dks(i1, k - i1 - 1)) * In)) * mu +
                A1 * go.col(id3(i1 - 1, k - i1, k - 1)) + A2 * go.col(id3(i1, k - i1 - 1, k - 1));
            Gn.block(0, n * id3(i1, k - i1, k), n, n) = tG;
            dks(i1, k - i1) = (Gn.block(0, n * id3(i1, k - i1, k), n, n).trace() + gn.col(id3(i1, k - i1, k)).dot(mu)) / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        //
        tG.noalias() = A1 * (dks(k - 1, 0) * In + Go.block(0, n * (k - 1), n, n));
        MatrixXx::Map(gn.col(k).data(), n, 1).noalias() =
            (tG - Go.block(0, n * (k - 1), n, n)
             - ((dks(k - 1, 0)) * In)) * mu + A1 * go.col(k - 1);
        Gn.block(0, n * k, n, n) = tG;
        dks(k, 0) = (Gn.block(0, n * k, n, n).trace() + gn.col(k).dot(mu)) / (2 * k);

        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) {
                for(Index i2 = 0; i2 <= k - i1; i2++) {
                    dks(i1, i2 + (k - i1 - i2) * (m + 1)) /= 1e10;
                }
            }
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = 0; i2 <= m - i1; i2++) {
            //         for(Index i3 = std::max(k - i1 - i2, 0); i3 <= m - i1 - i2; i3++) {
            //             lscf(i1, i2 + i3 * (m + 1)) -= LN1E10;
            //         }
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd  h3_ijk_mE(const MatrixBase<MatrixXd>& A1,
                             const MatrixBase<MatrixXd>& A2,
                             const MatrixBase<MatrixXd>& A3,
                             const VectorXd mu, const Index m,
                             ArrayXd& lscf, const double thr_margin,
                             int nthreads);
template ArrayXXl  h3_ijk_mE(const MatrixBase<MatrixXl>& A1,
                             const MatrixBase<MatrixXl>& A2,
                             const MatrixBase<MatrixXl>& A3,
                             const VectorXl mu, const Index m,
                             ArrayXl& lscf, const long double thr_margin,
                             int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::ArrayBase<Derived>& mu, const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (m + 1) * m / 2);
    ArrayXXx Gn = ArrayXXx::Zero(n, (m + 2) * (m + 1) / 2);
    ArrayXXx go = ArrayXXx::Zero(n, (m + 1) * m / 2);
    ArrayXXx gn = ArrayXXx::Zero(n, (m + 2) * (m + 1) / 2);
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (k + 1) / 2) = Gn.block(0, 0, n, k * (k + 1) / 2);
        go.block(0, 0, n, k * (k + 1) / 2) = gn.block(0, 0, n, k * (k + 1) / 2);

        tG = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col(0));
        gn.col(0) = (tG - Go.col(0)
             - (dks(0, (k - 1) * (m + 1)))) * mu + A3 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k * (m + 1)) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            tG = A2 * (dks(0, i2 - 1 + i3 * (m + 1)) + Go.col(id3(0, i2 - 1, k - 1))) +
                 A3 * (dks(0, i2 + (i3 - 1) * (m + 1)) + Go.col(id3(0, i2, k - 1)));
            gn.col(id3(0, i2, k)) =
                (tG - Go.col(id3(0, i2 - 1, k - 1)) - Go.col(id3(0, i2, k - 1))
                 - ((dks(0, i2 - 1 + i3 * (m + 1)) +
                    dks(0, i2 + (i3 - 1) * (m + 1))))) * mu +
                A2 * go.col(id3(0, i2 - 1, k - 1)) + A3 * go.col(id3(0, i2, k - 1));
            Gn.col(id3(0, i2, k)) = tG;
            dks(0, i2 + i3 * (m + 1)) = (Gn.col(id3(0, i2, k)).sum() + (mu * gn.col(id3(0, i2, k))).sum()) / (2 * k);
        }
        tG = A2 * (dks(0, k - 1) + Go.col(id3(0, k - 1, k - 1)));
        gn.col(id3(0, k, k)) = (tG - Go.col(id3(0, k - 1, k - 1))
             - ((dks(0, k - 1)))) * mu + A2 * go.col(id3(0, k - 1, k - 1));
        Gn.col(id3(0, k, k)) = tG;
        dks(0, k) = (Gn.col(id3(0, k, k)).sum() + (mu * gn.col(id3(0, k, k))).sum()) / (2 * k);
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            tG = A1 * (dks(i1 - 1, (k - i1) * (m + 1)) + Go.col(i1 - 1)) +
                 A3 * (dks(i1, (k - i1 - 1) * (m + 1)) + Go.col(i1));
            gn.col(i1) =
                (tG - Go.col(i1 - 1) - Go.col(i1)
                 - ((dks(i1 - 1, (k - i1) * (m + 1)) + dks(i1, (k - i1 - 1) * (m + 1))))) * mu +
                A1 * go.col(i1 - 1) + A3 * go.col(i1);
            Gn.col(i1) = tG;
            dks(i1, (k - i1) * (m + 1)) = (Gn.col(i1).sum() + (mu * gn.col(i1)).sum()) / (2 * k);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                tG = A1 * (dks(i1 - 1, i2 + i3 * (m + 1)) + Go.col(id3(i1 - 1, i2, k - 1))) +
                     A2 * (dks(i1, i2 - 1 + i3 * (m + 1)) + Go.col(id3(i1, i2 - 1, k - 1))) +
                     A3 * (dks(i1, i2 + (i3 - 1) * (m + 1)) + Go.col(id3(i1, i2, k - 1)));
                gn.col(id3(i1, i2, k)) =
                    (tG - Go.col(id3(i1 - 1, i2, k - 1)) -
                     Go.col(id3(i1, i2 - 1, k - 1)) - Go.col(id3(i1, i2, k - 1))
                     - ((dks(i1 - 1, i2 + i3 * (m + 1)) + dks(i1, i2 - 1 + i3 * (m + 1)) +
                        dks(i1, i2 + (i3 - 1) * (m + 1))))) * mu +
                    A1 * go.col(id3(i1 - 1, i2, k - 1)) + A2 * go.col(id3(i1, i2 - 1, k - 1)) + A3 * go.col(id3(i1, i2, k - 1));
                Gn.col(id3(i1, i2, k)) = tG;
                dks(i1, i2 + i3 * (m + 1)) = (Gn.col(id3(i1, i2, k)).sum() + (mu * gn.col(id3(i1, i2, k))).sum()) / (2 * k);
            }
            tG = A1 * (dks(i1 - 1, k - i1) + Go.col(id3(i1 - 1, k - i1, k - 1))) +
                 A2 * (dks(i1, k - i1 - 1) + Go.col(id3(i1, k - i1 - 1, k - 1)));
            gn.col(id3(i1, k - i1, k)) =
                (tG - Go.col(id3(i1 - 1, k - i1, k - 1)) -
                 Go.col(id3(i1, k - i1 - 1, k - 1))
                 - ((dks(i1 - 1, k - i1) + dks(i1, k - i1 - 1)))) * mu +
                A1 * go.col(id3(i1 - 1, k - i1, k - 1)) + A2 * go.col(id3(i1, k - i1 - 1, k - 1));
            Gn.col(id3(i1, k - i1, k)) = tG;
            dks(i1, k - i1) = (Gn.col(id3(i1, k - i1, k)).sum() + (mu * gn.col(id3(i1, k - i1, k))).sum()) / (2 * k);
        }
#ifdef _OPENMP
}
#endif
        tG = A1 * (dks(k - 1, 0) + Go.col(k - 1));
        gn.col(k) = (tG - Go.col(k - 1)
             - ((dks(k - 1, 0)))) * mu + A1 * go.col(k - 1);
        Gn.col(k) = tG;
        dks(k, 0) = (Gn.col(k).sum() + (mu * gn.col(k)).sum()) / (2 * k);

        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index i1 = 0; i1 <= k; i1++) {
                for(Index i2 = 0; i2 <= k - i1; i2++) {
                    dks(i1, i2 + (k - i1 - i2) * (m + 1)) /= 1e10;
                }
            }
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index i1 = 0; i1 <= m; i1++) {
            //     for(Index i2 = 0; i2 <= m - i1; i2++) {
            //         for(Index i3 = std::max(k - i1 - i2, 0); i3 <= m - i1 - i2; i3++) {
            //             lscf(i1, i2 + i3 * (m + 1)) -= LN1E10;
            //         }
            //     }
            // }
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd h3_ijk_vE(const ArrayBase<ArrayXd>& A1,
                            const ArrayBase<ArrayXd>& A2,
                            const ArrayBase<ArrayXd>& A3,
                            const ArrayBase<ArrayXd>& mu, const Index m,
                            ArrayXd& lscf, const double thr_margin,
                            int nthreads);
template ArrayXXl h3_ijk_vE(const ArrayBase<ArrayXl>& A1,
                            const ArrayBase<ArrayXl>& A2,
                            const ArrayBase<ArrayXl>& A3,
                            const ArrayBase<ArrayXl>& mu, const Index m,
                            ArrayXl& lscf, const long double thr_margin,
                            int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * (p + 1) * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (p + 1) * (m + 1));
    MatrixXx go = MatrixXx::Zero(n, (p + 1) * m);
    MatrixXx gn = MatrixXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        MatrixXx::Map(gn.col(i).data(), n, 1).noalias() = Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG.noalias() = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG - Go.block(0, 0, n, n) - (dks(0, k - 1) * In)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                           A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            MatrixXx::Map(gn.col(i).data(), n, 1).noalias() =
                (tG - Go.block(0, i * n, n, n) - (dks(i, (k - 1)) * In)) * mu +
                A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            tG.noalias() = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) +
                           A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
            MatrixXx::Map(gn.col(j * (p + 1)).data(), n, 1).noalias() =
                (tG - Go.block(0, j * n * (p + 1), n, n) - Go.block(0, (j - 1) * n * (p + 1), n, n) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.block(0, j * n * (p + 1), n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
            for(Index i = 1; i <= p; i++) {
                tG.noalias() = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                               A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                               A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                MatrixXx::Map(gn.col(j * (p + 1) + i).data(), n, 1).noalias() =
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
        tG.noalias() = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
        MatrixXx::Map(gn.col(k * (p + 1)).data(), n, 1).noalias() =
            (tG - Go.block(0, (k - 1) * n * (p + 1), n, n) -
             (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * go.col((k - 1) * (p + 1));
        Gn.block(0, k * n * (p + 1), n, n) = tG;
        dks(0, k * (m + 1)) = (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                           A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
            MatrixXx::Map(gn.col(k * (p + 1) + i).data(), n, 1).noalias() =
                (tG - Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n)
                 - (dks(i, (k - 1) * (m + 1)) * In)) * mu +
                A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.block(0, k * n * (p + 1) + i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() + gn.col(k * (p + 1) + i).dot(mu)) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index j = 0; j <= k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index j = 0; j < k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= LN1E10;
            // lscf.block(0, k * (m + 1), p + 1, (m + 1) * (m + 1 - k)) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd htil3_pjk_mE(const MatrixBase<MatrixXd>& A1,
                               const MatrixBase<MatrixXd>& A2,
                               const MatrixBase<MatrixXd>& A3,
                               const VectorXd mu, const Index m, const Index p,
                               ArrayXd& lscf, const double thr_margin,
                               int nthreads);
template ArrayXXl htil3_pjk_mE(const MatrixBase<MatrixXl>& A1,
                               const MatrixBase<MatrixXl>& A2,
                               const MatrixBase<MatrixXl>& A3,
                               const VectorXl mu, const Index m, const Index p,
                               ArrayXl& lscf, const long double thr_margin,
                               int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx Gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    ArrayXXx go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 500 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks(0, k - 1) + Go.col(0));
        gn.col(0) = (tG - Go.col(0) - (dks(0, k - 1))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                 A2 * (dks(i, k - 1) + Go.col(i));
            gn.col(i) = (tG - Go.col(i)
                         - (dks(i, (k - 1)))) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.col(i) = tG;
            dks(i, k) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1))) +
                 A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1)));
            gn.col(j * (p + 1)) =
                (tG - Go.col(j * (p + 1)) - Go.col((j - 1) * (p + 1)) - ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.col(j * (p + 1)) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
            for(Index i = 1; i <= p; i++) {
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
#ifdef _OPENMP
}
#endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1)));
        gn.col(k * (p + 1)) =
            (tG - Go.col((k - 1) * (p + 1)) - (dks(0, (k - 1) * (m + 1)))) * mu +
            A3 * go.col((k - 1) * (p + 1));
        Gn.col(k * (p + 1)) = tG;
        dks(0, k * (m + 1)) = (Gn.col(k * (p + 1)).sum() + (mu * gn.col(k * (p + 1))).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k * (m + 1)) + Gn.col(k * (p + 1) + i - 1)) +
                 A3 * (dks(i, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1) + i));
            gn.col(k * (p + 1) + i) = (tG - Go.col((k - 1) * (p + 1) + i)
                         - (dks(i, (k - 1) * (m + 1)))) * mu +
                        A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.col(k * (p + 1) + i) = tG;
            dks(i, k * (m + 1)) = (Gn.col(k * (p + 1) + i).sum() + (mu * gn.col(k * (p + 1) + i)).sum()) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index j = 0; j <= k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index j = 0; j < k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= LN1E10;
            // lscf.block(0, k * (m + 1), p + 1, (m + 1) * (m + 1 - k)) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd htil3_pjk_vE(const ArrayBase<ArrayXd>& A1,
                               const ArrayBase<ArrayXd>& A2,
                               const ArrayBase<ArrayXd>& A3,
                               const ArrayBase<ArrayXd>& mu,
                               const Index m, const Index p,
                               ArrayXd& lscf, const double thr_margin,
                               int nthreads);
template ArrayXXl htil3_pjk_vE(const ArrayBase<ArrayXl>& A1,
                               const ArrayBase<ArrayXl>& A2,
                               const ArrayBase<ArrayXl>& A3,
                               const ArrayBase<ArrayXl>& mu,
                               const Index m, const Index p,
                               ArrayXl& lscf, const long double thr_margin,
                               int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * (p + 1) * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (p + 1) * (m + 1));
    MatrixXx go = MatrixXx::Zero(n, (p + 1) * m);
    MatrixXx gn = MatrixXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
            A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        MatrixXx::Map(gn.col(i).data(), n, 1).noalias() =
            Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG.noalias() = A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG + Go.block(0, 0, n, n) + (dks(0, k - 1) * In)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                           A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
            MatrixXx::Map(gn.col(i).data(), n, 1).noalias() =
                (tG + Go.block(0, i * n, n, n) + (dks(i, (k - 1)) * In)) * mu +
                A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.block(0, i * n, n, n) = tG;
            dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            tG.noalias() = A2 * (dks(0, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) +
                           A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
            MatrixXx::Map(gn.col(j * (p + 1)).data(), n, 1).noalias() =
                (tG + Go.block(0, j * n * (p + 1), n, n) + Go.block(0, (j - 1) * n * (p + 1), n, n) +
                 ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))) * In)) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.block(0, j * n * (p + 1), n, n) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
            for(Index i = 1; i <= p; i++) {
                tG.noalias() = A1 * (dks(i - 1, (k - j) + j * (m + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                               A2 * (dks(i, (k - j - 1) + j * (m + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                               A3 * (dks(i, (k - j) + (j - 1) * (m + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                MatrixXx::Map(gn.col(j * (p + 1) + i).data(), n, 1).noalias() =
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
        tG.noalias() = A3 * (dks(0, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
        MatrixXx::Map(gn.col(k * (p + 1)).data(), n, 1).noalias() =
            (tG + Go.block(0, (k - 1) * n * (p + 1), n, n) +
             (dks(0, (k - 1) * (m + 1)) * In)) * mu + A3 * go.col((k - 1) * (p + 1));
        Gn.block(0, k * n * (p + 1), n, n) = tG;
        dks(0, k * (m + 1)) = (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * (dks(i - 1, k * (m + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                           A3 * (dks(i, (k - 1) * (m + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
            MatrixXx::Map(gn.col(k * (p + 1) + i).data(), n, 1).noalias() =
                (tG + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n)
                 + (dks(i, (k - 1) * (m + 1)) * In)) * mu +
                A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.block(0, k * n * (p + 1) + i * n, n, n) = tG;
            dks(i, k * (m + 1)) = (Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() + gn.col(k * (p + 1) + i).dot(mu)) / (2 * (k + i));
        }
        if(Gn.maxCoeff() > thr || gn.maxCoeff() > thr) {
            for(Index j = 0; j <= k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index j = 0; j < k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= LN1E10;
            // lscf.block(0, k * (m + 1), p + 1, (m + 1) * (m + 1 - k)) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd hhat3_pjk_mE(const MatrixBase<MatrixXd>& A1,
                               const MatrixBase<MatrixXd>& A2,
                               const MatrixBase<MatrixXd>& A3,
                               const VectorXd mu, const Index m, const Index p,
                               ArrayXd& lscf, const double thr_margin,
                               int nthreads);
template ArrayXXl hhat3_pjk_mE(const MatrixBase<MatrixXl>& A1,
                               const MatrixBase<MatrixXl>& A2,
                               const MatrixBase<MatrixXl>& A3,
                               const VectorXl mu, const Index m, const Index p,
                               ArrayXl& lscf, const long double thr_margin,
                               int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads == 0) nthreads = omp_get_num_procs() / 2;
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 1));
    dks(0, 0) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx Gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    ArrayXXx go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 500 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks(0, k - 1) + Go.col(0));
        gn.col(0) = (tG + Go.col(0) + (dks(0, k - 1))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG = A1 * (dks(i - 1, k) + Gn.col(i - 1)) +
                 A2 * (dks(i, k - 1) + Go.col(i));
            gn.col(i) = (tG + Go.col(i)
                         + (dks(i, (k - 1)))) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.col(i) = tG;
            dks(i, k) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
        }
#ifdef _OPENMP
#pragma omp parallel private(tG)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            tG = A2 * (dks(0, (k - j - 1) + j * (m + 1)) + Go.col(j * (p + 1))) + A3 * (dks(0, (k - j) + (j - 1) * (m + 1)) + Go.col((j - 1) * (p + 1)));
            gn.col(j * (p + 1)) =
                (tG + Go.col(j * (p + 1)) + Go.col((j - 1) * (p + 1)) +
                 ((dks(0, (k - j - 1) + j * (m + 1)) + dks(0, (k - j) + (j - 1) * (m + 1))))) * mu +
                A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
            Gn.col(j * (p + 1)) = tG;
            dks(0, (k - j) + j * (m + 1)) = (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
            for(Index i = 1; i <= p; i++) {
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
#ifdef _OPENMP
}
#endif
        tG = A3 * (dks(0, (k - 1) * (m + 1)) + Go.col((k - 1) * (p + 1)));
        gn.col(k * (p + 1)) =
            (tG + Go.col((k - 1) * (p + 1)) + (dks(0, (k - 1) * (m + 1)))) * mu +
            A3 * go.col((k - 1) * (p + 1));
        Gn.col(k * (p + 1)) = tG;
        dks(0, k * (m + 1)) = (Gn.col(k * (p + 1)).sum() + (mu * gn.col(k * (p + 1))).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
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
            for(Index j = 0; j <= k; j++) dks.col((k - j) + j * (m + 1)) /= 1e10;
            Gn /= 1e10;
            gn /= 1e10;
            // for(Index j = 0; j < k; j++) lscf.block(0, (k - j) + j * (m + 1), p + 1, m + 1 - (k - j)) -= LN1E10;
            // lscf.block(0, k * (m + 1), p + 1, (m + 1) * (m + 1 - k)) -= LN1E10;
            lscf.tail(m + 1 - k) -= LN1E10;
        }
    }
    return dks;
}
template ArrayXXd hhat3_pjk_vE(const ArrayBase<ArrayXd>& A1,
                               const ArrayBase<ArrayXd>& A2,
                               const ArrayBase<ArrayXd>& A3,
                               const ArrayBase<ArrayXd>& mu,
                               const Index m, const Index p,
                               ArrayXd& lscf, const double thr_margin,
                               int nthreads);
template ArrayXXl hhat3_pjk_vE(const ArrayBase<ArrayXl>& A1,
                               const ArrayBase<ArrayXl>& A2,
                               const ArrayBase<ArrayXl>& A3,
                               const ArrayBase<ArrayXl>& mu,
                               const Index m, const Index p,
                               ArrayXl& lscf, const long double thr_margin,
                               int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil3_pqr_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index p, const Eigen::Index q, const Eigen::Index r) { //, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf, const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const MatrixXx In = MatrixXx::Identity(n, n);
    const Index m = q + r;
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (q + 1) * (r + 1));
    dks(0, 0) = 1;
    // Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    const MatrixXx zeromat_n_np = MatrixXx::Zero(n, n * (p + 1));
    const MatrixXx zeromat_n_p = MatrixXx::Zero(n, p + 1);
    MatrixXx Go = MatrixXx::Zero(n, n * (p + 1) * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (p + 1) * (m + 1));
    MatrixXx go = MatrixXx::Zero(n, (p + 1) * m);
    MatrixXx gn = MatrixXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() = A1 * (dks(i - 1, 0) * In + Gn.block(0, (i - 1) * n, n, n));
        MatrixXx::Map(gn.col(i).data(), n, 1).noalias() = Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 500 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        if(k <= q) {
            MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =A2 * (dks(0, k - 1) * In + Go.block(0, 0, n, n));
            MatrixXx::Map(gn.col(0).data(), n, 1).noalias() = Gn.block(0, 0, n, n) * mu + A2 * go.col(0);
            dks(0, k) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
            for(Index i = 1; i <= p; i++) {
                MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
                    A1 * (dks(i - 1, k) * In + Gn.block(0, (i - 1) * n, n, n)) +
                    A2 * (dks(i, k - 1) * In + Go.block(0, i * n, n, n));
                MatrixXx::Map(gn.col(i).data(), n, 1).noalias() = Gn.block(0, i * n, n, n) * mu +
                            A1 * gn.col(i - 1) + A2 * go.col(i);
                dks(i, k) = (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
            }
        } else {
            Gn.block(0, 0, n, n * (p + 1)) = zeromat_n_p;
            gn.block(0, 0, n, p + 1) = zeromat_n_p;
        }
        for(Index j = 1; j < k; j++) {
            if(((k - j) <= q) && (j <= r)) {
                MatrixXx::Map(Gn.block(0, j * n * (p + 1), n, n).data(), n, n).noalias() =
                    A2 * (dks(0, (k - j - 1) + j * (q + 1)) * In + Go.block(0, j * n * (p + 1), n, n)) +
                    A3 * (dks(0, (k - j) + (j - 1) * (q + 1)) * In + Go.block(0, (j - 1) * n * (p + 1), n, n));
                MatrixXx::Map(gn.col(j * (p + 1)).data(), n, 1).noalias() =
                    Gn.block(0, j * n * (p + 1), n, n) * mu +
                    A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
                dks(0, (k - j) + j * (q + 1)) = (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
                for(Index i = 1; i <= p; i++) {
                    MatrixXx::Map(Gn.block(0, j * n * (p + 1) + i * n, n, n).data(), n, n).noalias() =
                        A1 * (dks(i - 1, (k - j) + j * (q + 1)) * In + Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n)) +
                        A2 * (dks(i, (k - j - 1) + j * (q + 1)) * In + Go.block(0, j * n * (p + 1) + i * n, n, n)) +
                        A3 * (dks(i, (k - j) + (j - 1) * (q + 1)) * In + Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n));
                    MatrixXx::Map(gn.col(j * (p + 1) + i).data(), n, 1).noalias() =
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
            MatrixXx::Map(Gn.block(0, k * n * (p + 1), n, n).data(), n, n).noalias() =
                A3 * (dks(0, (k - 1) * (q + 1)) * In + Go.block(0, (k - 1) * n * (p + 1), n, n));
            MatrixXx::Map(gn.col(k * (p + 1)).data(), n, 1).noalias() =
                Gn.block(0, k * n * (p + 1), n, n) * mu + A3 * go.col((k - 1) * (p + 1));
            dks(0, k * (q + 1)) = (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
            for(Index i = 1; i <= p; i++) {
                MatrixXx::Map(Gn.block(0, k * n * (p + 1) + i * n, n, n).data(), n, n).noalias() =
                    A1 * (dks(i - 1, k * (q + 1)) * In + Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n)) +
                    A3 * (dks(i, (k - 1) * (q + 1)) * In + Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n));
                MatrixXx::Map(gn.col(k * (p + 1) + i).data(), n, 1).noalias() =
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
        //     lscf -= LN1E10;
        // }
    }
    return dks;
}
template ArrayXXd dtil3_pqr_mE(const MatrixBase<MatrixXd>& A1,
                               const MatrixBase<MatrixXd>& A2,
                               const MatrixBase<MatrixXd>& A3,
                               const VectorXd mu,
                               const Index p, const Index q, const Index r);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil3_pqr_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index p, const Eigen::Index q, const Eigen::Index r) { //, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf, const typename Derived::Scalar thr_margin) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    const Index m = q + r;
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (q + 1) * (r + 1));
    dks(0, 0) = 1;
    // Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    const ArrayXXx zeromat_n_p = ArrayXXx::Zero(n, p + 1);
    ArrayXXx Go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx Gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    ArrayXXx go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    for(Index i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks(i - 1, 0) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks(i, 0) = (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    for(Index k = 1; k <= m; k++) {
        if(k % 10000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        if(k <= q) {
            Gn.col(0) = A2 * (dks(0, k - 1) + Go.col(0));
            gn.col(0) = Gn.col(0) * mu + A2 * go.col(0);
            dks(0, k) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
            for(Index i = 1; i <= p; i++) {
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
        for(Index j = 1; j < k; j++) {
            if(((k - j) <= q) && (j <= r)) {
                Gn.col(j * (p + 1)) =
                    A2 * (dks(0, (k - j - 1) + j * (q + 1)) + Go.col(j * (p + 1))) +
                    A3 * (dks(0, (k - j) + (j - 1) * (q + 1)) + Go.col((j - 1) * (p + 1)));
                gn.col(j * (p + 1)) =
                    Gn.col(j * (p + 1)) * mu +
                    A2 * go.col(j * (p + 1)) + A3 * go.col((j - 1) * (p + 1));
                dks(0, (k - j) + j * (q + 1)) = (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
                for(Index i = 1; i <= p; i++) {
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
            for(Index i = 1; i <= p; i++) {
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
        //     lscf -= LN1E10;
        // }
    }
    return dks;
}
template ArrayXXd dtil3_pqr_vE(const ArrayBase<ArrayXd>& A1,
                               const ArrayBase<ArrayXd>& A2,
                               const ArrayBase<ArrayXd>& A3,
                               const ArrayBase<ArrayXd>& mu,
                               const Index p, const Index q, const Index r);
