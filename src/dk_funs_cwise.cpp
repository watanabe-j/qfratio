#include "config.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <limits>
#include "dk_funs_cwise.h"

#ifdef _OPENMP
    #include <omp.h>
    // [[Rcpp::plugins(openmp)]]
#endif

#define LN1E10 M_LN10 * 10

using std::min;
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

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatXd;


template <typename Derived>
inline void update_scale_2D(Eigen::ArrayBase<Derived>& lscf,
                            const Eigen::Index i0, const Eigen::Index j0,
                            const Eigen::Index M) {
    lscf.ULTat(i0, j0, M) -= LN1E10;
    typename Derived::Scalar lscf0 = lscf.ULTat(i0, j0, M);
    Index ie = M - j0;
    Index je = M - i0;
    for(Index ii = i0 + 1; ii < M - j0; ii++) {
        if(lscf.ULTat(ii, j0, M) <= lscf0) {
            ie = ii;
            break;
        }
    }
    for(Index jj = j0 + 1; jj < M - i0; jj++) {
        if(lscf.ULTat(i0, jj, M) <= lscf0) {
            je = jj;
            break;
        }
    }
    lscf.ULTcol(j0, M).segment(i0, ie - i0) = lscf0;
    for(Index jj = j0 + 1; jj < je; jj++) {
        lscf.ULTcol(jj, M).segment(i0, min(ie, M - jj) - i0) = lscf0;
    }
}

template <typename Derived>
inline void update_scale_3D(Eigen::ArrayBase<Derived>& lscf,
                            const Eigen::Index i0, const Eigen::Index j0,
                            const Eigen::Index k0, const Eigen::Index M) {
    lscf.ULCat(i0, j0, k0, M) -= LN1E10;
    typename Derived::Scalar lscf0 = lscf.ULCat(i0, j0, k0, M);
    Index ie = M - j0 - k0;
    Index je = M - i0 - k0;
    Index ke = M - i0 - j0;
    for(Index ii = i0 + 1; ii < M - j0 - k0; ii++) {
        if(lscf.ULCat(ii, j0, k0, M) <= lscf0) {
            ie = ii;
            break;
        }
    }
    for(Index jj = j0 + 1; jj < M - i0 - k0; jj++) {
        if(lscf.ULCat(i0, jj, k0, M) <= lscf0) {
            je = jj;
            break;
        }
    }
    for(Index kk = k0 + 1; kk < M - i0 - j0; kk++) {
        if(lscf.ULCat(i0, j0, kk, M) <= lscf0) {
            ke = kk;
            break;
        }
    }
    lscf.ULCcol(j0, k0).segment(i0, ie - i0) = lscf0;
    for(Index jj = j0 + 1; jj < je; jj++) {
        lscf.ULCcol(jj, k0, M).segment(i0, min(ie, M - jj - k0) - i0) = lscf0;
    }
    for(Index kk = k0 + 1; kk < ke; kk++) {
        lscf.ULCcol(j0, kk, M).segment(i0, min(ie, M - j0 - kk) - i0) = lscf0;
        for(Index jj = j0 + 1; jj < je; jj++) {
            lscf.ULCcol(jj, kk, M).segment(i0, min(ie, M - jj - kk) - i0) = lscf0;
        }
    }
}

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


template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_d2_ij_mE(Eigen::Index i1, Eigen::Index k,
                              const Eigen::Index m, const Eigen::Index n, 
                              typename DerivedA::Scalar &thr,
                              Eigen::ArrayBase<DerivedA> &dks,
                              Eigen::ArrayBase<DerivedB> &lscf,
                              Eigen::MatrixBase<DerivedC> &Gn) {
    if(Gn.block(0, n * i1, n, n).maxCoeff() > thr) {
        dks.ULTat(i1, k - i1, m + 1) /= 1e10;
        Gn.block(0, n * i1, n, n) /= 1e10;
        update_scale_2D(lscf, i1, k - i1, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d2_ij_mEc(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) / 2);
    dks.ULTat(0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx Go = MatrixXx::Zero(n, n * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 1));
    Gn.block(0, 0, n, n).diagonal().array() += dks.ULTat(0, 0, m + 1);
    Scalar s1, s2;
    for(Index k = 1; k <= m; k++) {
        if(k % 200 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k) = Gn.block(0, 0, n, n * k);

        MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =
            A2 * Go.block(0, 0, n, n);
        dks.ULTat(0, k, m + 1) = Gn.block(0, 0, n, n).trace() / (2 * k);
        Gn.block(0, 0, n, n).diagonal().array() += dks.ULTat(0, k, m + 1);
        scale_in_d2_ij_mE(0, k, m, n, thr, dks, lscf, Gn);
#ifdef _OPENMP
#pragma omp parallel private(s1, s2)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULTat(i1, k - i1 - 1, m + 1) - lscf.ULTat(i1 - 1, k - i1, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULTat(i1 - 1, k - i1, m + 1) - lscf.ULTat(i1, k - i1 - 1, m + 1)));
            MatrixXx::Map(Gn.block(0, n * i1, n, n).data(), n, n).noalias() =
                s1 * A1 * Go.block(0, n * (i1 - 1), n, n) +
                s2 * A2 * Go.block(0, n * i1, n, n);
            dks.ULTat(i1, k - i1, m + 1) = Gn.block(0, n * i1, n, n).trace() / (2 * k);
            Gn.block(0, n * i1, n, n).diagonal().array() += dks.ULTat(i1, k - i1, m + 1);
            scale_in_d2_ij_mE(i1, k, m, n, thr, dks, lscf, Gn);
        }
#ifdef _OPENMP
}
#endif
        MatrixXx::Map(Gn.block(0, n * k, n, n).data(), n, n).noalias() =
            A1 * Go.block(0, n * (k - 1), n, n);
        dks.ULTat(k, 0, m + 1) = Gn.block(0, n * k, n, n).trace() / (2 * k);
        Gn.block(0, n * k, n, n).diagonal().array() += dks.ULTat(k, 0, m + 1);
        scale_in_d2_ij_mE(k, k, m, n, thr, dks, lscf, Gn);
    }
    return dks;
}
template ArrayXd d2_ij_mEc(const MatrixBase<MatrixXd>& A1,
                           const DiagMatXd& A2,
                           const Index m, ArrayXd& lscf,
                           const double thr_margin, int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_d2_ij_vE(Eigen::Index i1, Eigen::Index k,
                              const Eigen::Index &m, const Eigen::Index &n, 
                              typename DerivedA::Scalar &thr,
                              Eigen::ArrayBase<DerivedA> &dks,
                              Eigen::ArrayBase<DerivedB> &lscf,
                              Eigen::ArrayBase<DerivedC> &Gn) {
    if(Gn.col(i1).maxCoeff() > thr) {
        dks.ULTat(i1, k - i1, m + 1) /= 1e10;
        Gn.col(i1) /= 1e10;
        update_scale_2D(lscf, i1, k - i1, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d2_ij_vEc(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) / 2);
    dks.ULTat(0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx Go = ArrayXXx::Zero(n, m);
    ArrayXXx Gn = ArrayXXx::Zero(n, m + 1);
    Scalar s1, s2;
    for(Index k = 1; k <= m; k++) {
        if(k % 2000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k) = Gn.block(0, 0, n, k);

        Gn.col(0) = A2 * (dks.ULTat(0, k - 1, m + 1) + Go.col(0));
        dks.ULTat(0, k, m + 1) = Gn.col(0).sum() / (2 * k);
        scale_in_d2_ij_vE(0, k, m, n, thr, dks, lscf, Gn);
#ifdef _OPENMP
#pragma omp parallel private(s1, s2)
        {
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULTat(i1, k - i1 - 1, m + 1) - lscf.ULTat(i1 - 1, k - i1, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULTat(i1 - 1, k - i1, m + 1) - lscf.ULTat(i1, k - i1 - 1, m + 1)));
            Gn.col(i1) = s1 * A1 * (dks.ULTat(i1 - 1, k - i1, m + 1) + Go.col(i1 - 1)) +
                         s2 * A2 * (dks.ULTat(i1, k - i1 - 1, m + 1) + Go.col(i1));
            dks.ULTat(i1, k - i1, m + 1) = Gn.col(i1).sum() / (2 * k);
            scale_in_d2_ij_vE(i1, k, m, n, thr, dks, lscf, Gn);
        }
#ifdef _OPENMP
}
#endif
        Gn.col(k) = A1 * (dks.ULTat(k - 1, 0, m + 1) + Go.col(k - 1));
        dks.ULTat(k, 0, m + 1) = Gn.col(k).sum() / (2 * k);
        scale_in_d2_ij_vE(k, k, m, n, thr, dks, lscf, Gn);
    }
    return dks;
}
template ArrayXd d2_ij_vEc(const ArrayBase<ArrayXd>& A1,
                           const ArrayBase<ArrayXd>& A2,
                           const Index m, ArrayXd& lscf,
                           const double thr_margin, int nthreads);


template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
inline void scale_in_h2_ij_mE(Eigen::Index i1, Eigen::Index k,
                              const Eigen::Index &m, const Eigen::Index &n, 
                              typename DerivedA::Scalar &thr,
                              Eigen::ArrayBase<DerivedA> &dks,
                              Eigen::ArrayBase<DerivedB> &lscf,
                              Eigen::MatrixBase<DerivedC> &Gn,
                              Eigen::MatrixBase<DerivedD> &gn) {
    if(Gn.block(0, n * i1, n, n).maxCoeff() > thr || gn.col(i1).maxCoeff() > thr) {
        dks.ULTat(i1, k - i1, m + 1) /= 1e10;
        Gn.block(0, n * i1, n, n) /= 1e10;
        gn.col(i1) /= 1e10;
        update_scale_2D(lscf, i1, k - i1, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h2_ij_mEc(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
         const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) / 2);
    dks.ULTat(0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 1));
    MatrixXx go = MatrixXx::Zero(n, m);
    MatrixXx gn = MatrixXx::Zero(n, m + 1);
    Gn.block(0, 0, n, n).diagonal().array() += dks.ULTat(0, 0, m + 1);
    Scalar s1, s2;
    for(Index k = 1; k <= m; k++) {
        if(k % 200 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k) = Gn.block(0, 0, n, n * k);
        go.block(0, 0, n, k) = gn.block(0, 0, n, k);

        tG.noalias() = A2 * Go.block(0, 0, n, n);
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG - Go.block(0, 0, n, n)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks.ULTat(0, k, m + 1) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        Gn.block(0, 0, n, n).diagonal().array() += dks.ULTat(0, k, m + 1);
        scale_in_h2_ij_mE(0, k, m, n, thr, dks, lscf, Gn, gn);
#ifdef _OPENMP
#pragma omp parallel private(tG, s1, s2)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULTat(i1, k - i1 - 1, m + 1) - lscf.ULTat(i1 - 1, k - i1, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULTat(i1 - 1, k - i1, m + 1) - lscf.ULTat(i1, k - i1 - 1, m + 1)));
            tG.noalias() = s1 * A1 * Go.block(0, n * (i1 - 1), n, n) +
                           s2 * A2 * Go.block(0, n * i1, n, n);
            MatrixXx::Map(gn.col(i1).data(), n, 1).noalias() =
                (tG - s1 * Go.block(0, n * (i1 - 1), n, n) -
                 s2 * Go.block(0, n * i1, n, n)) * mu +
                s1 * A1 * go.col(i1 - 1) + s2 * A2 * go.col(i1);
            Gn.block(0, n * i1, n, n) = tG;
            dks.ULTat(i1, k - i1, m + 1) = (Gn.block(0, n * i1, n, n).trace() + gn.col(i1).dot(mu)) / (2 * k);
            Gn.block(0, n * i1, n, n).diagonal().array() += dks.ULTat(i1, k - i1, m + 1);
            scale_in_h2_ij_mE(i1, k, m, n, thr, dks, lscf, Gn, gn);
        }
#ifdef _OPENMP
}
#endif
        tG.noalias() = A1 * Go.block(0, n * (k - 1), n, n);
        MatrixXx::Map(gn.col(k).data(), n, 1).noalias() =
            (tG - Go.block(0, n * (k - 1), n, n)) * mu + A1 * go.col(k - 1);
        Gn.block(0, n * k, n, n) = tG;
        dks.ULTat(k, 0, m + 1) = (Gn.block(0, n * k, n, n).trace() + gn.col(k).dot(mu)) / (2 * k);
        Gn.block(0, n * k, n, n).diagonal().array() += dks.ULTat(k, 0, m + 1);
        scale_in_h2_ij_mE(k, k, m, n, thr, dks, lscf, Gn, gn);
    }
    return dks;
}
template ArrayXd h2_ij_mEc(const MatrixBase<MatrixXd>& A1,
                           const DiagMatXd& A2,
                           const VectorXd& mu, const Index m,
                           ArrayXd& lscf, const double thr_margin,
                           int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_h2_ij_vE(Eigen::Index i1, Eigen::Index k,
                              const Eigen::Index &m, const Eigen::Index &n, 
                              typename DerivedA::Scalar &thr,
                              Eigen::ArrayBase<DerivedA> &dks,
                              Eigen::ArrayBase<DerivedB> &lscf,
                              Eigen::ArrayBase<DerivedC> &Gn,
                              Eigen::ArrayBase<DerivedC> &gn) {
    if(Gn.col(i1).maxCoeff() > thr || gn.col(i1).maxCoeff() > thr) {
        dks.ULTat(i1, k - i1, m + 1) /= 1e10;
        Gn.col(i1) /= 1e10;
        gn.col(i1) /= 1e10;
        update_scale_2D(lscf, i1, k - i1, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h2_ij_vEc(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::ArrayBase<Derived>& mu,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) / 2);
    dks.ULTat(0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, m);
    ArrayXXx Gn = ArrayXXx::Zero(n, m + 1);
    ArrayXXx go = ArrayXXx::Zero(n, m);
    ArrayXXx gn = ArrayXXx::Zero(n, m + 1);
    Scalar s1, s2;
    for(Index k = 1; k <= m; k++) {
        if(k % 2000 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k) = Gn.block(0, 0, n, k);
        go.block(0, 0, n, k) = gn.block(0, 0, n, k);

        tG = A2 * (dks.ULTat(0, k - 1, m + 1) + Go.col(0));
        gn.col(0) = (tG - Go.col(0)
             - ((dks.ULTat(0, k - 1, m + 1)))) * mu + A2 * go.col(0);
        Gn.col(0) = tG;
        dks.ULTat(0, k, m + 1) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        scale_in_h2_ij_vE(0, k, m, n, thr, dks, lscf, Gn, gn);
#ifdef _OPENMP
#pragma omp parallel private(tG, s1, s2)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULTat(i1, k - i1 - 1, m + 1) - lscf.ULTat(i1 - 1, k - i1, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULTat(i1 - 1, k - i1, m + 1) - lscf.ULTat(i1, k - i1 - 1, m + 1)));
            tG = s1 * A1 * (dks.ULTat(i1 - 1, k - i1, m + 1) + Go.col(i1 - 1)) +
                 s2 * A2 * (dks.ULTat(i1, k - i1 - 1, m + 1) + Go.col(i1));
            gn.col(i1) =
                (tG - s1 * Go.col(i1 - 1) -
                 s2 * Go.col(i1)
                 - ((s1 * dks.ULTat(i1 - 1, k - i1, m + 1) + s2 * dks.ULTat(i1, k - i1 - 1, m + 1)))) * mu +
                s1 * A1 * go.col(i1 - 1) + s2 * A2 * go.col(i1);
            Gn.col(i1) = tG;
            dks.ULTat(i1, k - i1, m + 1) = (Gn.col(i1).sum() + (mu * gn.col(i1)).sum()) / (2 * k);
            scale_in_h2_ij_vE(i1, k, m, n, thr, dks, lscf, Gn, gn);
        }
#ifdef _OPENMP
}
#endif
        tG = A1 * (dks.ULTat(k - 1, 0, m + 1) + Go.col(k - 1));
        gn.col(k) = (tG - Go.col(k - 1)
             - ((dks.ULTat(k - 1, 0, m + 1)))) * mu + A1 * go.col(k - 1);
        Gn.col(k) = tG;
        dks.ULTat(k, 0, m + 1) = (Gn.col(k).sum() + (mu * gn.col(k)).sum()) / (2 * k);
        scale_in_h2_ij_vE(k, k, m, n, thr, dks, lscf, Gn, gn);
    }
    return dks;
}
template ArrayXd h2_ij_vEc(const ArrayBase<ArrayXd>& A1,
                           const ArrayBase<ArrayXd>& A2,
                           const ArrayBase<ArrayXd>& mu,
                           const Index m, ArrayXd& lscf,
                           const double thr_margin, int nthreads);


template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_d3_ijk_mE(Eigen::Index i1, Eigen::Index i2, Eigen::Index k,
                               const Eigen::Index m, const Eigen::Index n, 
                               typename DerivedA::Scalar &thr,
                               Eigen::ArrayBase<DerivedA> &dks,
                               Eigen::ArrayBase<DerivedB> &lscf,
                               Eigen::MatrixBase<DerivedC> &Gn) {
    Index i3 = k - i1 - i2;
    if(Gn.block(0, n * id3(i1, i2, k), n, n).maxCoeff() > thr) {
        dks.ULCat(i1, i2, i3, m + 1) /= 1e10;
        Gn.block(0, n * id3(i1, i2, k), n, n) /= 1e10;
        update_scale_3D(lscf, i1, i2, i3, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d3_ijk_mEc(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    dks.ULCat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx Go = MatrixXx::Zero(n, n * (m + 1) * m / 2);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 2) * (m + 1) / 2);
    Gn.block(0, 0, n, n).diagonal().array() += dks.ULCat(0, 0, 0, m + 1);
    Scalar s1, s2, s3, min_lscf;
    for(Index k = 1; k <= m; k++) {
        if(k % 50 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k * (k + 1) / 2) = Gn.block(0, 0, n, n * k * (k + 1) / 2);

        MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =
            A3 * Go.block(0, 0, n, n);
        dks.ULCat(0, 0, k, m + 1) = Gn.block(0, 0, n, n).trace() / (2 * k);
        Gn.block(0, 0, n, n).diagonal().array() += dks.ULCat(0, 0, k, m + 1);
        scale_in_d3_ijk_mE(0, 0, k, m, n, thr, dks, lscf, Gn);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            s2 = exp(min<Scalar>(0, lscf.ULCat(0, i2, i3 - 1, m + 1) - lscf.ULCat(0, i2 - 1, i3, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(0, i2 - 1, i3, m + 1) - lscf.ULCat(0, i2, i3 - 1, m + 1)));
            MatrixXx::Map(Gn.block(0, n * id3(0, i2, k), n, n).data(), n, n).noalias() =
                s2 * A2 * Go.block(0, n * id3(0, i2 - 1, k - 1), n, n) +
                s3 * A3 * Go.block(0, n * id3(0, i2, k - 1), n, n);
            dks.ULCat(0, i2, i3, m + 1) = Gn.block(0, n * id3(0, i2, k), n, n).trace() / (2 * k);
            Gn.block(0, n * id3(0, i2, k), n, n).diagonal().array() += dks.ULCat(0, i2, i3, m + 1);
            scale_in_d3_ijk_mE(0, i2, k, m, n, thr, dks, lscf, Gn);
        }
        MatrixXx::Map(Gn.block(0, n * id3(0, k, k), n, n).data(), n, n).noalias() =
            A2 * Go.block(0, n * id3(0, k - 1, k - 1), n, n);
        dks.ULCat(0, k, 0, m + 1) = Gn.block(0, n * id3(0, k, k), n, n).trace() / (2 * k);
        Gn.block(0, n * id3(0, k, k), n, n).diagonal().array() += dks.ULCat(0, k, 0, m + 1);
        scale_in_d3_ijk_mE(0, k, k, m, n, thr, dks, lscf, Gn);
#ifdef _OPENMP
#pragma omp parallel private(s1, s2, s3, min_lscf)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, 0, k - i1 - 1, m + 1) - lscf.ULCat(i1 - 1, 0, k - i1, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, 0, k - i1, m + 1) - lscf.ULCat(i1, 0, k - i1 - 1, m + 1)));
            MatrixXx::Map(Gn.block(0, n * i1, n, n).data(), n, n).noalias() =
                s1 * A1 * Go.block(0, n * (i1 - 1), n, n) +
                s3 * A3 * Go.block(0, n * i1, n, n);
            dks.ULCat(i1, 0, k - i1, m + 1) = Gn.block(0, n * i1, n, n).trace() / (2 * k);
            Gn.block(0, n * i1, n, n).diagonal().array() += dks.ULCat(i1, 0, k - i1, m + 1);
            scale_in_d3_ijk_mE(i1, 0, k, m, n, thr, dks, lscf, Gn);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                min_lscf = min<Scalar>({lscf.ULCat(i1 - 1, i2, i3, m + 1),
                                        lscf.ULCat(i1, i2 - 1, i3, m + 1),
                                        lscf.ULCat(i1, i2, i3 - 1, m + 1)});
                s1 = exp(min_lscf - lscf.ULCat(i1 - 1, i2, i3, m + 1));
                s2 = exp(min_lscf - lscf.ULCat(i1, i2 - 1, i3, m + 1));
                s3 = exp(min_lscf - lscf.ULCat(i1, i2, i3 - 1, m + 1));
                MatrixXx::Map(Gn.block(0, n * id3(i1, i2, k), n, n).data(), n, n).noalias() =
                    s1 * A1 * Go.block(0, n * id3(i1 - 1, i2, k - 1), n, n) +
                    s2 * A2 * Go.block(0, n * id3(i1, i2 - 1, k - 1), n, n) +
                    s3 * A3 * Go.block(0, n * id3(i1, i2, k - 1), n, n);
                dks.ULCat(i1, i2, i3, m + 1) = Gn.block(0, n * id3(i1, i2, k), n, n).trace() / (2 * k);
                Gn.block(0, n * id3(i1, i2, k), n, n).diagonal().array() += dks.ULCat(i1, i2, i3, m + 1);
                scale_in_d3_ijk_mE(i1, i2, k, m, n, thr, dks, lscf, Gn);
            }
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, k - i1 - 1, 0, m + 1) - lscf.ULCat(i1 - 1, k - i1, 0, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, k - i1, 0, m + 1) - lscf.ULCat(i1, k - i1 - 1, 0, m + 1)));
            MatrixXx::Map(Gn.block(0, n * id3(i1, k - i1, k), n, n).data(), n, n).noalias() =
                s1 * A1 * Go.block(0, n * id3(i1 - 1, k - i1, k - 1), n, n) +
                s2 * A2 * Go.block(0, n * id3(i1, k - i1 - 1, k - 1), n, n);
            dks.ULCat(i1, k - i1, 0, m + 1) = Gn.block(0, n * id3(i1, k - i1, k), n, n).trace() / (2 * k);
            Gn.block(0, n * id3(i1, k - i1, k), n, n).diagonal().array() += dks.ULCat(i1, k - i1, 0, m + 1);
            scale_in_d3_ijk_mE(i1, k - i1, k, m, n, thr, dks, lscf, Gn);
        }
#ifdef _OPENMP
}
#endif
        //
        MatrixXx::Map(Gn.block(0, n * k, n, n).data(), n, n).noalias() =
            A1 * Go.block(0, n * (k - 1), n, n);
        dks.ULCat(k, 0, 0, m + 1) = Gn.block(0, n * k, n, n).trace() / (2 * k);
        Gn.block(0, n * k, n, n).diagonal().array() += dks.ULCat(k, 0, 0, m + 1);
        scale_in_d3_ijk_mE(k, 0, k, m, n, thr, dks, lscf, Gn);
    }
    return dks;
}
template ArrayXd d3_ijk_mEc(const MatrixBase<MatrixXd>& A1,
                            const DiagMatXd& A2,
                            const MatrixBase<MatrixXd>& A3,
                            const Index m, ArrayXd& lscf,
                            const double thr_margin, int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_d3_ijk_vE(Eigen::Index i1, Eigen::Index i2, Eigen::Index k,
                               const Eigen::Index m, const Eigen::Index n, 
                               typename DerivedA::Scalar &thr,
                               Eigen::ArrayBase<DerivedA> &dks,
                               Eigen::ArrayBase<DerivedB> &lscf,
                               Eigen::ArrayBase<DerivedC> &Gn) {
    Index i3 = k - i1 - i2;
    if(Gn.col(id3(i1, i2, k)).maxCoeff() > thr) {
        dks.ULCat(i1, i2, i3, m + 1) /= 1e10;
        Gn.col(id3(i1, i2, k)) /= 1e10;
        update_scale_3D(lscf, i1, i2, i3, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d3_ijk_vEc(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    dks.ULCat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (m + 1) * m / 2);
    ArrayXXx Gn = ArrayXXx::Zero(n, (m + 2) * (m + 1) / 2);
    Scalar s1, s2, s3, min_lscf;
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (k + 1) / 2) = Gn.block(0, 0, n, k * (k + 1) / 2);

        Gn.col(0) = A3 * (dks.ULCat(0, 0, k - 1, m + 1) + Go.col(0));
        dks.ULCat(0, 0, k, m + 1) = Gn.col(0).sum() / (2 * k);
        scale_in_d3_ijk_vE(0, 0, k, m, n, thr, dks, lscf, Gn);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            s2 = exp(min<Scalar>(0, lscf.ULCat(0, i2, i3 - 1, m + 1) - lscf.ULCat(0, i2 - 1, i3, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(0, i2 - 1, i3, m + 1) - lscf.ULCat(0, i2, i3 - 1, m + 1)));
            Gn.col(id3(0, i2, k)) =
                s2 * A2 * (dks.ULCat(0, i2 - 1, i3, m + 1) + Go.col(id3(0, i2 - 1, k - 1))) +
                s3 * A3 * (dks.ULCat(0, i2, i3 - 1, m + 1) + Go.col(id3(0, i2, k - 1)));
            dks.ULCat(0, i2, i3, m + 1) = Gn.col(id3(0, i2, k)).sum() / (2 * k);
            scale_in_d3_ijk_vE(0, i2, k, m, n, thr, dks, lscf, Gn);
        }
        Gn.col(id3(0, k, k)) = A2 * (dks.ULCat(0, k - 1, 0, m + 1) + Go.col(id3(0, k - 1, k - 1)));
        dks.ULCat(0, k, 0, m + 1) = Gn.col(id3(0, k, k)).sum() / (2 * k);
        scale_in_d3_ijk_vE(0, k, k, m, n, thr, dks, lscf, Gn);
#ifdef _OPENMP
#pragma omp parallel private(s1, s2, s3, min_lscf)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, 0, k - i1 - 1, m + 1) - lscf.ULCat(i1 - 1, 0, k - i1, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, 0, k - i1, m + 1) - lscf.ULCat(i1, 0, k - i1 - 1, m + 1)));
            Gn.col(i1) =
                s1 * A1 * (dks.ULCat(i1 - 1, 0, k - i1, m + 1) + Go.col(i1 - 1)) +
                s3 * A3 * (dks.ULCat(i1, 0, k - i1 - 1, m + 1) + Go.col(i1));
            dks.ULCat(i1, 0, k - i1, m + 1) = Gn.col(i1).sum() / (2 * k);
            scale_in_d3_ijk_vE(i1, 0, k, m, n, thr, dks, lscf, Gn);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                min_lscf = min<Scalar>({lscf.ULCat(i1 - 1, i2, i3, m + 1),
                                        lscf.ULCat(i1, i2 - 1, i3, m + 1),
                                        lscf.ULCat(i1, i2, i3 - 1, m + 1)});
                s1 = exp(min_lscf - lscf.ULCat(i1 - 1, i2, i3, m + 1));
                s2 = exp(min_lscf - lscf.ULCat(i1, i2 - 1, i3, m + 1));
                s3 = exp(min_lscf - lscf.ULCat(i1, i2, i3 - 1, m + 1));
                Gn.col(id3(i1, i2, k)) =
                    s1 * A1 * (dks.ULCat(i1 - 1, i2, i3, m + 1) + Go.col(id3(i1 - 1, i2, k - 1))) +
                    s2 * A2 * (dks.ULCat(i1, i2 - 1, i3, m + 1) + Go.col(id3(i1, i2 - 1, k - 1))) +
                    s3 * A3 * (dks.ULCat(i1, i2, i3 - 1, m + 1) + Go.col(id3(i1, i2, k - 1)));
                dks.ULCat(i1, i2, i3, m + 1) = Gn.col(id3(i1, i2, k)).sum() / (2 * k);
                scale_in_d3_ijk_vE(i1, i2, k, m, n, thr, dks, lscf, Gn);
            }
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, k - i1 - 1, 0, m + 1) - lscf.ULCat(i1 - 1, k - i1, 0, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, k - i1, 0, m + 1) - lscf.ULCat(i1, k - i1 - 1, 0, m + 1)));
            Gn.col(id3(i1, k - i1, k)) =
                s1 * A1 * (dks.ULCat(i1 - 1, k - i1, 0, m + 1) + Go.col(id3(i1 - 1, k - i1, k - 1))) +
                s2 * A2 * (dks.ULCat(i1, k - i1 - 1, 0, m + 1) + Go.col(id3(i1, k - i1 - 1, k - 1)));
            dks.ULCat(i1, k - i1, 0, m + 1) = Gn.col(id3(i1, k - i1, k)).sum() / (2 * k);
            scale_in_d3_ijk_vE(i1, k - i1, k, m, n, thr, dks, lscf, Gn);
        }
#ifdef _OPENMP
}
#endif
        Gn.col(k) = A1 * (dks.ULCat(k - 1, 0, 0, m + 1) + Go.col(k - 1));
        dks.ULCat(k, 0, 0, m + 1) = Gn.col(k).sum() / (2 * k);
        scale_in_d3_ijk_vE(k, 0, k, m, n, thr, dks, lscf, Gn);
    }
    return dks;
}
template ArrayXd d3_ijk_vEc(const ArrayBase<ArrayXd>& A1,
                            const ArrayBase<ArrayXd>& A2,
                            const ArrayBase<ArrayXd>& A3,
                            const Index m, ArrayXd& lscf,
                            const double thr_margin, int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_d3_pjk_mE(Eigen::Index j, Eigen::Index k,
                               const Eigen::Index m, const Eigen::Index n,
                               const Eigen::Index p,
                               typename DerivedA::Scalar &thr,
                               Eigen::ArrayBase<DerivedA> &dks,
                               Eigen::ArrayBase<DerivedB> &lscf,
                               Eigen::MatrixBase<DerivedC> &Gn) {
    if(Gn.block(0, j * n * (p + 1), n, n * (p + 1)).maxCoeff() > thr) {
        dks.ULScol(k - j, j, m + 1) /= 1e10;
        Gn.block(0, j * n * (p + 1), n, n * (p + 1)) /= 1e10;
        update_scale_2D(lscf, k - j, j, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Index m,
          const Eigen::Index p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 2) / 2);
    dks.ULSat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx Go = MatrixXx::Zero(n, n * (p + 1) * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (p + 1) * (m + 1));
    Scalar s2, s3;
    Gn.block(0, 0, n, n).diagonal().array() = dks.ULSat(0, 0, 0, m + 1);
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
            A1 * Gn.block(0, (i - 1) * n, n, n);
        dks.ULSat(i, 0, 0, m + 1) = Gn.block(0, i * n, n, n).trace() / (2 * i);
        Gn.block(0, i * n, n, n).diagonal().array() += dks.ULSat(i, 0, 0, m + 1);
    }
    scale_in_d3_pjk_mE(0, 0, m, n, p, thr, dks, lscf, Gn);
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        MatrixXx::Map(Gn.block(0, 0, n, n).data(), n, n).noalias() =
            A2 * Go.block(0, 0, n, n);
        dks.ULSat(0, k, 0, m + 1) = Gn.block(0, 0, n, n).trace() / (2 * k);
        Gn.block(0, 0, n, n).diagonal().array() += dks.ULSat(0, k, 0, m + 1);
        for(Index i = 1; i <= p; i++) {
            MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
                A1 * Gn.block(0, (i - 1) * n, n, n) +
                A2 * Go.block(0, i * n, n, n);
            dks.ULSat(i, k, 0, m + 1) = Gn.block(0, i * n, n, n).trace() / (2 * (k + i));
            Gn.block(0, i * n, n, n).diagonal().array() += dks.ULSat(i, k, 0, m + 1);
        }
        scale_in_d3_pjk_mE(0, k, m, n, p, thr, dks, lscf, Gn);
#ifdef _OPENMP
#pragma omp parallel private(s2, s3)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            s2 = exp(min<Scalar>(0, lscf.ULTat(k - j, j - 1, m + 1) - lscf.ULTat(k - j - 1, j, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULTat(k - j - 1, j, m + 1) - lscf.ULTat(k - j, j - 1, m + 1)));
            MatrixXx::Map(Gn.block(0, j * n * (p + 1), n, n).data(), n, n).noalias() =
                s2 * A2 * Go.block(0, j * n * (p + 1), n, n) +
                s3 * A3 * Go.block(0, (j - 1) * n * (p + 1), n, n);
            dks.ULSat(0, k - j, j, m + 1) = Gn.block(0, j * n * (p + 1), n, n).trace() / (2 * k);
            Gn.block(0, j * n * (p + 1), n, n).diagonal().array() += dks.ULSat(0, k - j, j, m + 1);
            for(Index i = 1; i <= p; i++) {
                MatrixXx::Map(Gn.block(0, j * n * (p + 1) + i * n, n, n).data(), n, n).noalias() =
                         A1 * Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n) +
                    s2 * A2 * Go.block(0, j * n * (p + 1) + i * n, n, n) +
                    s3 * A3 * Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n);
                dks.ULSat(i, k - j, j, m + 1) = Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() / (2 * (k + i));
                Gn.block(0, j * n * (p + 1) + i * n, n, n).diagonal().array() += dks.ULSat(i, k - j, j, m + 1);
            }
            scale_in_d3_pjk_mE(j, k, m, n, p, thr, dks, lscf, Gn);
        }
#ifdef _OPENMP
}
#endif
        MatrixXx::Map(Gn.block(0, k * n * (p + 1), n, n).data(), n, n).noalias() =
            A3 * Go.block(0, (k - 1) * n * (p + 1), n, n);
        dks.ULSat(0, 0, k, m + 1) = Gn.block(0, k * n * (p + 1), n, n).trace() / (2 * k);
        Gn.block(0, k * n * (p + 1), n, n).diagonal().array() += dks.ULSat(0, 0, k, m + 1);
        for(Index i = 1; i <= p; i++) {
            MatrixXx::Map(Gn.block(0, k * n * (p + 1) + i * n, n, n).data(), n, n).noalias() =
                A1 * Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n) +
                A3 * Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n);
            dks.ULSat(i, 0, k, m + 1) = Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() / (2 * (k + i));
            Gn.block(0, k * n * (p + 1) + i * n, n, n).diagonal().array() += dks.ULSat(i, 0, k, m + 1);
        }
        scale_in_d3_pjk_mE(k, k, m, n, p, thr, dks, lscf, Gn);
    }
    return dks;
}
template ArrayXXd d3_pjk_mEc(const MatrixBase<MatrixXd>& A1,
                const DiagMatXd& A2, const MatrixBase<MatrixXd>& A3,
                const Index m, const Index p, ArrayXd& lscf,
                const double thr_margin, int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_d3_pjk_vE(Eigen::Index j, Eigen::Index k,
                               const Eigen::Index m, const Eigen::Index n,
                               const Eigen::Index p,
                               typename DerivedA::Scalar &thr,
                               Eigen::ArrayBase<DerivedA> &dks,
                               Eigen::ArrayBase<DerivedB> &lscf,
                               Eigen::ArrayBase<DerivedC> &Gn) {
    if(Gn.block(0, j * (p + 1), n, p + 1).maxCoeff() > thr) {
        dks.ULScol(k - j, j, m + 1) /= 1e10;
        Gn.block(0, j * (p + 1), n, p + 1) /= 1e10;
        update_scale_2D(lscf, k - j, j, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::Index m, const Eigen::Index p,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 2) / 2);
    dks.ULSat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx Gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    Scalar s2, s3;
    for(Index i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks.ULSat(i - 1, 0, 0, m + 1) + Gn.col(i - 1));
        dks.ULSat(i, 0, 0, m + 1) = Gn.col(i).sum() / (2 * i);
    }
    scale_in_d3_pjk_vE(0, 0, m, n, p, thr, dks, lscf, Gn);
    for(Index k = 1; k <= m; k++) {
        if(k % 500 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        Gn.col(0) = A2 * (dks.ULSat(0, k - 1, 0, m + 1) + Go.col(0));
        dks.ULSat(0, k, 0, m + 1) = Gn.col(0).sum() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            Gn.col(i) = A1 * (dks.ULSat(i - 1, k, 0, m + 1) + Gn.col(i - 1)) +
                        A2 * (dks.ULSat(i, k - 1, 0, m + 1) + Go.col(i));
            dks.ULSat(i, k, 0, m + 1) = Gn.col(i).sum() / (2 * (k + i));
        }
        scale_in_d3_pjk_vE(0, k, m, n, p, thr, dks, lscf, Gn);
#ifdef _OPENMP
#pragma omp parallel private(s2, s3)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            s2 = exp(min<Scalar>(0, lscf.ULTat(k - j, j - 1, m + 1) - lscf.ULTat(k - j - 1, j, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULTat(k - j - 1, j, m + 1) - lscf.ULTat(k - j, j - 1, m + 1)));
            Gn.col(j * (p + 1)) =
                s2 * A2 * (dks.ULSat(0, k - j - 1, j, m + 1) + Go.col(j * (p + 1))) +
                s3 * A3 * (dks.ULSat(0, k - j, j - 1, m + 1) + Go.col((j - 1) * (p + 1)));
            dks.ULSat(0, k - j, j, m + 1) = Gn.col(j * (p + 1)).sum() / (2 * k);
            for(Index i = 1; i <= p; i++) {
                Gn.col(j * (p + 1) + i) =
                         A1 * (dks.ULSat(i - 1, k - j, j, m + 1) + Gn.col(j * (p + 1) + (i - 1))) +
                    s2 * A2 * (dks.ULSat(i, k - j - 1, j, m + 1) + Go.col(j * (p + 1) + i)) +
                    s3 * A3 * (dks.ULSat(i, k - j, j - 1, m + 1) + Go.col((j - 1) * (p + 1) + i));
                dks.ULSat(i, k - j, j, m + 1) = Gn.col(j * (p + 1) + i).sum() / (2 * (k + i));
            }
            scale_in_d3_pjk_vE(j, k, m, n, p, thr, dks, lscf, Gn);
        }
#ifdef _OPENMP
}
#endif
        Gn.col(k * (p + 1)) = A3 * (dks.ULSat(0, 0, k - 1, m + 1) + Go.col((k - 1) * (p + 1)));
        dks.ULSat(0, 0, k, m + 1) = Gn.col(k * (p + 1)).sum() / (2 * k);
        for(Index i = 1; i <= p; i++) {
            Gn.col(k * (p + 1) + i) =
                A1 * (dks.ULSat(i - 1, 0, k, m + 1) + Gn.col(k * (p + 1) + i - 1)) +
                A3 * (dks.ULSat(i, 0, k - 1, m + 1) + Go.col((k - 1) * (p + 1) + i));
            dks.ULSat(i, 0, k, m + 1) = Gn.col(k * (p + 1) + i).sum() / (2 * (k + i));
        }
        scale_in_d3_pjk_vE(k, k, m, n, p, thr, dks, lscf, Gn);
    }
    return dks;
}
template ArrayXXd d3_pjk_vEc(const ArrayBase<ArrayXd>& A1,
                            const ArrayBase<ArrayXd>& A2,
                            const ArrayBase<ArrayXd>& A3,
                            const Index m, const Index p, ArrayXd& lscf,
                            const double thr_margin, int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
inline void scale_in_h3_ijk_mE(Eigen::Index i1, Eigen::Index i2, Eigen::Index k,
                               const Eigen::Index m, const Eigen::Index n,
                               typename DerivedA::Scalar &thr,
                               Eigen::ArrayBase<DerivedA> &dks,
                               Eigen::ArrayBase<DerivedB> &lscf,
                               Eigen::MatrixBase<DerivedC> &Gn,
                               Eigen::MatrixBase<DerivedD> &gn) {
    Index i3 = k - i1 - i2;
    if(Gn.block(0, n * id3(i1, i2, k), n, n).maxCoeff() > thr || gn.col(id3(i1, i2, k)).maxCoeff() > thr) {
        dks.ULCat(i1, i2, i3, m + 1) /= 1e10;
        Gn.block(0, n * id3(i1, i2, k), n, n) /= 1e10;
        gn.col(id3(i1, i2, k)) /= 1e10;
        update_scale_3D(lscf, i1, i2, i3, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h3_ijk_mEc(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    dks.ULCat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * (m + 1) * m / 2);
    MatrixXx Gn = MatrixXx::Zero(n, n * (m + 2) * (m + 1) / 2);
    MatrixXx go = MatrixXx::Zero(n, (m + 1) * m / 2);
    MatrixXx gn = MatrixXx::Zero(n, (m + 2) * (m + 1) / 2);
    Gn.block(0, 0, n, n).diagonal().array() += dks.ULCat(0, 0, 0, m + 1);
    Scalar s1, s2, s3, min_lscf;
    for(Index k = 1; k <= m; k++) {
        if(k % 50 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, n * k * (k + 1) / 2) = Gn.block(0, 0, n, n * k * (k + 1) / 2);
        go.block(0, 0, n, k * (k + 1) / 2) = gn.block(0, 0, n, k * (k + 1) / 2);

        tG.noalias() = A3 * Go.block(0, 0, n, n);
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG - Go.block(0, 0, n, n)) * mu + A3 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks.ULCat(0, 0, k, m + 1) = (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        Gn.block(0, 0, n, n).diagonal().array() += dks.ULCat(0, 0, k, m + 1);
        scale_in_h3_ijk_mE(0, 0, k, m, n, thr, dks, lscf, Gn, gn);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            s2 = exp(min<Scalar>(0, lscf.ULCat(0, i2, i3 - 1, m + 1) - lscf.ULCat(0, i2 - 1, i3, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(0, i2 - 1, i3, m + 1) - lscf.ULCat(0, i2, i3 - 1, m + 1)));
            tG.noalias() = s2 * A2 * Go.block(0, n * id3(0, i2 - 1, k - 1), n, n) +
                           s3 * A3 * Go.block(0, n * id3(0, i2, k - 1), n, n);
            MatrixXx::Map(gn.col(id3(0, i2, k)).data(), n, 1).noalias() =
                (tG -
                 s2 * Go.block(0, n * id3(0, i2 - 1, k - 1), n, n) - s3 * Go.block(0, n * id3(0, i2, k - 1), n, n)) * mu +
                s2 * A2 * go.col(id3(0, i2 - 1, k - 1)) + s3 * A3 * go.col(id3(0, i2, k - 1));
            Gn.block(0, n * id3(0, i2, k), n, n) = tG;
            dks.ULCat(0, i2, i3, m + 1) = (Gn.block(0, n * id3(0, i2, k), n, n).trace() + gn.col(id3(0, i2, k)).dot(mu)) / (2 * k);
            Gn.block(0, n * id3(0, i2, k), n, n).diagonal().array() += dks.ULCat(0, i2, i3, m + 1);
            scale_in_h3_ijk_mE(0, i2, k, m, n, thr, dks, lscf, Gn, gn);
        }
        tG.noalias() = A2 * Go.block(0, n * id3(0, k - 1, k - 1), n, n);
        MatrixXx::Map(gn.col(id3(0, k, k)).data(), n, 1).noalias() = (tG - Go.block(0, n * id3(0, k - 1, k - 1), n, n)) * mu + A2 * go.col(id3(0, k - 1, k - 1));
        Gn.block(0, n * id3(0, k, k), n, n) = tG;
        dks.ULCat(0, k, 0, m + 1) = (Gn.block(0, n * id3(0, k, k), n, n).trace() + gn.col(id3(0, k, k)).dot(mu)) / (2 * k);
        Gn.block(0, n * id3(0, k, k), n, n).diagonal().array() += dks.ULCat(0, k, 0, m + 1);
        scale_in_h3_ijk_mE(0, k, k, m, n, thr, dks, lscf, Gn, gn);
#ifdef _OPENMP
#pragma omp parallel private(tG, s1, s2, s3, min_lscf)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, 0, k - i1 - 1, m + 1) - lscf.ULCat(i1 - 1, 0, k - i1, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, 0, k - i1, m + 1) - lscf.ULCat(i1, 0, k - i1 - 1, m + 1)));
            tG.noalias() = s1 * A1 * Go.block(0, n * (i1 - 1), n, n) +
                           s3 * A3 * Go.block(0, n * i1, n, n);
            MatrixXx::Map(gn.col(i1).data(), n, 1).noalias() =
                (tG - s1 * Go.block(0, n * (i1 - 1), n, n) - s3 * Go.block(0, n * i1, n, n)) * mu +
                s1 * A1 * go.col(i1 - 1) + s3 * A3 * go.col(i1);
            Gn.block(0, n * i1, n, n) = tG;
            dks.ULCat(i1, 0, k - i1, m + 1) = (Gn.block(0, n * i1, n, n).trace() + gn.col(i1).dot(mu)) / (2 * k);
            Gn.block(0, n * i1, n, n).diagonal().array() += dks.ULCat(i1, 0, k - i1, m + 1);
            scale_in_h3_ijk_mE(i1, 0, k, m, n, thr, dks, lscf, Gn, gn);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                min_lscf = min<Scalar>({lscf.ULCat(i1 - 1, i2, i3, m + 1), 
                                        lscf.ULCat(i1, i2 - 1, i3, m + 1), 
                                        lscf.ULCat(i1, i2, i3 - 1, m + 1)});
                s1 = exp(min_lscf - lscf.ULCat(i1 - 1, i2, i3, m + 1));
                s2 = exp(min_lscf - lscf.ULCat(i1, i2 - 1, i3, m + 1));
                s3 = exp(min_lscf - lscf.ULCat(i1, i2, i3 - 1, m + 1));
                tG.noalias() = s1 * A1 * Go.block(0, n * id3(i1 - 1, i2, k - 1), n, n) +
                               s2 * A2 * Go.block(0, n * id3(i1, i2 - 1, k - 1), n, n) +
                               s3 * A3 * Go.block(0, n * id3(i1, i2, k - 1), n, n);
                MatrixXx::Map(gn.col(id3(i1, i2, k)).data(), n, 1).noalias() =
                    (tG - s1 * Go.block(0, n * id3(i1 - 1, i2, k - 1), n, n) -
                     s2 * Go.block(0, n * id3(i1, i2 - 1, k - 1), n, n) - s3 * Go.block(0, n * id3(i1, i2, k - 1), n, n)) * mu +
                    s1 * A1 * go.col(id3(i1 - 1, i2, k - 1)) + s2 * A2 * go.col(id3(i1, i2 - 1, k - 1)) + s3 * A3 * go.col(id3(i1, i2, k - 1));
                Gn.block(0, n * id3(i1, i2, k), n, n) = tG;
                dks.ULCat(i1, i2, i3, m + 1) = (Gn.block(0, n * id3(i1, i2, k), n, n).trace() + gn.col(id3(i1, i2, k)).dot(mu)) / (2 * k);
                Gn.block(0, n * id3(i1, i2, k), n, n).diagonal().array() += dks.ULCat(i1, i2, i3, m + 1);
                scale_in_h3_ijk_mE(i1, i2, k, m, n, thr, dks, lscf, Gn, gn);
            }
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, k - i1 - 1, 0, m + 1) - lscf.ULCat(i1 - 1, k - i1, 0, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, k - i1, 0, m + 1) - lscf.ULCat(i1, k - i1 - 1, 0, m + 1)));
            tG.noalias() = s1 * A1 * Go.block(0, n * id3(i1 - 1, k - i1, k - 1), n, n) +
                           s2 * A2 * Go.block(0, n * id3(i1, k - i1 - 1, k - 1), n, n);
            MatrixXx::Map(gn.col(id3(i1, k - i1, k)).data(), n, 1).noalias() =
                (tG - s1 * Go.block(0, n * id3(i1 - 1, k - i1, k - 1), n, n) -
                 s2 * Go.block(0, n * id3(i1, k - i1 - 1, k - 1), n, n)) * mu +
                s1 * A1 * go.col(id3(i1 - 1, k - i1, k - 1)) + s2 * A2 * go.col(id3(i1, k - i1 - 1, k - 1));
            Gn.block(0, n * id3(i1, k - i1, k), n, n) = tG;
            dks.ULCat(i1, k - i1, 0, m + 1) = (Gn.block(0, n * id3(i1, k - i1, k), n, n).trace() + gn.col(id3(i1, k - i1, k)).dot(mu)) / (2 * k);
            Gn.block(0, n * id3(i1, k - i1, k), n, n).diagonal().array() += dks.ULCat(i1, k - i1, 0, m + 1);
            scale_in_h3_ijk_mE(i1, k - i1, k, m, n, thr, dks, lscf, Gn, gn);
        }
#ifdef _OPENMP
}
#endif
        //
        tG.noalias() = A1 * Go.block(0, n * (k - 1), n, n);
        MatrixXx::Map(gn.col(k).data(), n, 1).noalias() =
            (tG - Go.block(0, n * (k - 1), n, n)) * mu + A1 * go.col(k - 1);
        Gn.block(0, n * k, n, n) = tG;
        dks.ULCat(k, 0, 0, m + 1) = (Gn.block(0, n * k, n, n).trace() + gn.col(k).dot(mu)) / (2 * k);
        Gn.block(0, n * k, n, n).diagonal().array() += dks.ULCat(k, 0, 0, m + 1);
        scale_in_h3_ijk_mE(k, 0, k, m, n, thr, dks, lscf, Gn, gn);
    }
    return dks;
}
template ArrayXd  h3_ijk_mEc(const MatrixBase<MatrixXd>& A1,
                             const DiagMatXd& A2,
                             const MatrixBase<MatrixXd>& A3,
                             const VectorXd mu, const Index m,
                             ArrayXd& lscf, const double thr_margin,
                             int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_h3_ijk_vE(Eigen::Index i1, Eigen::Index i2, Eigen::Index k,
                               const Eigen::Index m, const Eigen::Index n,
                               typename DerivedA::Scalar &thr,
                               Eigen::ArrayBase<DerivedA> &dks,
                               Eigen::ArrayBase<DerivedB> &lscf,
                               Eigen::ArrayBase<DerivedC> &Gn,
                               Eigen::ArrayBase<DerivedC> &gn) {
    Index i3 = k - i1 - i2;
    if(Gn.col(id3(i1, i2, k)).maxCoeff() > thr || gn.col(id3(i1, i2, k)).maxCoeff() > thr) {
        dks.ULCat(i1, i2, i3, m + 1) /= 1e10;
        Gn.col(id3(i1, i2, k)) /= 1e10;
        gn.col(id3(i1, i2, k)) /= 1e10;
        update_scale_3D(lscf, i1, i2, i3, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h3_ijk_vEc(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::ArrayBase<Derived>& mu, const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXx dks = ArrayXx::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    dks.ULCat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (m + 1) * m / 2);
    ArrayXXx Gn = ArrayXXx::Zero(n, (m + 2) * (m + 1) / 2);
    ArrayXXx go = ArrayXXx::Zero(n, (m + 1) * m / 2);
    ArrayXXx gn = ArrayXXx::Zero(n, (m + 2) * (m + 1) / 2);
    Scalar s1, s2, s3, min_lscf;
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (k + 1) / 2) = Gn.block(0, 0, n, k * (k + 1) / 2);
        go.block(0, 0, n, k * (k + 1) / 2) = gn.block(0, 0, n, k * (k + 1) / 2);

        tG = A3 * (dks.ULCat(0, 0, k - 1, m + 1) + Go.col(0));
        gn.col(0) = (tG - Go.col(0)
             - (dks.ULCat(0, 0, k - 1, m + 1))) * mu + A3 * go.col(0);
        Gn.col(0) = tG;
        dks.ULCat(0, 0, k, m + 1) = (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        scale_in_h3_ijk_vE(0, 0, k, m, n, thr, dks, lscf, Gn, gn);
        for(Index i2 = 1; i2 < k; i2++) {
            Index i3 = k - i2;
            s2 = exp(min<Scalar>(0, lscf.ULCat(0, i2, i3 - 1, m + 1) - lscf.ULCat(0, i2 - 1, i3, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(0, i2 - 1, i3, m + 1) - lscf.ULCat(0, i2, i3 - 1, m + 1)));
            tG = s2 * A2 * (dks.ULCat(0, i2 - 1, i3, m + 1) + Go.col(id3(0, i2 - 1, k - 1))) +
                 s3 * A3 * (dks.ULCat(0, i2, i3 - 1, m + 1) + Go.col(id3(0, i2, k - 1)));
            gn.col(id3(0, i2, k)) =
                (tG - s2 * Go.col(id3(0, i2 - 1, k - 1)) - s3 * Go.col(id3(0, i2, k - 1))
                 - ((s2 * dks.ULCat(0, i2 - 1, i3, m + 1) +
                    s3 * dks.ULCat(0, i2, i3 - 1, m + 1)))) * mu +
                s2 * A2 * go.col(id3(0, i2 - 1, k - 1)) + s3 * A3 * go.col(id3(0, i2, k - 1));
            Gn.col(id3(0, i2, k)) = tG;
            dks.ULCat(0, i2, i3, m + 1) = (Gn.col(id3(0, i2, k)).sum() + (mu * gn.col(id3(0, i2, k))).sum()) / (2 * k);
            scale_in_h3_ijk_vE(0, i2, k, m, n, thr, dks, lscf, Gn, gn);
        }
        tG = A2 * (dks.ULCat(0, k - 1, 0, m + 1) + Go.col(id3(0, k - 1, k - 1)));
        gn.col(id3(0, k, k)) = (tG - Go.col(id3(0, k - 1, k - 1))
             - ((dks.ULCat(0, k - 1, 0, m + 1)))) * mu + A2 * go.col(id3(0, k - 1, k - 1));
        Gn.col(id3(0, k, k)) = tG;
        dks.ULCat(0, k, 0, m + 1) = (Gn.col(id3(0, k, k)).sum() + (mu * gn.col(id3(0, k, k))).sum()) / (2 * k);
        scale_in_h3_ijk_vE(0, k, k, m, n, thr, dks, lscf, Gn, gn);
#ifdef _OPENMP
#pragma omp parallel private(tG, s1, s2, s3, min_lscf)
{
#pragma omp for
#endif
        for(Index i1 = 1; i1 < k; i1++) {
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, 0, k - i1 - 1, m + 1) - lscf.ULCat(i1 - 1, 0, k - i1, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, 0, k - i1, m + 1) - lscf.ULCat(i1, 0, k - i1 - 1, m + 1)));
            tG = s1 * A1 * (dks.ULCat(i1 - 1, 0, k - i1, m + 1) + Go.col(i1 - 1)) +
                 s3 * A3 * (dks.ULCat(i1, 0, k - i1 - 1, m + 1) + Go.col(i1));
            gn.col(i1) =
                (tG - s1 * Go.col(i1 - 1) - s3 * Go.col(i1)
                 - ((s1 * dks.ULCat(i1 - 1, 0, k - i1, m + 1) + s3 * dks.ULCat(i1, 0, k - i1 - 1, m + 1)))) * mu +
                s1 * A1 * go.col(i1 - 1) + s3 * A3 * go.col(i1);
            Gn.col(i1) = tG;
            dks.ULCat(i1, 0, k - i1, m + 1) = (Gn.col(i1).sum() + (mu * gn.col(i1)).sum()) / (2 * k);
            scale_in_h3_ijk_vE(i1, 0, k, m, n, thr, dks, lscf, Gn, gn);
            for(Index i2 = 1; i2 < k - i1; i2++) {
                Index i3 = k - i1 - i2;
                min_lscf = min<Scalar>({lscf.ULCat(i1 - 1, i2, i3, m + 1),
                                        lscf.ULCat(i1, i2 - 1, i3, m + 1),
                                        lscf.ULCat(i1, i2, i3 - 1, m + 1)});
                s1 = exp(min_lscf - lscf.ULCat(i1 - 1, i2, i3, m + 1));
                s2 = exp(min_lscf - lscf.ULCat(i1, i2 - 1, i3, m + 1));
                s3 = exp(min_lscf - lscf.ULCat(i1, i2, i3 - 1, m + 1));
                tG = s1 * A1 * (dks.ULCat(i1 - 1, i2, i3, m + 1) + Go.col(id3(i1 - 1, i2, k - 1))) +
                     s2 * A2 * (dks.ULCat(i1, i2 - 1, i3, m + 1) + Go.col(id3(i1, i2 - 1, k - 1))) +
                     s3 * A3 * (dks.ULCat(i1, i2, i3 - 1, m + 1) + Go.col(id3(i1, i2, k - 1)));
                gn.col(id3(i1, i2, k)) =
                    (tG - s1 * Go.col(id3(i1 - 1, i2, k - 1)) -
                     s2 * Go.col(id3(i1, i2 - 1, k - 1)) - s3 * Go.col(id3(i1, i2, k - 1))
                     - ((s1 * dks.ULCat(i1 - 1, i2, i3, m + 1) + s2 * dks.ULCat(i1, i2 - 1, i3, m + 1) +
                        s3 * dks.ULCat(i1, i2, i3 - 1, m + 1)))) * mu +
                    s1 * A1 * go.col(id3(i1 - 1, i2, k - 1)) + s2 * A2 * go.col(id3(i1, i2 - 1, k - 1)) + s3 * A3 * go.col(id3(i1, i2, k - 1));
                Gn.col(id3(i1, i2, k)) = tG;
                dks.ULCat(i1, i2, i3, m + 1) = (Gn.col(id3(i1, i2, k)).sum() + (mu * gn.col(id3(i1, i2, k))).sum()) / (2 * k);
                scale_in_h3_ijk_vE(i1, i2, k, m, n, thr, dks, lscf, Gn, gn);
            }
            s1 = exp(min<Scalar>(0, lscf.ULCat(i1, k - i1 - 1, 0, m + 1) - lscf.ULCat(i1 - 1, k - i1, 0, m + 1)));
            s2 = exp(min<Scalar>(0, lscf.ULCat(i1 - 1, k - i1, 0, m + 1) - lscf.ULCat(i1, k - i1 - 1, 0, m + 1)));
            tG = s1 * A1 * (dks.ULCat(i1 - 1, k - i1, 0, m + 1) + Go.col(id3(i1 - 1, k - i1, k - 1))) +
                 s2 * A2 * (dks.ULCat(i1, k - i1 - 1, 0, m + 1) + Go.col(id3(i1, k - i1 - 1, k - 1)));
            gn.col(id3(i1, k - i1, k)) =
                (tG - s1 * Go.col(id3(i1 - 1, k - i1, k - 1)) -
                 s2 * Go.col(id3(i1, k - i1 - 1, k - 1))
                 - ((s1 * dks.ULCat(i1 - 1, k - i1, 0, m + 1) + s2 * dks.ULCat(i1, k - i1 - 1, 0, m + 1)))) * mu +
                s1 * A1 * go.col(id3(i1 - 1, k - i1, k - 1)) + s2 * A2 * go.col(id3(i1, k - i1 - 1, k - 1));
            Gn.col(id3(i1, k - i1, k)) = tG;
            dks.ULCat(i1, k - i1, 0, m + 1) = (Gn.col(id3(i1, k - i1, k)).sum() + (mu * gn.col(id3(i1, k - i1, k))).sum()) / (2 * k);
            scale_in_h3_ijk_vE(i1, k - i1, k, m, n, thr, dks, lscf, Gn, gn);
        }
#ifdef _OPENMP
}
#endif
        tG = A1 * (dks.ULCat(k - 1, 0, 0, m + 1) + Go.col(k - 1));
        gn.col(k) = (tG - Go.col(k - 1)
             - ((dks.ULCat(k - 1, 0, 0, m + 1)))) * mu + A1 * go.col(k - 1);
        Gn.col(k) = tG;
        dks.ULCat(k, 0, 0, m + 1) = (Gn.col(k).sum() + (mu * gn.col(k)).sum()) / (2 * k);
        scale_in_h3_ijk_vE(k, 0, k, m, n, thr, dks, lscf, Gn, gn);
    }
    return dks;
}
template ArrayXd h3_ijk_vEc(const ArrayBase<ArrayXd>& A1,
                            const ArrayBase<ArrayXd>& A2,
                            const ArrayBase<ArrayXd>& A3,
                            const ArrayBase<ArrayXd>& mu, const Index m,
                            ArrayXd& lscf, const double thr_margin,
                            int nthreads);

template <typename DerivedA, typename DerivedB, typename DerivedC, typename DerivedD>
inline void scale_in_hxxx3_pjk_mE(Eigen::Index j, Eigen::Index k,
                                  const Eigen::Index m, const Eigen::Index n,
                                  const Eigen::Index p,
                                  typename DerivedA::Scalar &thr,
                                  Eigen::ArrayBase<DerivedA> &dks,
                                  Eigen::ArrayBase<DerivedB> &lscf,
                                  Eigen::MatrixBase<DerivedC> &Gn,
                                  Eigen::MatrixBase<DerivedD> &gn) {
    if(Gn.block(0, j * n * (p + 1), n, n * (p + 1)).maxCoeff() > thr || gn.block(0, j * (p + 1), n, p + 1).maxCoeff() > thr) {
        dks.ULScol(k - j, j, m + 1) /= 1e10;
        Gn.block(0, j * n * (p + 1), n, n * (p + 1)) /= 1e10;
        gn.block(0, j * (p + 1), n, p + 1) /= 1e10;
        update_scale_2D(lscf, k - j, j, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hxxx3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads,
             const typename Derived::Scalar c) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Matrix<Scalar, Dynamic, Dynamic> MatrixXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 2) / 2);
    dks.ULSat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    MatrixXx tG(n, n);
    MatrixXx Go = MatrixXx::Zero(n, n * (p + 1) * m);
    MatrixXx Gn = MatrixXx::Zero(n, n * (p + 1) * (m + 1));
    MatrixXx go = MatrixXx::Zero(n, (p + 1) * m);
    MatrixXx gn = MatrixXx::Zero(n, (p + 1) * (m + 1));
    Gn.block(0, 0, n, n).diagonal().array() += dks.ULSat(0, 0, 0, m + 1);
    Scalar s2, s3;
    for(Index i = 1; i <= p; i++) {
        MatrixXx::Map(Gn.block(0, i * n, n, n).data(), n, n).noalias() =
            A1 * Gn.block(0, (i - 1) * n, n, n);
        MatrixXx::Map(gn.col(i).data(), n, 1).noalias() =
            Gn.block(0, i * n, n, n) * mu + A1 * gn.col(i - 1);
        dks.ULSat(i, 0, 0, m + 1) =
            (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * i);
        Gn.block(0, i * n, n, n).diagonal().array() += dks.ULSat(i, 0, 0, m + 1);
    }
    scale_in_hxxx3_pjk_mE(0, 0, m, n, p, thr, dks, lscf, Gn, gn);
    for(Index k = 1; k <= m; k++) {
        if(k % 100 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * n * (p + 1)) = Gn.block(0, 0, n, k * n * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG.noalias() = A2 * Go.block(0, 0, n, n);
        MatrixXx::Map(gn.col(0).data(), n, 1).noalias() =
            (tG + c * Go.block(0, 0, n, n)) * mu + A2 * go.col(0);
        Gn.block(0, 0, n, n) = tG;
        dks.ULSat(0, k, 0, m + 1) =
            (Gn.block(0, 0, n, n).trace() + gn.col(0).dot(mu)) / (2 * k);
        Gn.block(0, 0, n, n).diagonal().array() += dks.ULSat(0, k, 0, m + 1);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * Gn.block(0, (i - 1) * n, n, n) +
                           A2 * Go.block(0, i * n, n, n);
            MatrixXx::Map(gn.col(i).data(), n, 1).noalias() =
                (tG + c * Go.block(0, i * n, n, n)) * mu +
                A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.block(0, i * n, n, n) = tG;
            dks.ULSat(i, k, 0, m + 1) =
                (Gn.block(0, i * n, n, n).trace() + gn.col(i).dot(mu)) / (2 * (k + i));
            Gn.block(0, i * n, n, n).diagonal().array() += dks.ULSat(i, k, 0, m + 1);
        }
        scale_in_hxxx3_pjk_mE(0, k, m, n, p, thr, dks, lscf, Gn, gn);
#ifdef _OPENMP
#pragma omp parallel private(tG, s2, s3)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            s2 = exp(min<Scalar>(0, lscf.ULTat(k - j, j - 1, m + 1) - lscf.ULTat(k - j - 1, j, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULTat(k - j - 1, j, m + 1) - lscf.ULTat(k - j, j - 1, m + 1)));
            tG.noalias() = s2 * A2 * Go.block(0, j * n * (p + 1), n, n) +
                           s3 * A3 * Go.block(0, (j - 1) * n * (p + 1), n, n);
            MatrixXx::Map(gn.col(j * (p + 1)).data(), n, 1).noalias() =
                (tG + c * (s2 * Go.block(0, j * n * (p + 1), n, n) + s3 * Go.block(0, (j - 1) * n * (p + 1), n, n))) * mu +
                s2 * A2 * go.col(j * (p + 1)) + s3 * A3 * go.col((j - 1) * (p + 1));
            Gn.block(0, j * n * (p + 1), n, n) = tG;
            dks.ULSat(0, k - j, j, m + 1) =
                (Gn.block(0, j * n * (p + 1), n, n).trace() + gn.col(j * (p + 1)).dot(mu)) / (2 * k);
            Gn.block(0, j * n * (p + 1), n, n).diagonal().array() += dks.ULSat(0, k - j, j, m + 1);
            for(Index i = 1; i <= p; i++) {
                tG.noalias() =      A1 * Gn.block(0, j * n * (p + 1) + (i - 1) * n, n, n) +
                               s2 * A2 * Go.block(0, j * n * (p + 1) + i * n, n, n) +
                               s3 * A3 * Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n);
                MatrixXx::Map(gn.col(j * (p + 1) + i).data(), n, 1).noalias() =
                    (tG + c * (s2 * Go.block(0, j * n * (p + 1) + i * n, n, n) +
                               s3 * Go.block(0, (j - 1) * n * (p + 1) + i * n, n, n))) * mu +
                    A1 * gn.col(j * (p + 1) + i - 1) + s2 * A2 * go.col(j * (p + 1) + i) + s3 * A3 * go.col((j - 1) * (p + 1) + i);
                Gn.block(0, j * n * (p + 1) + i * n, n, n) = tG;
                dks.ULSat(i, k - j, j, m + 1) =
                    (Gn.block(0, j * n * (p + 1) + i * n, n, n).trace() + gn.col(j * (p + 1) + i).dot(mu)) / (2 * (k + i));
                Gn.block(0, j * n * (p + 1) + i * n, n, n).diagonal().array() += dks.ULSat(i, k - j, j, m + 1);
            }
            scale_in_hxxx3_pjk_mE(j, k, m, n, p, thr, dks, lscf, Gn, gn);
        }
#ifdef _OPENMP
}
#endif
        tG.noalias() = A3 * Go.block(0, (k - 1) * n * (p + 1), n, n);
        MatrixXx::Map(gn.col(k * (p + 1)).data(), n, 1).noalias() =
            (tG + c * Go.block(0, (k - 1) * n * (p + 1), n, n)) * mu + A3 * go.col((k - 1) * (p + 1));
        Gn.block(0, k * n * (p + 1), n, n) = tG;
        dks.ULSat(0, 0, k, m + 1) =
            (Gn.block(0, k * n * (p + 1), n, n).trace() + gn.col(k * (p + 1)).dot(mu)) / (2 * k);
        Gn.block(0, k * n * (p + 1), n, n).diagonal().array() += dks.ULSat(0, 0, k, m + 1);
        for(Index i = 1; i <= p; i++) {
            tG.noalias() = A1 * Gn.block(0, k * n * (p + 1) + (i - 1) * n, n, n) +
                           A3 * Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n);
            MatrixXx::Map(gn.col(k * (p + 1) + i).data(), n, 1).noalias() =
                (tG + c * Go.block(0, (k - 1) * n * (p + 1) + i * n, n, n)) * mu +
                A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.block(0, k * n * (p + 1) + i * n, n, n) = tG;
            dks.ULSat(i, 0, k, m + 1) =
                (Gn.block(0, k * n * (p + 1) + i * n, n, n).trace() + gn.col(k * (p + 1) + i).dot(mu)) / (2 * (k + i));
            Gn.block(0, k * n * (p + 1) + i * n, n, n).diagonal().array() += dks.ULSat(i, 0, k, m + 1);
        }
        scale_in_hxxx3_pjk_mE(k, k, m, n, p, thr, dks, lscf, Gn, gn);
    }
    return dks;
}

template <typename DerivedA, typename DerivedB, typename DerivedC>
inline void scale_in_hxxx3_pjk_vE(Eigen::Index j, Eigen::Index k,
                                  const Eigen::Index m, const Eigen::Index n,
                                  const Eigen::Index p,
                                  typename DerivedA::Scalar &thr,
                                  Eigen::ArrayBase<DerivedA> &dks,
                                  Eigen::ArrayBase<DerivedB> &lscf,
                                  Eigen::ArrayBase<DerivedC> &Gn,
                                  Eigen::ArrayBase<DerivedC> &gn) {
    if(Gn.block(0, j * (p + 1), n, p + 1).maxCoeff() > thr || gn.block(0, j * (p + 1), n, p + 1).maxCoeff() > thr) {
        dks.ULScol(k - j, j, m + 1) /= 1e10;
        Gn.block(0, j * (p + 1), n, p + 1) /= 1e10;
        gn.block(0, j * (p + 1), n, p + 1) /= 1e10;
        update_scale_2D(lscf, k - j, j, m + 1);
    }
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hxxx3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads,
             const typename Derived::Scalar c) {
#ifdef _OPENMP
    if(nthreads <= 0) nthreads = std::max(omp_get_num_procs() / 2, 1);
    omp_set_num_threads(nthreads);
#endif
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    const Index n = A1.rows();
    ArrayXXx dks = ArrayXXx::Zero(p + 1, (m + 1) * (m + 2) / 2);
    dks.ULSat(0, 0, 0, m + 1) = 1;
    Scalar thr = std::numeric_limits<Scalar>::max() / thr_margin / Scalar(n);
    ArrayXx tG(n);
    ArrayXXx Go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx Gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    ArrayXXx go = ArrayXXx::Zero(n, (p + 1) * m);
    ArrayXXx gn = ArrayXXx::Zero(n, (p + 1) * (m + 1));
    Scalar s2, s3;
    for(Index i = 1; i <= p; i++) {
        Gn.col(i) = A1 * (dks.ULSat(i - 1, 0, 0, m + 1) + Gn.col(i - 1));
        gn.col(i) = Gn.col(i) * mu + A1 * gn.col(i - 1);
        dks.ULSat(i, 0, 0, m + 1) =
            (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * i);
    }
    scale_in_hxxx3_pjk_vE(0, 0, m, n, p, thr, dks, lscf, Gn, gn);
    for(Index k = 1; k <= m; k++) {
        if(k % 500 == 0) {
            Rcpp::checkUserInterrupt();
        }
        Go.block(0, 0, n, k * (p + 1)) = Gn.block(0, 0, n, k * (p + 1));
        go.block(0, 0, n, k * (p + 1)) = gn.block(0, 0, n, k * (p + 1));
        tG = A2 * (dks.ULSat(0, k - 1, 0, m + 1) + Go.col(0));
        gn.col(0) =
            (tG + c * Go.col(0) + c * dks.ULSat(0, k - 1, 0, m + 1)) * mu +
            A2 * go.col(0);
        Gn.col(0) = tG;
        dks.ULSat(0, k, 0, m + 1) =
            (Gn.col(0).sum() + (mu * gn.col(0)).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG = A1 * (dks.ULSat(i - 1, k, 0, m + 1) + Gn.col(i - 1)) +
                 A2 * (dks.ULSat(i, k - 1, 0, m + 1) + Go.col(i));
            gn.col(i) = (tG + c * Go.col(i)
                         + c * dks.ULSat(i, k - 1, 0, m + 1)) * mu +
                        A1 * gn.col(i - 1) + A2 * go.col(i);
            Gn.col(i) = tG;
            dks.ULSat(i, k, 0, m + 1) =
                (Gn.col(i).sum() + (mu * gn.col(i)).sum()) / (2 * (k + i));
        }
        scale_in_hxxx3_pjk_vE(0, k, m, n, p, thr, dks, lscf, Gn, gn);
#ifdef _OPENMP
#pragma omp parallel private(tG, s2, s3)
{
#pragma omp for
#endif
        for(Index j = 1; j < k; j++) {
            s2 = exp(min<Scalar>(0, lscf.ULTat(k - j, j - 1, m + 1) - lscf.ULTat(k - j - 1, j, m + 1)));
            s3 = exp(min<Scalar>(0, lscf.ULTat(k - j - 1, j, m + 1) - lscf.ULTat(k - j, j - 1, m + 1)));
            tG = s2 * A2 * (dks.ULSat(0, k - j - 1, j, m + 1) + Go.col(j * (p + 1))) +
                 s3 * A3 * (dks.ULSat(0, k - j, j - 1, m + 1) + Go.col((j - 1) * (p + 1)));
            gn.col(j * (p + 1)) =
                (tG + c * (s2 * Go.col(j * (p + 1)) + s3 * Go.col((j - 1) * (p + 1)))
                    + c * (s2 * dks.ULSat(0, k - j - 1, j, m + 1) + s3 * dks.ULSat(0, k - j, j - 1, m + 1))) * mu +
                s2 * A2 * go.col(j * (p + 1)) + s3 * A3 * go.col((j - 1) * (p + 1));
            Gn.col(j * (p + 1)) = tG;
            dks.ULSat(0, k - j, j, m + 1) =
                (Gn.col(j * (p + 1)).sum() + (mu * gn.col(j * (p + 1))).sum()) / (2 * k);
            for(Index i = 1; i <= p; i++) {
                tG =      A1 * (dks.ULSat(i - 1, k - j, j, m + 1) + Gn.col(j * (p + 1) + i - 1)) +
                     s2 * A2 * (dks.ULSat(i, k - j - 1, j, m + 1) + Go.col(j * (p + 1) + i)) +
                     s3 * A3 * (dks.ULSat(i, k - j, j - 1, m + 1) + Go.col((j - 1) * (p + 1) + i));
                gn.col(j * (p + 1) + i) =
                    (tG + c * (s2 * Go.col(j * (p + 1) + i) + s3 * Go.col((j - 1) * (p + 1) + i))
                        + c * (s2 * dks.ULSat(i, k - j - 1, j, m + 1) + s3 * dks.ULSat(i, k - j, j - 1, m + 1))) * mu +
                    A1 * gn.col(j * (p + 1) + i - 1) + s2 * A2 * go.col(j * (p + 1) + i) + s3 * A3 * go.col((j - 1) * (p + 1) + i);
                Gn.col(j * (p + 1) + i) = tG;
                dks.ULSat(i, k - j, j, m + 1) =
                    (Gn.col(j * (p + 1) + i).sum() + (mu * gn.col(j * (p + 1) + i)).sum()) / (2 * (k + i));
            }
            scale_in_hxxx3_pjk_vE(j, k, m, n, p, thr, dks, lscf, Gn, gn);
        }
#ifdef _OPENMP
}
#endif
        tG = A3 * (dks.ULSat(0, 0, k - 1, m + 1) + Go.col((k - 1) * (p + 1)));
        gn.col(k * (p + 1)) =
            (tG + c * Go.col((k - 1) * (p + 1))
                + c * dks.ULSat(0, 0, k - 1, m + 1)) * mu +
            A3 * go.col((k - 1) * (p + 1));
        Gn.col(k * (p + 1)) = tG;
        dks.ULSat(0, 0, k, m + 1) =
            (Gn.col(k * (p + 1)).sum() + (mu * gn.col(k * (p + 1))).sum()) / (2 * k);
        for(Index i = 1; i <= p; i++) {
            tG = A1 * (dks.ULSat(i - 1, 0, k, m + 1) + Gn.col(k * (p + 1) + i - 1)) +
                 A3 * (dks.ULSat(i, 0, k - 1, m + 1) + Go.col((k - 1) * (p + 1) + i));
            gn.col(k * (p + 1) + i) =
                (tG + c * Go.col((k - 1) * (p + 1) + i)
                    + c * dks.ULSat(i, 0, k - 1, m + 1)) * mu +
                A1 * gn.col(k * (p + 1) + i - 1) + A3 * go.col((k - 1) * (p + 1) + i);
            Gn.col(k * (p + 1) + i) = tG;
            dks.ULSat(i, 0, k, m + 1) = (Gn.col(k * (p + 1) + i).sum() + (mu * gn.col(k * (p + 1) + i)).sum()) / (2 * (k + i));
        }
        scale_in_hxxx3_pjk_vE(k, k, m, n, p, thr, dks, lscf, Gn, gn);
    }
    return dks;
}

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
    return hxxx3_pjk_mEc(A1, A2, A3, mu, m, p, lscf, thr_margin, nthreads, -1);
}
template ArrayXXd htil3_pjk_mEc(const MatrixBase<MatrixXd>& A1,
                               const DiagMatXd& A2,
                               const MatrixBase<MatrixXd>& A3,
                               const VectorXd mu, const Index m, const Index p,
                               ArrayXd& lscf, const double thr_margin,
                               int nthreads);


// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
    return hxxx3_pjk_vEc(A1, A2, A3, mu, m, p, lscf, thr_margin, nthreads, -1);
}
template ArrayXXd htil3_pjk_vEc(const ArrayBase<ArrayXd>& A1,
                               const ArrayBase<ArrayXd>& A2,
                               const ArrayBase<ArrayXd>& A3,
                               const ArrayBase<ArrayXd>& mu,
                               const Index m, const Index p, ArrayXd& lscf,
                               const double thr_margin, int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
    return hxxx3_pjk_mEc(A1, A2, A3, mu, m, p, lscf, thr_margin, nthreads, 1);
}
template ArrayXXd hhat3_pjk_mEc(const MatrixBase<MatrixXd>& A1,
                               const DiagMatXd& A2,
                               const MatrixBase<MatrixXd>& A3,
                               const VectorXd mu, const Index m, const Index p,
                               ArrayXd& lscf, const double thr_margin,
                               int nthreads);

// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads) {
    return hxxx3_pjk_vEc(A1, A2, A3, mu, m, p, lscf, thr_margin, nthreads, 1);
}
template ArrayXXd hhat3_pjk_vEc(const ArrayBase<ArrayXd>& A1,
                               const ArrayBase<ArrayXd>& A2,
                               const ArrayBase<ArrayXd>& A3,
                               const ArrayBase<ArrayXd>& mu,
                               const Index m, const Index p, ArrayXd& lscf,
                               const double thr_margin, int nthreads);
