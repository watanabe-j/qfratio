#ifndef DK_FUNS_H
#define DK_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// #include "qfratio_types.h"

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d1_i_vE(const Eigen::ArrayBase<Derived>& L, const int m, 
        Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
        const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d1_i_mE(const Eigen::MatrixBase<Derived>& A, const int m,
        Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
        const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
dtil1_i_vE(const Eigen::ArrayBase<Derived>& L,
           const Eigen::ArrayBase<Derived>& mu, const int m,
           Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
           const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
dtil1_i_mE(const Eigen::MatrixBase<Derived> &A,
           const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> &mu,
           const int m,
           Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> &lscf,
           const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
arl_vE(const Eigen::ArrayBase<Derived>& L, const Eigen::ArrayBase<Derived>& D,
       const int m, const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
arl_mE(const Eigen::MatrixBase<Derived>& A,
       const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
       const int m, const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_pj_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2, const int m, const int p,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_pj_vE(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2, const int m, const int p,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2,
         const int m, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_vE(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::ArrayBase<Derived>& mu,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil2_pj_mE(const Eigen::MatrixBase<Derived>& A1,
            const Eigen::MatrixBase<Derived>& A2,
            const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
            const int m, const int p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil2_pj_vE(const Eigen::ArrayBase<Derived>& A1,
            const Eigen::ArrayBase<Derived>& A2,
            const Eigen::ArrayBase<Derived>& mu, const int m, const int p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat2_pj_mE(const Eigen::MatrixBase<Derived>& A1,
            const Eigen::MatrixBase<Derived>& A2,
            const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
            const int m, const int p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat2_pj_vE(const Eigen::ArrayBase<Derived>& A1,
            const Eigen::ArrayBase<Derived>& A2,
            const Eigen::ArrayBase<Derived>& mu,
            const int m, const int p,
            Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
            const typename Derived::Scalar thr_margin);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil2_pq_mE(const Eigen::MatrixBase<Derived>& A1,
            const Eigen::MatrixBase<Derived>& A2,
            const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
            const int p, const int q);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil2_pq_vE(const Eigen::ArrayBase<Derived>& A1,
            const Eigen::ArrayBase<Derived>& A2,
            const Eigen::ArrayBase<Derived>& mu, const int p, const int q);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> &lscf,
          const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mE(const Eigen::MatrixBase<Derived>& A1, const Eigen::MatrixBase<Derived>& A2, const Eigen::MatrixBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2, const Eigen::ArrayBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::ArrayBase<Derived>& mu, const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil3_pqr_mE(const Eigen::MatrixBase<Derived> &A1,
             const Eigen::MatrixBase<Derived> &A2,
             const Eigen::MatrixBase<Derived> &A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int p, const int q, const int r);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
dtil3_pqr_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int p, const int q, const int r);

#endif
