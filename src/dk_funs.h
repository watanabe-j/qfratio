#ifndef DK_FUNS_H
#define DK_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// #include "qfratio_types.h"

Eigen::ArrayXd d1_i_vE(const Eigen::ArrayXd& L, const int m, Eigen::ArrayXd& lscf);
Eigen::ArrayXd d1_i_mE(const Eigen::MatrixXd& A, const int m, Eigen::ArrayXd& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
dtil1_i_vE(const Eigen::ArrayBase<Derived>& L,
           const Eigen::ArrayBase<Derived>& mu, const int m,
           Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
dtil1_i_mE(const Eigen::MatrixBase<Derived>& A,
           const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
           const int m,
           Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);

Eigen::ArrayXXd arl_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& D, const int m);
Eigen::ArrayXXd arl_mE(const Eigen::MatrixXd& A, const Eigen::VectorXd& mu, const int m);

Eigen::ArrayXXd d2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const int m, const int p, Eigen::ArrayXd& lscf);
Eigen::ArrayXXd d2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const int m, const int p, Eigen::ArrayXd& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2,
         const int m, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_vE(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::ArrayBase<Derived>& mu,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);

Eigen::ArrayXXd htil2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXd& lscf);
Eigen::ArrayXXd htil2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXd& lscf);
Eigen::ArrayXXd hhat2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXd& lscf);
Eigen::ArrayXXd hhat2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXd& lscf);
Eigen::ArrayXXd dtil2_pq_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::VectorXd& mu, const int p, const int q);
Eigen::ArrayXXd dtil2_pq_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int p, const int q);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mE(const Eigen::MatrixBase<Derived>& A1, const Eigen::MatrixBase<Derived>& A2, const Eigen::MatrixBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2, const Eigen::ArrayBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::ArrayBase<Derived>& mu, const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf); // , int nthreads);

Eigen::ArrayXXd dtil3_pqr_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::MatrixXd& A3, const Eigen::VectorXd mu,
    const int p, const int q, const int r);
Eigen::ArrayXXd dtil3_pqr_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& A3, const Eigen::ArrayXd& mu,
    const int p, const int q, const int r);

#endif
