#ifndef DK_FUNS_H
#define DK_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

// #include "qfratio_types.h"

Eigen::ArrayXd d1_i_vE(const Eigen::ArrayXd& L, const int m, Eigen::ArrayXd& lscf);
Eigen::ArrayXd d1_i_mE(const Eigen::MatrixXd& A, const int m, Eigen::ArrayXd& lscf);
Eigen::ArrayXd dtil1_i_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& mu,
    const int m, Eigen::ArrayXd& lscf);
Eigen::ArrayXd dtil1_i_mE(const Eigen::MatrixXd& A, const Eigen::VectorXd& mu,
    const int m, Eigen::ArrayXd& lscf);
Eigen::ArrayXXd arl_vE(const Eigen::ArrayXd& L, const Eigen::ArrayXd& D, const int m);
Eigen::ArrayXXd arl_mE(const Eigen::MatrixXd& A, const Eigen::VectorXd& mu, const int m);

Eigen::ArrayXXd d2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const int m, const int p, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd d2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const int m, const int p, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd d2_ij_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const int m, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd d2_ij_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const int m, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd h2_ij_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::VectorXd& mu, const int m, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd h2_ij_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int m, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd htil2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
         const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd htil2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd hhat2_pj_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::VectorXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd hhat2_pj_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int m, const int p, Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd dtil2_pq_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::VectorXd& mu, const int p, const int q);
Eigen::ArrayXXd dtil2_pq_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& mu, const int p, const int q);

Eigen::ArrayXXd d3_ijk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::MatrixXd& A3, const int m, Eigen::ArrayXXd& lscf, int nthreads);
Eigen::ArrayXXd d3_ijk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& A3, const int m, Eigen::ArrayXXd& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mE(const Eigen::MatrixBase<Derived>& A1, const Eigen::MatrixBase<Derived>& A2, const Eigen::MatrixBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2, const Eigen::ArrayBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf);

Eigen::ArrayXXd h3_ijk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::MatrixXd& A3, const Eigen::VectorXd& mu, const int m,
    Eigen::ArrayXXd& lscf, int nthreads);
Eigen::ArrayXXd h3_ijk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& A3, const Eigen::ArrayXd& mu, const int m, Eigen::ArrayXXd& lscf); // , int nthreads);
Eigen::ArrayXXd htil3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::MatrixXd& A3, const Eigen::VectorXd& mu, const int m,
    const int p, Eigen::ArrayXXd& lscf, int nthreads);
Eigen::ArrayXXd htil3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& A3, const Eigen::ArrayXd& mu, const int m,
    const int p, Eigen::ArrayXXd& lscf); // , int nthreads);
Eigen::ArrayXXd hhat3_pjk_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::MatrixXd& A3, const Eigen::VectorXd& mu, const int m,
    const int p, Eigen::ArrayXXd& lscf, int nthreads);
Eigen::ArrayXXd hhat3_pjk_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& A3, const Eigen::ArrayXd& mu, const int m,
    const int p, Eigen::ArrayXXd& lscf); // , int nthreads);
Eigen::ArrayXXd dtil3_pqr_mE(const Eigen::MatrixXd& A1, const Eigen::MatrixXd& A2,
    const Eigen::MatrixXd& A3, const Eigen::VectorXd mu,
    const int p, const int q, const int r);
Eigen::ArrayXXd dtil3_pqr_vE(const Eigen::ArrayXd& A1, const Eigen::ArrayXd& A2,
    const Eigen::ArrayXd& A3, const Eigen::ArrayXd& mu,
    const int p, const int q, const int r);

#endif
