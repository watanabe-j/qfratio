#ifndef DK_FUNS_CWISE_H
#define DK_FUNS_CWISE_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d2_ij_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2,
         const int m, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_mE(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::MatrixBase<Derived>& A2,
         const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h2_ij_vE(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::ArrayBase<Derived>& mu,
         const int m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf,
          int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mE(const Eigen::MatrixBase<Derived>& A1, const Eigen::MatrixBase<Derived>& A2, const Eigen::MatrixBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vE(const Eigen::ArrayBase<Derived>& A1, const Eigen::ArrayBase<Derived>& A2, const Eigen::ArrayBase<Derived>& A3,
          const int m, const int p, Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_mE(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::MatrixBase<Derived>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
          const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf,
          int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
h3_ijk_vE(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::ArrayBase<Derived>& mu, const int m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf,
             int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf); // , int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_mE(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::MatrixBase<Derived>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf,
             int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_vE(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const int m, const int p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>& lscf); // , int nthreads);


#endif
