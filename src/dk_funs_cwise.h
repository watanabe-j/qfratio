#ifndef DK_FUNS_CWISE_H
#define DK_FUNS_CWISE_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d2_ij_mEc(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d2_ij_vEc(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h2_ij_mEc(const Eigen::MatrixBase<Derived>& A1,
         const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
         const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1>& mu,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h2_ij_vEc(const Eigen::ArrayBase<Derived>& A1,
         const Eigen::ArrayBase<Derived>& A2,
         const Eigen::ArrayBase<Derived>& mu,
         const Eigen::Index m,
         Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
         const typename Derived::Scalar thr_margin, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d3_ijk_mEc(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
d3_ijk_vEc(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Index m, const Eigen::Index p,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
d3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::Index m, const Eigen::Index p,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h3_ijk_mEc(const Eigen::MatrixBase<Derived>& A1,
          const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
          const Eigen::MatrixBase<Derived>& A3,
          const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
          const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
h3_ijk_vEc(const Eigen::ArrayBase<Derived>& A1,
          const Eigen::ArrayBase<Derived>& A2,
          const Eigen::ArrayBase<Derived>& A3,
          const Eigen::ArrayBase<Derived>& mu, const Eigen::Index m,
          Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
          const typename Derived::Scalar thr_margin, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
htil3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_mEc(const Eigen::MatrixBase<Derived>& A1,
             const Eigen::DiagonalMatrix<typename Derived::Scalar, Eigen::Dynamic>& A2,
             const Eigen::MatrixBase<Derived>& A3,
             const Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads);
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hhat3_pjk_vEc(const Eigen::ArrayBase<Derived>& A1,
             const Eigen::ArrayBase<Derived>& A2,
             const Eigen::ArrayBase<Derived>& A3,
             const Eigen::ArrayBase<Derived>& mu,
             const Eigen::Index m, const Eigen::Index p,
             Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf,
             const typename Derived::Scalar thr_margin, int nthreads);


#endif
