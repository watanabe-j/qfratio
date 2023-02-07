#ifndef HGS_FUNS_H
#define HGS_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// #include "qfratio_types.h"

// template <typename T>
// Eigen::Array<T, Eigen::Dynamic, 1> get_lrf(const T a, const Eigen::Index n);
Eigen::ArrayXd get_lrf(const double a, const Eigen::Index n);
Eigen::Array<long double, Eigen::Dynamic, 1> get_lrf(const long double a, const Eigen::Index n);

template <typename Derived>
void set_cumsum(const Eigen::DenseBase<Derived>& Data, Eigen::DenseBase<Derived>& Out);

template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1>  get_sign_rf(const T a, const Eigen::Index n);

template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1> get_sign_rfp1(const T a, const Eigen::Index n);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_1dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1, 
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_2dEc(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf);


template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_3dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1,
        const typename Derived::Scalar a2, const typename Derived::Scalar a3,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_3dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar a3, const typename Derived::Scalar b,
        const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf);


template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiagE(const Eigen::ArrayBase<Derived>& X);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiagEc(const Eigen::ArrayBase<Derived>& X);

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiag3DE(const Eigen::ArrayBase<Derived>& X);

#endif
