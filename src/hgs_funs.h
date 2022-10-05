#ifndef HGS_FUNS_H
#define HGS_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// #include "qfratio_types.h"

// template <typename T>
// Eigen::Array<T, Eigen::Dynamic, 1> get_lrf(const T a, const int n);
Eigen::ArrayXd get_lrf(const double a, const int n);
Eigen::Array<long double, Eigen::Dynamic, 1> get_lrf(const long double a, const int n);

template <typename Derived>
void set_cumsum(const Eigen::DenseBase<Derived>& Data, Eigen::DenseBase<Derived>& Out);

template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1>  get_sign_rf(const T a, const int n);

template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1> get_sign_rfp1(const T a, const int n);

Eigen::ArrayXd hgs_1dE(const Eigen::ArrayXd& dks,
    const double a1,
    const double b, const double lconst, const Eigen::ArrayXd& lscf);

// template <typename Derived>
// Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
// hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
//         const typename Derived::Scalar a1, const typename Derived::Scalar a2,
//         const double b, const typename Derived::Scalar lconst,
//         const Eigen::ArrayBase<Derived>& lscf);

Eigen::ArrayXXd
hgs_2dE(const Eigen::ArrayXXd& dks,
        const double a1, const double a2, const double b,
        // const double lconst, const Eigen::ArrayXXd& lscf);
        const double lconst, const Eigen::ArrayXd& lscf);
Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>
hgs_2dE(const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& dks,
        const long double a1, const long double a2, const long double b,
        const long double lconst,
        // const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& lscf);
        const Eigen::Array<long double, Eigen::Dynamic, 1>& lscf);

// template <typename Derived>
// Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
// hgs_3dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1,
//         const typename Derived::Scalar a2, const typename Derived::Scalar a3,
//         const double b, const typename Derived::Scalar lconst,
//         const Eigen::ArrayBase<Derived>& lscf);

Eigen::ArrayXXd
hgs_3dE(const Eigen::ArrayXXd& dks,
        const double a1, const double a2, const double a3,
        // const double b, const double lconst, const Eigen::ArrayXXd& lscf);
        const double b, const double lconst, const Eigen::ArrayXd& lscf);
Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>
hgs_3dE(const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& dks,
        const long double a1, const long double a2, const long double a3,
        const long double b, const long double lconst,
        // const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& lscf);
        const Eigen::Array<long double, Eigen::Dynamic, 1>& lscf);
#endif
