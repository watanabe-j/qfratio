#ifndef HGS_FUNS_H
#define HGS_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
// #include "qfratio_types.h"

Eigen::ArrayXd get_lrf(const double a, const int n);
void set_cumsum(const Eigen::ArrayXd& Data, Eigen::ArrayXd& Out);
// void set_cumprod_sign(const double a, const int n, Eigen::ArrayXd& Out);
Eigen::ArrayXd get_sign_rf(const double a, const int n);
Eigen::ArrayXd get_sign_rfp1(const double a, const int n);

Eigen::ArrayXd hgs_1dE(const Eigen::ArrayXd& dks,
    const double a1,
    const double b, const double lconst, const Eigen::ArrayXd& lscf);
Eigen::ArrayXXd hgs_2dE(const Eigen::ArrayXXd& dks,
    const double a1, const double a2,
    const double b, const double lconst, const Eigen::ArrayXXd& lscf);
Eigen::ArrayXXd hgs_3dE(const Eigen::ArrayXXd& dks,
    const double a1, const double a2, const double a3,
    const double b, const double lconst, const Eigen::ArrayXXd& lscf);

#endif
