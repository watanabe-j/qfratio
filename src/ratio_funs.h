#ifndef RATIO_FUNS_H
#define RATIO_FUNS_H

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

Eigen::ArrayXd sum_counterdiagE(const Eigen::ArrayXXd & X);

Eigen::ArrayXd sum_counterdiag3DE(const Eigen::ArrayXXd & X);

#endif
