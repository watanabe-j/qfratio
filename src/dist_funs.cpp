#include "config.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

// // These are to use gsl
// #include <RcppGSL.h>
// // [[Rcpp::depends(RcppGSL)]]
// #include <gsl/gsl_sf_hyperg.h>

#include "dk_funs.h"
#include "hgs_funs.h"

using Eigen::log;
using Eigen::abs;
using Eigen::ArrayXd;
using Eigen::Index;


Eigen::ArrayXd get_lgm(const double a, const Eigen::Index n) {
    return ArrayXd::LinSpaced(n, a, a + n - 1).lgamma();
}

//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}
//'
// [[Rcpp::export]]
SEXP A1B1_E(const Eigen::ArrayXd D1, const Eigen::ArrayXd D2,
            const Eigen::Index m, const double thr_margin = 100) {
    Index n1 = D1.size();
    Index n2 = D2.size();
    double n1_ = n1;
    double n2_ = n2;
    double m_ = m;
    Index size_out = (m + 1) * (m + 2) / 2;
    ArrayXd D1i = D1.cwiseInverse();
    ArrayXd D2i = D2.cwiseInverse();
    double sumtrDi = D1i.sum() + D2i.sum();
    ArrayXd D1d = D1i / sumtrDi;
    ArrayXd D2d = D2i / sumtrDi;
    double trD1d = D1d.sum();
    double trD2d = D2d.sum();
    ArrayXd D1h = ArrayXd::Constant(n1, trD1d) - D1d;
    ArrayXd D2h = ArrayXd::Constant(n2, trD2d) - D2d;
    ArrayXd lscf1 = ArrayXd::Zero(m + 1);
    ArrayXd lscf2 = ArrayXd::Zero(m + 1);
    ArrayXd dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
    ArrayXd dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
    ArrayXd ansmat = ArrayXd::Zero(size_out);
    ArrayXd Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXd Aldenj = get_lgm(n2_ / 2, m + 1);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULTcol(k, m + 1) +=
            // log(abs(dk1.head(m + 1 - k))) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
            // log(abs(dk2(k))) - Aldenj(k) - lscf2(k);
            log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
            log(dk2(k)) - Aldenj(k) - lscf2(k);
        double lnum = std::lgamma((n1_ + n2_) / 2 + k);
        for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) += lnum;
    }
    ansmat += (log(D1d).sum() + log(D2d).sum()) / 2;
    ansmat = exp(ansmat);
    // for(Index k = 0; k <= m; k++) {
    //     ansmat.ULTcol(k, m + 1) *= sign(dk1.head(m + 1 - k)) * signd(dk2(k));
    // }

    ArrayXd ls = ArrayXd::LinSpaced(m + 1, 0, m_);
    ArrayXd a1s(size_out);
    ArrayXd bs(size_out);
    for(Index k = 0; k <= m; k++) {
        double k_ = k;
        a1s.ULTcol(k, m + 1) = ls.head(m + 1 - k) + k_ + (n1_ + n2_) / 2;
        bs.ULTcol(k, m + 1) = ls.head(m + 1 - k) + n1_ / 2 + 1;
    }

    // // Using the C++ library GSL with RcppGSL; this is
    // // not quite portable as the library need to be installed separately.
    // ArrayXd hgres(size_out);
    // for(Index i = 0; i < size_out; i++) {
    //     hgres(i) = gsl_sf_hyperg_2F1(a1s(i), 1, bs(i), trD1d);
    // }

    // Calling gsl::hyperg_2F1 from R
    Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
    Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
    Eigen::Map<ArrayXd> hgres(hgv.begin(), size_out);

    ansmat *= hgres;

    ArrayXd ansseq = sum_counterdiagE(ansmat);
    return Rcpp::List::create(
        Rcpp::Named("ansseq") = ansseq);
}
