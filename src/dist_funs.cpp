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
#include "dk_funs_cwise.h"

using Eigen::log;
using Eigen::abs;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Index;

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatXd;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, 1> VectorXl;
typedef Eigen::DiagonalMatrix<long double, Eigen::Dynamic> DiagMatXl;


Eigen::ArrayXd get_lgm(const double a, const Eigen::Index n) {
    return ArrayXd::LinSpaced(n, a, a + n - 1).lgamma();
}

ArrayXl get_lgm(const long double a, const Eigen::Index n) {
    ArrayXl ans(n);
    for(Index i = 0; i < n; i++) ans[i] = std::lgammal(a + i);
    return ans;
}


//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}, double
//'
// [[Rcpp::export]]
SEXP A1B1_Ed(const Eigen::ArrayXd D1, const Eigen::ArrayXd D2,
             const Eigen::ArrayXd mu1, const Eigen::ArrayXd mu2,
             const Eigen::Index m, const double thr_margin = 100,
             int nthreads = 0, const double tol_zero = 2.2e-14) {
    Index n1 = D1.size();
    Index n2 = D2.size();
    double n1_ = n1;
    double n2_ = n2;
    double m_ = m;
    bool central = is_zero_E(mu1, tol_zero) && is_zero_E(mu2, tol_zero);
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
    ArrayXd dk1(m + 1);
    ArrayXd dk2(m + 1);
    ArrayXd ansmat = ArrayXd::Zero(size_out);
    bool diminished = false;
    ArrayXd Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXd Aldenj = get_lgm(n2_ / 2, m + 1);
    if(central) {
        dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
        dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
    } else {
        VectorXd mu1d = D1d.sqrt() * mu1;
        VectorXd mu2d = D2d.sqrt() * mu2;
        MatrixXd mu1mat = mu1d * mu1d.transpose() / 2;
        MatrixXd mu2mat = mu2d * mu2d.transpose() / 2;
        DiagMatXd D1hmat = D1h.matrix().asDiagonal();
        DiagMatXd D2hmat = D2h.matrix().asDiagonal();
        ArrayXd dkm1 = d2_ij_mE(mu1mat, D1hmat, m, lscf1, thr_margin, nthreads);
        ArrayXd dkm2 = d2_ij_mE(mu2mat, D2hmat, m, lscf2, thr_margin, nthreads);
        dkm1 = log(dkm1);
        dkm2 = log(dkm2);
        ArrayXd seqlrf_1_2 = get_lrf(0.5, m + 1);
        for(Index k = 0; k <= m; k++) {
            dkm1.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            dkm2.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
        }
        dkm1 = exp(dkm1);
        dkm2 = exp(dkm2);
        dk1 = sum_counterdiagE(dkm1);
        dk2 = sum_counterdiagE(dkm2);
        diminished = ((lscf1 < 0).any() && dkm1.cwiseEqual(0).any()) ||
                     ((lscf2 < 0).any() && dkm2.cwiseEqual(0).any());
    }
    for(Index k = 0; k <= m; k++) {
        double k_ = k;
        ansmat.ULTcol(k, m + 1) +=
            log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
            log(dk2(k)) - Aldenj(k) - lscf2(k);
        double lnum = std::lgamma((n1_ + n2_) / 2 + k_);
        for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) += lnum;
    }
    ansmat += (log(D1d).sum() + log(D2d).sum()) / 2;
    ansmat -= (mu1.matrix().squaredNorm() + mu2.matrix().squaredNorm()) / 2;
    ansmat = exp(ansmat);

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
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}, long double
//'
// [[Rcpp::export]]
SEXP A1B1_El(const Eigen::Array<long double, Eigen::Dynamic, 1> D1,
             const Eigen::Array<long double, Eigen::Dynamic, 1> D2,
             const Eigen::Array<long double, Eigen::Dynamic, 1> mu1,
             const Eigen::Array<long double, Eigen::Dynamic, 1> mu2,
             const Eigen::Index m, const long double thr_margin = 100,
             int nthreads = 0, const long double tol_zero = 2.2e-14) {
    Index n1 = D1.size();
    Index n2 = D2.size();
    long double n1_ = n1;
    long double n2_ = n2;
    long double m_ = m;
    bool central = is_zero_E(mu1, tol_zero) && is_zero_E(mu2, tol_zero);
    Index size_out = (m + 1) * (m + 2) / 2;
    ArrayXl D1i = D1.cwiseInverse();
    ArrayXl D2i = D2.cwiseInverse();
    long double sumtrDi = D1i.sum() + D2i.sum();
    ArrayXl D1d = D1i / sumtrDi;
    ArrayXl D2d = D2i / sumtrDi;
    long double trD1d = D1d.sum();
    long double trD2d = D2d.sum();
    ArrayXl D1h = ArrayXl::Constant(n1, trD1d) - D1d;
    ArrayXl D2h = ArrayXl::Constant(n2, trD2d) - D2d;
    ArrayXl lscf1 = ArrayXl::Zero(m + 1);
    ArrayXl lscf2 = ArrayXl::Zero(m + 1);
    ArrayXl dk1(m + 1);
    ArrayXl dk2(m + 1);
    ArrayXl ansmat = ArrayXl::Zero(size_out);
    bool diminished = false;
    ArrayXl Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXl Aldenj = get_lgm(n2_ / 2, m + 1);
    if(central) {
        dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
        dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
    } else {
        VectorXl mu1d = D1d.sqrt() * mu1;
        VectorXl mu2d = D2d.sqrt() * mu2;
        MatrixXl mu1mat = mu1d * mu1d.transpose() / 2;
        MatrixXl mu2mat = mu2d * mu2d.transpose() / 2;
        DiagMatXl D1hmat = D1h.matrix().asDiagonal();
        DiagMatXl D2hmat = D2h.matrix().asDiagonal();
        ArrayXl dkm1 = d2_ij_mE(mu1mat, D1hmat, m, lscf1, thr_margin, nthreads);
        ArrayXl dkm2 = d2_ij_mE(mu2mat, D2hmat, m, lscf2, thr_margin, nthreads);
        dkm1 = log(dkm1);
        dkm2 = log(dkm2);
        ArrayXl seqlrf_1_2 = get_lrf((long double)(0.5), m + 1);
        for(Index k = 0; k <= m; k++) {
            dkm1.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            dkm2.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
        }
        dkm1 = exp(dkm1);
        dkm2 = exp(dkm2);
        dk1 = sum_counterdiagE(dkm1);
        dk2 = sum_counterdiagE(dkm2);
        diminished = ((lscf1 < 0).any() && dkm1.cwiseEqual(0).any()) ||
                     ((lscf2 < 0).any() && dkm2.cwiseEqual(0).any());
    }
    for(Index k = 0; k <= m; k++) {
        long double k_ = k;
        ansmat.ULTcol(k, m + 1) +=
            log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
            log(dk2(k)) - Aldenj(k) - lscf2(k);
        long double lnum = std::lgamma((n1_ + n2_) / 2 + k_);
        for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) += lnum;
    }
    ansmat += (log(D1d).sum() + log(D2d).sum()) / 2;
    ansmat -= (mu1.matrix().squaredNorm() + mu2.matrix().squaredNorm()) / 2;
    ansmat = exp(ansmat);

    ArrayXl ls = ArrayXl::LinSpaced(m + 1, 0, m_);
    ArrayXl a1s(size_out);
    ArrayXl bs(size_out);
    for(Index k = 0; k <= m; k++) {
        long double k_ = k;
        a1s.ULTcol(k, m + 1) = ls.head(m + 1 - k) + k_ + (n1_ + n2_) / 2;
        bs.ULTcol(k, m + 1) = ls.head(m + 1 - k) + n1_ / 2 + 1;
    }

    // // Using the C++ library GSL with RcppGSL; this is
    // // not quite portable as the library need to be installed separately.
    // ArrayXl hgres(size_out);
    // for(Index i = 0; i < size_out; i++) {
    //     hgres(i) = (long double)(gsl_sf_hyperg_2F1(a1s(i), 1, bs(i), trD1d));
    // }

    // Calling gsl::hyperg_2F1 from R
    Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
    Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
    Eigen::Map<ArrayXd> hgresd(hgv.begin(), size_out);
    ArrayXl hgres = hgresd.cast<long double>();

    ansmat *= hgres;

    ArrayXl ansseq = sum_counterdiagE(ansmat);
    return Rcpp::List::create(
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("diminished") = diminished);
}


// Coefficient-wise scaling will not work properly
//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP A1B1_Ec(const Eigen::ArrayXd D1, const Eigen::ArrayXd D2,
             const Eigen::ArrayXd mu1, const Eigen::ArrayXd mu2,
             const Eigen::Index m, const double thr_margin = 100,
             int nthreads = 0, const double tol_zero = 2.2e-14) {
    Index n1 = D1.size();
    Index n2 = D2.size();
    double n1_ = n1;
    double n2_ = n2;
    double m_ = m;
    bool central = is_zero_E(mu1, tol_zero) && is_zero_E(mu2, tol_zero);
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
    ArrayXd lscf1;
    ArrayXd lscf2;
    ArrayXd dk1(m + 1);
    ArrayXd dk2(m + 1);
    ArrayXd ansmat = ArrayXd::Zero(size_out);
    bool diminished = false;
    ArrayXd Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXd Aldenj = get_lgm(n2_ / 2, m + 1);
    if(central) {
        lscf1 = ArrayXd::Zero(m + 1);
        lscf2 = ArrayXd::Zero(m + 1);
        dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
        dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
        for(Index k = 0; k <= m; k++) {
            double k_ = k;
            ansmat.ULTcol(k, m + 1) +=
                log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
                log(dk2(k)) - Aldenj(k) - lscf2(k);
            double lnum = std::lgamma((n1_ + n2_) / 2 + k_);
            for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) += lnum;
        }
    } else {
        lscf1 = ArrayXd::Zero(size_out);
        lscf2 = ArrayXd::Zero(size_out);
        VectorXd mu1d = D1d.sqrt() * mu1;
        VectorXd mu2d = D2d.sqrt() * mu2;
        MatrixXd mu1mat = mu1d * mu1d.transpose() / 2;
        MatrixXd mu2mat = mu2d * mu2d.transpose() / 2;
        DiagMatXd D1hmat = D1h.matrix().asDiagonal();
        DiagMatXd D2hmat = D2h.matrix().asDiagonal();
        ArrayXd dkm1 = d2_ij_mEc(mu1mat, D1hmat, m, lscf1, thr_margin, nthreads);
        ArrayXd dkm2 = d2_ij_mEc(mu2mat, D2hmat, m, lscf2, thr_margin, nthreads);
        dkm1 = log(dkm1);
        dkm2 = log(dkm2);
        ArrayXd seqlrf_1_2 = get_lrf(0.5, m + 1);
        for(Index k = 0; k <= m; k++) {
            dkm1.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            dkm2.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
        }
        dkm1 -= lscf1;
        dkm2 -= lscf2;
        dkm1 = exp(dkm1);
        dkm2 = exp(dkm2);
        dk1 = sum_counterdiagE(dkm1);
        dk2 = sum_counterdiagE(dkm2);
        diminished = ((lscf1 < 0).any() && dkm1.cwiseEqual(0).any()) ||
                     ((lscf2 < 0).any() && dkm2.cwiseEqual(0).any());
        for(Index k = 0; k <= m; k++) {
            double k_ = k;
            ansmat.ULTcol(k, m + 1) +=
                log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) +
                log(dk2(k)) - Aldenj(k);
            double lnum = std::lgamma((n1_ + n2_) / 2 + k_);
            for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) += lnum;
        }
    }
    ansmat += (log(D1d).sum() + log(D2d).sum()) / 2;
    ansmat -= (mu1.matrix().squaredNorm() + mu2.matrix().squaredNorm()) / 2;
    ansmat = exp(ansmat);

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
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("diminished") = diminished);
}
