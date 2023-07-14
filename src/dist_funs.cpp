#include "config.h"
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

// These are to use gsl
#include <RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

#include "dk_funs.h"
#include "hgs_funs.h"
#include "dk_funs_cwise.h"

using Eigen::exp;
using Eigen::log;
using Eigen::abs;
using Eigen::ArrayXi;
using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
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

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
get_subset(const Eigen::ArrayBase<Derived>& X, const Eigen::ArrayXi& cond) {
    Index n = X.size();
    Index outsize = cond.sum();
    Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> out(outsize);
    for(Index i = 0, j = 0; i < n; i++) {
        if(cond(i)) {
            out(j) = X(i);
            j++;
        }
    }
    return out;
}

void check_hgstatus(Eigen::ArrayXi& hgstatus, const bool stop_on_error) {
    if(hgstatus.any()) {
        std::string errmsg = "problem in gsl_sf_hyperg_2F1_e():";
        bool eunimpl = hgstatus.cwiseEqual(24).any();
        bool eovrflw = hgstatus.cwiseEqual(16).any();
        bool emaxiter = hgstatus.cwiseEqual(11).any();
        bool edom = hgstatus.cwiseEqual(1).any();
        bool eother = !(eunimpl || eovrflw || emaxiter || edom);
        if(eunimpl) {
            errmsg += "\n  evaluation failed due to singularity";
        }
        if(eovrflw) {
            errmsg += "\n  numerical overflow encountered";
        }
        if(emaxiter) {
            errmsg += "\n  max iteration reached";
        }
        if(edom) {
            errmsg += "\n  parameter outside acceptable domain";
        }
        if(eother) {
            errmsg += "\n  unexpected kind of error";
        }
        if(stop_on_error) {
            Rcpp::stop(errmsg);
        } else {
            Rcpp::warning(errmsg);
        }
    }
}

//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}, double
//'
// [[Rcpp::export]]
SEXP p_A1B1_Ed(const double quantile,
               const Eigen::MatrixXd A, const Eigen::MatrixXd B,
               const Eigen::ArrayXd mu,
               const Eigen::Index m, const bool stop_on_error,
               const double thr_margin = 100,
               int nthreads = 0, const double tol_zero = 2.2e-14) {
    ArrayXd ansseq = ArrayXd::Zero(m + 1);
    bool diminished = false;
    bool exact = false;
    MatrixXd A_qB =  A - quantile * B;
    SelfAdjointEigenSolver<MatrixXd> eigA_qB(A_qB);
    const ArrayXd L = eigA_qB.eigenvalues();
    const ArrayXi ind_pos = (L >  tol_zero).cast<int>();
    const ArrayXi ind_neg = (L < -tol_zero).cast<int>();
    const ArrayXd D1 =  get_subset(L, ind_pos);
    const ArrayXd D2 = -get_subset(L, ind_neg);
    Index n1 = D1.size();
    Index n2 = D2.size();
    if(n2 == 0) {
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("diminished") = diminished,
            Rcpp::Named("exact")      = exact);
    }
    if(n1 == 0) {
        ansseq(0) = 1;
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("diminished") = diminished,
            Rcpp::Named("exact")      = exact);
    }
    const MatrixXd U = eigA_qB.eigenvectors();
    const ArrayXd nu = U.transpose() * mu.matrix();
    const ArrayXd mu1 = get_subset(nu, ind_pos);
    const ArrayXd mu2 = get_subset(nu, ind_neg);
    double n1_ = n1;
    double n2_ = n2;
    double m_ = m;
    bool n1_is_1 = (n1 == 1);
    bool n2_is_1 = (n2 == 1);
    bool cent_mu1 = is_zero_E(mu1, tol_zero);
    bool cent_mu2 = is_zero_E(mu2, tol_zero);
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
    ArrayXd dk1 = ArrayXd::Zero(m + 1);
    ArrayXd dk2 = ArrayXd::Zero(m + 1);
    ArrayXd Alnum = get_lgm((n1_ + n2_) / 2, m + 1);
    ArrayXd Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXd Aldenj = get_lgm(n2_ / 2, m + 1);
    if(cent_mu1) {
        if(n1_is_1) {
            dk1(0) = 1;
        } else {
            dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
        }
    } else {
        VectorXd mu1d = D1d.sqrt() * mu1;
        MatrixXd mu1mat = mu1d * mu1d.transpose() / 2;
        ArrayXd seqlrf_1_2 = get_lrf(0.5, m + 1);
        if(n1_is_1) {
            dk1 = d1_i_mE(mu1mat, m, lscf1, thr_margin);
            dk1 = log(dk1);
            dk1 -= seqlrf_1_2;
            dk1 = exp(dk1);
        } else {
            DiagMatXd D1hmat = D1h.matrix().asDiagonal();
            ArrayXd dkm1 = d2_ij_mE(mu1mat, D1hmat, m, lscf1, thr_margin, nthreads);
            dkm1 = log(dkm1);
            for(Index k = 0; k <= m; k++) {
                dkm1.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            }
            dkm1 = exp(dkm1);
            dk1 = sum_counterdiagE(dkm1);
            diminished = diminished || ((lscf1 < 0).any() && dkm1.cwiseEqual(0).any());
        }
    }
    if(cent_mu2) {
        if(n2_is_1) {
            dk2(0) = 1;
        } else {
            dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
        }
    } else {
        VectorXd mu2d = D2d.sqrt() * mu2;
        MatrixXd mu2mat = mu2d * mu2d.transpose() / 2;
        ArrayXd seqlrf_1_2 = get_lrf(0.5, m + 1);
        if(n2_is_1) {
            dk2 = d1_i_mE(mu2mat, m, lscf2, thr_margin);
            dk2 = log(dk2);
            dk2 -= seqlrf_1_2;
            dk2 = exp(dk2);
        } else {
            DiagMatXd D2hmat = D2h.matrix().asDiagonal();
            ArrayXd dkm2 = d2_ij_mE(mu2mat, D2hmat, m, lscf2, thr_margin, nthreads);
            dkm2 = log(dkm2);
            for(Index k = 0; k <= m; k++) {
                dkm2.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            }
            dkm2 = exp(dkm2);
            dk2 = sum_counterdiagE(dkm2);
            diminished = diminished || ((lscf2 < 0).any() && dkm2.cwiseEqual(0).any());
        }
    }
    if((n1_is_1 && cent_mu1) || (n2_is_1 && cent_mu2)) {
        ArrayXd a1s = ArrayXd::LinSpaced(m + 1, (n1_ + n2_) / 2, (n1_ + n2_) / 2 + m_);
        ArrayXd bs(m + 1);

        if(n1_is_1 && cent_mu1) {
            ansseq += log(dk2) - Aldenj - lscf2 + Alnum;
            ansseq -= Aldeni(0) + lscf1(0);
            ansseq -= (mu2.matrix().squaredNorm()) / 2;
            bs = ArrayXd::Constant(m + 1, n1_ / 2 + 1);
        } else if(n2_is_1 && cent_mu2) {
            ansseq += log(dk1) - Aldeni - lscf1 + Alnum;
            ansseq -= Aldenj(0) + lscf2(0);
            ansseq -= (mu1.matrix().squaredNorm()) / 2;
            bs = ArrayXd::LinSpaced(m + 1, n1_ / 2 + 1, n1_ / 2 + 1 + m_);
        }
        ansseq += (log(D1d).sum() + log(D2d).sum()) / 2;
        ansseq = exp(ansseq);

        // Using the C++ library GSL with RcppGSL
        ArrayXd hgres(m + 1);
        ArrayXi hgstatus(m + 1);
        gsl_sf_result hgtmp;
        gsl_set_error_handler_off();
        for(Index i = 0; i < m + 1; i++) {
            hgstatus(i) = gsl_sf_hyperg_2F1_e(a1s(i), 1, bs(i), trD1d, &hgtmp);
            hgres(i) = hgtmp.val;
        }
        check_hgstatus(hgstatus, stop_on_error);

        // // Calling gsl::hyperg_2F1 from R
        // Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
        // Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
        // Eigen::Map<ArrayXd> hgres(hgv.begin(), m + 1);

        ansseq *= hgres;
    } else {
        Index size_out = (m + 1) * (m + 2) / 2;
        ArrayXd ansmat = ArrayXd::Zero(size_out);
        for(Index k = 0; k <= m; k++) {
            ansmat.ULTcol(k, m + 1) += Alnum.tail(m + 1 - k) +
                log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
                log(dk2(k)) - Aldenj(k) - lscf2(k);
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

        // Using the C++ library GSL with RcppGSL
        ArrayXd hgres(size_out);
        ArrayXi hgstatus(size_out);
        gsl_sf_result hgtmp;
        gsl_set_error_handler_off();
        for(Index i = 0; i < size_out; i++) {
            hgstatus(i) = gsl_sf_hyperg_2F1_e(a1s(i), 1, bs(i), trD1d, &hgtmp);
            hgres(i) = hgtmp.val;
        }
        check_hgstatus(hgstatus, stop_on_error);

        // // Calling gsl::hyperg_2F1 from R
        // Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
        // Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
        // Eigen::Map<ArrayXd> hgres(hgv.begin(), size_out);

        ansmat *= hgres;

        ansseq = sum_counterdiagE(ansmat);
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("diminished") = diminished,
        Rcpp::Named("exact")      = exact);
}


//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}, long double
//'
// [[Rcpp::export]]
SEXP p_A1B1_El(const long double quantile,
               const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
               const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> B,
               const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
               const Eigen::Index m, const bool stop_on_error,
               const long double thr_margin = 100,
               int nthreads = 0, const long double tol_zero = 2.2e-14) {
    ArrayXl ansseq = ArrayXl::Zero(m + 1);
    bool diminished = false;
    bool exact = false;
    MatrixXl A_qB =  A - quantile * B;
    SelfAdjointEigenSolver<MatrixXl> eigA_qB(A_qB);
    const ArrayXl L = eigA_qB.eigenvalues();
    const ArrayXi ind_pos = (L >  tol_zero).cast<int>();
    const ArrayXi ind_neg = (L < -tol_zero).cast<int>();
    const ArrayXl D1 =  get_subset(L, ind_pos);
    const ArrayXl D2 = -get_subset(L, ind_neg);
    Index n1 = D1.size();
    Index n2 = D2.size();
    if(n2 == 0) {
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("diminished") = diminished,
            Rcpp::Named("exact")      = exact);
    }
    if(n1 == 0) {
        ansseq(0) = 1;
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("diminished") = diminished,
            Rcpp::Named("exact")      = exact);
    }
    const MatrixXl U = eigA_qB.eigenvectors();
    const ArrayXl nu = U.transpose() * mu.matrix();
    const ArrayXl mu1 = get_subset(nu, ind_pos);
    const ArrayXl mu2 = get_subset(nu, ind_neg);
    long double n1_ = n1;
    long double n2_ = n2;
    long double m_ = m;
    bool n1_is_1 = (n1 == 1);
    bool n2_is_1 = (n2 == 1);
    bool cent_mu1 = is_zero_E(mu1, tol_zero);
    bool cent_mu2 = is_zero_E(mu2, tol_zero);
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
    ArrayXl dk1 = ArrayXl::Zero(m + 1);
    ArrayXl dk2 = ArrayXl::Zero(m + 1);
    ArrayXl Alnum = get_lgm((n1_ + n2_) / 2, m + 1);
    ArrayXl Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXl Aldenj = get_lgm(n2_ / 2, m + 1);
    if(cent_mu1) {
        if(n1_is_1) {
            dk1(0) = 1;
        } else {
            dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
        }
    } else {
        VectorXl mu1d = D1d.sqrt() * mu1;
        MatrixXl mu1mat = mu1d * mu1d.transpose() / 2;
        ArrayXl seqlrf_1_2 = get_lrf((long double)(0.5), m + 1);
        if(n1_is_1) {
            dk1 = d1_i_mE(mu1mat, m, lscf1, thr_margin);
            dk1 = log(dk1);
            dk1 -= seqlrf_1_2;
            dk1 = exp(dk1);
        } else {
            DiagMatXl D1hmat = D1h.matrix().asDiagonal();
            ArrayXl dkm1 = d2_ij_mE(mu1mat, D1hmat, m, lscf1, thr_margin, nthreads);
            dkm1 = log(dkm1);
            for(Index k = 0; k <= m; k++) {
                dkm1.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            }
            dkm1 = exp(dkm1);
            dk1 = sum_counterdiagE(dkm1);
            diminished = diminished || ((lscf1 < 0).any() && dkm1.cwiseEqual(0).any());
        }
    }
    if(cent_mu2) {
        if(n2_is_1) {
            dk2(0) = 1;
        } else {
            dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
        }
    } else {
        VectorXl mu2d = D2d.sqrt() * mu2;
        MatrixXl mu2mat = mu2d * mu2d.transpose() / 2;
        ArrayXl seqlrf_1_2 = get_lrf((long double)(0.5), m + 1);
        if(n2_is_1) {
            dk2 = d1_i_mE(mu2mat, m, lscf2, thr_margin);
            dk2 = log(dk2);
            dk2 -= seqlrf_1_2;
            dk2 = exp(dk2);
        } else {
            DiagMatXl D2hmat = D2h.matrix().asDiagonal();
            ArrayXl dkm2 = d2_ij_mE(mu2mat, D2hmat, m, lscf2, thr_margin, nthreads);
            dkm2 = log(dkm2);
            for(Index k = 0; k <= m; k++) {
                dkm2.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            }
            dkm2 = exp(dkm2);
            dk2 = sum_counterdiagE(dkm2);
            diminished = diminished || ((lscf2 < 0).any() && dkm2.cwiseEqual(0).any());
        }
    }
    if((n1_is_1 && cent_mu1) || (n2_is_1 && cent_mu2)) {
        ArrayXl a1s = ArrayXl::LinSpaced(m + 1, (n1_ + n2_) / 2, (n1_ + n2_) / 2 + m_);
        ArrayXl bs(m + 1);

        if(n1_is_1 && cent_mu1) {
            ansseq += log(dk2) - Aldenj - lscf2 + Alnum;
            ansseq -= Aldeni(0) + lscf1(0);
            ansseq -= (mu2.matrix().squaredNorm()) / 2;
            bs = ArrayXl::Constant(m + 1, n1_ / 2 + 1);
        } else if(n2_is_1 && cent_mu2) {
            ansseq += log(dk1) - Aldeni - lscf1 + Alnum;
            ansseq -= Aldenj(0) + lscf2(0);
            ansseq -= (mu1.matrix().squaredNorm()) / 2;
            bs = ArrayXl::LinSpaced(m + 1, n1_ / 2 + 1, n1_ / 2 + 1 + m_);
        }
        ansseq += (log(D1d).sum() + log(D2d).sum()) / 2;
        ansseq = exp(ansseq);

        // Using the C++ library GSL with RcppGSL
        ArrayXl hgres(m + 1);
        ArrayXi hgstatus(m + 1);
        gsl_sf_result hgtmp;
        gsl_set_error_handler_off();
        for(Index i = 0; i < m + 1; i++) {
            hgstatus(i) = gsl_sf_hyperg_2F1_e(a1s(i), 1, bs(i), trD1d, &hgtmp);
            hgres(i) = (long double)(hgtmp.val);
        }
        check_hgstatus(hgstatus, stop_on_error);

        // // Calling gsl::hyperg_2F1 from R
        // Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
        // Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
        // Eigen::Map<ArrayXd> hgresd(hgv.begin(), m + 1);
        // ArrayXl hgres = hgresd.cast<long double>();

        ansseq *= hgres;
    } else {
        Index size_out = (m + 1) * (m + 2) / 2;
        ArrayXl ansmat = ArrayXl::Zero(size_out);
        for(Index k = 0; k <= m; k++) {
            ansmat.ULTcol(k, m + 1) += Alnum.tail(m + 1 - k) +
                log(dk1.head(m + 1 - k)) - Aldeni.head(m + 1 - k) - lscf1.head(m + 1 - k) +
                log(dk2(k)) - Aldenj(k) - lscf2(k);
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

        // Using the C++ library GSL with RcppGSL
        ArrayXl hgres(size_out);
        ArrayXi hgstatus(size_out);
        gsl_sf_result hgtmp;
        gsl_set_error_handler_off();
        for(Index i = 0; i < size_out; i++) {
            hgstatus(i) = gsl_sf_hyperg_2F1_e(a1s(i), 1, bs(i), trD1d, &hgtmp);
            hgres(i) = (long double)(hgtmp.val);
        }
        check_hgstatus(hgstatus, stop_on_error);

        // // Calling gsl::hyperg_2F1 from R
        // Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
        // Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
        // Eigen::Map<ArrayXd> hgresd(hgv.begin(), size_out);
        // ArrayXl hgres = hgresd.cast<long double>();

        ansmat *= hgres;

        ansseq = sum_counterdiagE(ansmat);
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("diminished") = diminished,
        Rcpp::Named("exact")      = exact);
}


// Coefficient-wise scaling will not work properly
//' @describeIn qfrm_cpp
//'   \code{pqfm_A1B1()}, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP p_A1B1_Ec(const double quantile,
               const Eigen::MatrixXd A, const Eigen::MatrixXd B,
               const Eigen::ArrayXd mu,
               const Eigen::Index m, const bool stop_on_error,
               const double thr_margin = 100,
               int nthreads = 0, const double tol_zero = 2.2e-14) {
    ArrayXd ansseq = ArrayXd::Zero(m + 1);
    bool diminished = false;
    bool exact = false;
    MatrixXd A_qB =  A - quantile * B;
    SelfAdjointEigenSolver<MatrixXd> eigA_qB(A_qB);
    const ArrayXd L = eigA_qB.eigenvalues();
    const ArrayXi ind_pos = (L >  tol_zero).cast<int>();
    const ArrayXi ind_neg = (L < -tol_zero).cast<int>();
    const ArrayXd D1 =  get_subset(L, ind_pos);
    const ArrayXd D2 = -get_subset(L, ind_neg);
    Index n1 = D1.size();
    Index n2 = D2.size();
    if(n2 == 0) {
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("diminished") = diminished,
            Rcpp::Named("exact")      = exact);
    }
    if(n1 == 0) {
        ansseq(0) = 1;
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("diminished") = diminished,
            Rcpp::Named("exact")      = exact);
    }
    const MatrixXd U = eigA_qB.eigenvectors();
    const ArrayXd nu = U.transpose() * mu.matrix();
    const ArrayXd mu1 = get_subset(nu, ind_pos);
    const ArrayXd mu2 = get_subset(nu, ind_neg);
    double n1_ = n1;
    double n2_ = n2;
    double m_ = m;
    bool n1_is_1 = (n1 == 1);
    bool n2_is_1 = (n2 == 1);
    bool cent_mu1 = is_zero_E(mu1, tol_zero);
    bool cent_mu2 = is_zero_E(mu2, tol_zero);
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
    ArrayXd dk1 = ArrayXd::Zero(m + 1);
    ArrayXd dk2 = ArrayXd::Zero(m + 1);
    Index size_out = (m + 1) * (m + 2) / 2;
    ArrayXd Alnum = get_lgm((n1_ + n2_) / 2, m + 1);
    ArrayXd Aldeni = get_lgm(n1_ / 2 + 1, m + 1);
    ArrayXd Aldenj = get_lgm(n2_ / 2, m + 1);
    if(cent_mu1) {
        if(n1_is_1) {
            // In this case dk1 = 0 except 0th; take log for later
            dk1 = -INFINITY;
            dk1(0) = 0;
        } else {
            lscf1 = ArrayXd::Zero(m + 1);
            dk1 = d1_i_vE(D1h, m, lscf1, thr_margin);
            dk1 = log(dk1);
            dk1 -= lscf1;
        }
    } else {
        VectorXd mu1d = D1d.sqrt() * mu1;
        MatrixXd mu1mat = mu1d * mu1d.transpose() / 2;
        ArrayXd seqlrf_1_2 = get_lrf(0.5, m + 1);
        if(n1_is_1) {
            lscf1 = ArrayXd::Zero(m + 1);
            dk1 = d1_i_mE(mu1mat, m, lscf1, thr_margin);
            dk1 = log(dk1);
            dk1 -= seqlrf_1_2 + lscf1;
        } else {
            lscf1 = ArrayXd::Zero(size_out);
            DiagMatXd D1hmat = D1h.matrix().asDiagonal();
            ArrayXd dkm1 = d2_ij_mEc(mu1mat, D1hmat, m, lscf1, thr_margin, nthreads);
            dkm1 = log(dkm1);
            for(Index k = 0; k <= m; k++) {
                dkm1.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            }
            dkm1 -= lscf1;
            dkm1 = exp(dkm1);
            dk1 = sum_counterdiagE(dkm1);
            dk1 = log(dk1);
            diminished = diminished || (((lscf1 < 0) && dkm1.cwiseEqual(0)).any());
        }
    }
    if(cent_mu2) {
        if(n2_is_1) {
            dk2 = log(dk2);
            dk2(0) = 0;
        } else {
            lscf2 = ArrayXd::Zero(m + 1);
            dk2 = d1_i_vE(D2h, m, lscf2, thr_margin);
            dk2 = log(dk2);
            dk2 -= lscf2;
        }
    } else {
        VectorXd mu2d = D2d.sqrt() * mu2;
        MatrixXd mu2mat = mu2d * mu2d.transpose() / 2;
        ArrayXd seqlrf_1_2 = get_lrf(0.5, m + 1);
        if(n2_is_1) {
            lscf2 = ArrayXd::Zero(m + 1);
            dk2 = d1_i_mE(mu2mat, m, lscf2, thr_margin);
            dk2 = log(dk2);
            dk2 -= seqlrf_1_2 + lscf2;
        } else {
            lscf2 = ArrayXd::Zero(size_out);
            DiagMatXd D2hmat = D2h.matrix().asDiagonal();
            ArrayXd dkm2 = d2_ij_mEc(mu2mat, D2hmat, m, lscf2, thr_margin, nthreads);
            dkm2 = log(dkm2);
            for(Index k = 0; k <= m; k++) {
                dkm2.ULTcol(k, m + 1) -= seqlrf_1_2.head(m + 1 - k);
            }
            dkm2 -= lscf2;
            dkm2 = exp(dkm2);
            dk2 = sum_counterdiagE(dkm2);
            dk2 = log(dk2);
            diminished = diminished || (((lscf2 < 0) && dkm2.cwiseEqual(0)).any());
        }
    }
    if((n1_is_1 && cent_mu1) || (n2_is_1 && cent_mu2)) {
        ArrayXd a1s = ArrayXd::LinSpaced(m + 1, (n1_ + n2_) / 2, (n1_ + n2_) / 2 + m_);
        ArrayXd bs(m + 1);

        if(n1_is_1 && cent_mu1) {
            ansseq += dk2 - Aldenj + Alnum;
            ansseq -= Aldeni(0); // Or dk1(0)
            ansseq -= (mu2.matrix().squaredNorm()) / 2;
            bs = ArrayXd::Constant(m + 1, n1_ / 2 + 1);
        } else if(n2_is_1 && cent_mu2) {
            ansseq += dk1 - Aldeni + Alnum;
            ansseq -= Aldenj(0);
            ansseq -= (mu1.matrix().squaredNorm()) / 2;
            bs = ArrayXd::LinSpaced(m + 1, n1_ / 2 + 1, n1_ / 2 + 1 + m_);
        }
        ansseq += (log(D1d).sum() + log(D2d).sum()) / 2;
        ansseq = exp(ansseq);

        // Using the C++ library GSL with RcppGSL
        ArrayXd hgres(m + 1);
        ArrayXi hgstatus(m + 1);
        gsl_sf_result hgtmp;
        gsl_set_error_handler_off();
        for(Index i = 0; i < m + 1; i++) {
            hgstatus(i) = gsl_sf_hyperg_2F1_e(a1s(i), 1, bs(i), trD1d, &hgtmp);
            hgres(i) = hgtmp.val;
        }
        check_hgstatus(hgstatus, stop_on_error);

        // // Calling gsl::hyperg_2F1 from R
        // Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
        // Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
        // Eigen::Map<ArrayXd> hgres(hgv.begin(), m + 1);

        ansseq *= hgres;
    } else {
        ArrayXd ansmat = ArrayXd::Zero(size_out);
        for(Index k = 0; k <= m; k++) {
            ansmat.ULTcol(k, m + 1) += Alnum.tail(m + 1 - k) +
                dk1.head(m + 1 - k) - Aldeni.head(m + 1 - k) +
                dk2(k) - Aldenj(k);
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

        // Using the C++ library GSL with RcppGSL
        ArrayXd hgres(size_out);
        ArrayXi hgstatus(size_out);
        gsl_sf_result hgtmp;
        gsl_set_error_handler_off();
        for(Index i = 0; i < size_out; i++) {
            hgstatus(i) = gsl_sf_hyperg_2F1_e(a1s(i), 1, bs(i), trD1d, &hgtmp);
            hgres(i) = hgtmp.val;
        }
        check_hgstatus(hgstatus, stop_on_error);

        // // Calling gsl::hyperg_2F1 from R
        // Rcpp::Function hyperg_2F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_2F1"];
        // Rcpp::NumericVector hgv = hyperg_2F1(a1s, 1, bs, trD1d);
        // Eigen::Map<ArrayXd> hgres(hgv.begin(), size_out);

        ansmat *= hgres;

        ansseq = sum_counterdiagE(ansmat);
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("diminished") = diminished,
        Rcpp::Named("exact")      = exact);
}


//' @describeIn qfrm_cpp
//'   \code{dqfm_A1I1()}
//'
// [[Rcpp::export]]
SEXP d_A1I1_Ed(const double quantile, const Eigen::ArrayXd LA,
               const Eigen::Index m, const double thr_margin = 100) {
    const double n_ = double(LA.size());
    const double L1 = LA.maxCoeff();
    const double Ls = LA.minCoeff();
    ArrayXd ansseq = ArrayXd::Zero(m);
    bool exact = false;
    if(quantile >= L1 || quantile <= Ls) {
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("exact")  = exact);
    }
    const double f = (quantile - Ls) / (L1 - quantile);
    const double n1_ = double((LA == L1).cast<int>().sum());
    const double ns_ = double((LA == Ls).cast<int>().sum());
    if(n1_ + ns_ == n_) {
        ansseq(0) = std::pow(f, n1_ / 2.0 - 1.0) / std::pow(1 + f, n_ / 2.0) /
                    std::beta(n1_ / 2.0, (n_ - n1_) / 2.0);
        ansseq(0) *= (L1 - Ls) / std::pow(L1 - quantile, 2.0);
        exact = true;
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("exact")  = exact);
    }
    const ArrayXd LA_for_psi = get_subset(LA, ((LA != L1) && (LA != Ls)).cast<int>());
    const ArrayXd psi = (LA_for_psi - Ls) / (L1 - LA_for_psi);
    ArrayXi ind_r1 = (psi > f).cast<int>();
    ArrayXd psi_r1 = get_subset(psi, ind_r1);
    ArrayXd psi_r2 = get_subset(psi, 1 - ind_r1);
    double pr = ind_r1.sum() + n1_;
    if((pr == n1_) || (pr == n_ - ns_)) {
        ArrayXd D;
        double nt_;
        ArrayXd lscf = ArrayXd::Zero(m + 1);
        if(pr == n1_) {
            D = psi / f;
            nt_ = n1_;
        } else {
            D = f / psi;
            nt_ = ns_;
        }
        ArrayXd dks = d1_i_vE(D, m, lscf, thr_margin);
        ansseq = hgs_1dE(dks, (2 - nt_) / 2, (n_ - nt_) / 2, 0, lscf);
    } else {
        Index size_out = (m + 1) * (m + 2) / 2;
        ArrayXd D1 = f / psi_r1;
        ArrayXd D2 = psi_r2 / f;
        ArrayXd lscf1 = ArrayXd::Zero(m + 1);
        ArrayXd lscf2 = ArrayXd::Zero(m + 1);
        ArrayXd dk1 = d1_i_vE(D1, m, lscf1, thr_margin);
        ArrayXd dk2 = d1_i_vE(D2, m, lscf2, thr_margin);
        dk1 = dk1.log();
        dk2 = dk2.log();
        double alpha = pr / 2 - 1;
        double beta = (n_ - pr) / 2 - 1;
        ArrayXd j_minus_k = ArrayXd::LinSpaced(2 * m + 1, -m, m);
        ArrayXd ordmat(size_out);
        for(Index k = 0; k <= m; k++) {
            ordmat.ULTcol(k, m + 1) = j_minus_k.segment(m - k, m + 1 - k);
        }
        ArrayXd ansmat = ArrayXd::Zero(size_out);
        for(Index k = 0; k <= m; k++) {
            ansmat.ULTcol(k, m + 1) += dk1.head(m + 1 - k) - lscf1.head(m + 1 - k) +
                                       dk2(k) - lscf2(k);
        }
        ansmat -= (alpha + 1 + ordmat).lgamma() + (beta + 1 - ordmat).lgamma() -
                  std::lgamma(alpha + 1) - std::lgamma(beta + 1);
        ansmat = ansmat.exp();
        ArrayXi ord_mod2 = (ordmat - (2 * (ordmat / 2).floor())).cast<int>();
        ansmat = ((ordmat > 0) && (ordmat <= beta) && (ord_mod2 == 1)).select(-ansmat, ansmat);
        ansmat = ((ordmat < 0) && (ordmat >= -alpha) && (ord_mod2 == 1)).select(-ansmat, ansmat);
        if(is_int_like(beta)) {
            ansmat = (ordmat > beta).select(0, ansmat);
        } else if(int(std::floor(beta)) % 2 == 0) {
            ansmat = (ordmat > beta).select(-ansmat, ansmat);
        }
        if(is_int_like(alpha)) {
            ansmat = (ordmat < -alpha).select(0, ansmat);
        } else if(int(std::floor(alpha)) % 2 == 0) {
            ansmat = (ordmat < -alpha).select(-ansmat, ansmat);
        }
        ansseq = sum_counterdiagE(ansmat);
    }
    ansseq *= exp(std::lgamma(n_ / 2) - std::lgamma(pr / 2) - std::lgamma((n_ - pr) / 2) +
                  (-psi_r1.log().sum() + (1 + psi).log().sum() +
                   (pr - 2) * log(f) - n_ * log(1 + f)) / 2);
    ansseq *= (L1 - Ls) / std::pow(L1 - quantile, 2);
    return Rcpp::List::create(
        Rcpp::Named("ansseq") = ansseq,
        Rcpp::Named("exact")  = exact);
}

struct broda_params {
    Eigen::ArrayXd L;
    Eigen::MatrixXd H;
    Eigen::ArrayXd nu;
};

double broda_fun(double u, void *p){
    struct broda_params *params = (struct broda_params *)p;
    const ArrayXd L = (params->L);
    const MatrixXd H = (params->H);
    const ArrayXd nu = (params->nu);
    double out;
    ArrayXd a = L * u;
    ArrayXd b = a.pow(2.0);
    ArrayXd c = 1.0 + b;
    ArrayXd theta = nu.pow(2.0);
    double beta = (a.atan() + theta * a / c).sum() / 2.0;
    double gamma = std::exp((theta * b / c).sum() / 2.0 + c.log().sum() / 4.0);
    ArrayXd Fi = c.inverse();
    VectorXd Finu = Fi * nu;
    double rho = (H * Fi.matrix().asDiagonal()).trace() +
                 Finu.transpose() *
                 (H - a.matrix().asDiagonal() * H * a.matrix().asDiagonal()) *
                 Finu;
    double delta = (H * (L * Fi).matrix().asDiagonal()).trace() +
                   2.0 * Finu.transpose() * H * L.matrix().asDiagonal() * Finu;
    out = (rho * std::cos(beta) - u * delta * std::sin(beta)) / gamma / 2.0;
    return out;
}

//' @describeIn qfrm_cpp
//'   \code{dqfm_broda()}
//'
// [[Rcpp::export]]
SEXP d_broda_Ed(const double quantile,
                const Eigen::MatrixXd A, const Eigen::MatrixXd B,
                const Eigen::ArrayXd mu,
                double autoscale_args, bool stop_on_error, double tol_zero,
                double epsabs, double epsrel, int limit) {
    MatrixXd A_qB =  A - quantile * B;
    SelfAdjointEigenSolver<MatrixXd> eigA_qB(A_qB);
    ArrayXd L = eigA_qB.eigenvalues();
    if((L >= -tol_zero).all() || (L <= tol_zero).all()) {
        return Rcpp::List::create(
            Rcpp::Named("value") = 0,
            Rcpp::Named("abs.error") = 0);
    }
    const MatrixXd U = eigA_qB.eigenvectors();
    const ArrayXd nu = U.transpose() * mu.matrix();
    MatrixXd H = U.transpose() * B * U;
    if(autoscale_args > 0.0) {
        double Labsmax = L.abs().maxCoeff() / autoscale_args;
        L /= Labsmax;
        H /= Labsmax;
    }
    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);
    double value, error;
    int status;
    struct broda_params params;
    params.L = L;
    params.H = H;
    params.nu = nu;
    gsl_function F;
    F.function = &broda_fun;
    F.params = &params;
    status = gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if(status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        if(stop_on_error) {
            Rcpp::stop(errmsg);
        } else {
            Rcpp::warning(errmsg);
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("value") = value,
        Rcpp::Named("abs.error") = error);
}

double Kder_fun(const Eigen::ArrayXd& Xii, const Eigen::ArrayXd& L,
                const ArrayXd& theta, double j) {
    double out = ((L * Xii).pow(j) * (1.0 + j * theta * Xii)).sum();
    out *= std::pow(2.0, j - 1.0) * std::tgamma(j);
    return out;
}

struct mgf_params {
    Eigen::ArrayXd L;
    Eigen::ArrayXd theta;
};

double Kp1_gslfun(double s, void *p) {
    struct mgf_params *params = (struct mgf_params *)p;
    const ArrayXd L = (params->L);
    const ArrayXd theta = (params->theta);
    ArrayXd Xii = (1.0 - 2.0 * s * L).inverse();
    double out = Kder_fun(Xii, L, theta, 1.0);
    return out;
}

// double Kp2_gslfun(double s, void *p) {
//     struct mgf_params *params = (struct mgf_params *)p;
//     const ArrayXd L = (params->L);
//     const ArrayXd theta = (params->theta);
//     ArrayXd Xii = (1.0 - 2.0 * s * L).inverse();
//     double out = Kder_fun(Xii, L, theta, 2.0);
//     return out;
// }

// void Kp1_Kp2_gslfun(double s, void *p, double *y, double *dy) {
//     struct mgf_params *params = (struct mgf_params *)p;
//     const ArrayXd L = (params->L);
//     const ArrayXd theta = (params->theta);
//     ArrayXd Xii = (1.0 - 2.0 * s * L).inverse();
//     *y = Kder_fun(Xii, L, theta, 1.0);
//     *dy = Kder_fun(Xii, L, theta, 2.0);
// }

double Kx_fun(double s, const Eigen::ArrayXd& L,
              const ArrayXd& theta, const Eigen::ArrayXd& Xii) {
    double out = (Xii.log() / 2.0 + s * L * theta * Xii).sum();
    return out;
}

double Mx_fun(double s, const Eigen::ArrayXd& L,
              const ArrayXd& theta, const Eigen::ArrayXd& Xii) {
    double out = Kx_fun(s, L, theta, Xii);
    return std::exp(out);
}

double J_fun(const Eigen::ArrayXd& Xii, const Eigen::ArrayXd& L,
             const MatrixXd& H, const VectorXd& Xiimu) {
    double out;
    out = (Xii * H.diagonal().array()).sum() + Xiimu.transpose() * H * Xiimu;
    return out;
}

double Jp1_fun(const Eigen::ArrayXd& Xii, const Eigen::ArrayXd& L,
               const MatrixXd& H, const VectorXd& Xiimu) {
    double out;
    out = 2.0 * (Xii.pow(2.0) * L * H.diagonal().array()).sum() +
          4.0 * Xiimu.transpose() * (L * Xii).matrix().asDiagonal() * H * Xiimu;
    return out;
}

double Jp2_fun(const Eigen::ArrayXd& Xii, const Eigen::ArrayXd& L,
               const MatrixXd& H, const VectorXd& Xiimu) {
    ArrayXd XiiL = Xii * L;
    double out;
    out =  8.0 * (XiiL.pow(2.0) * Xii * H.diagonal().array()).sum() +
          16.0 * Xiimu.transpose() * XiiL.pow(2.0).matrix().asDiagonal() * H * Xiimu +
           8.0 * Xiimu.transpose() * XiiL.matrix().asDiagonal() * H * XiiL.matrix().asDiagonal() * Xiimu;
    return out;
}

int butler_spa_root_find(double& s,
                         const Eigen::ArrayXd& L, const Eigen::ArrayXd& theta,
                         double epsabs, double epsrel, int maxiter,
                         bool stop_on_error) {
    gsl_set_error_handler_off();
    double s_lo = 0.5 / L.minCoeff() + epsabs;
    double s_hi = 0.5 / L.maxCoeff() - epsabs;
    struct mgf_params params;
    params.L = L;
    params.theta = theta;
    // root finding with Brent method
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *solver = gsl_root_fsolver_alloc(T);
    gsl_function F;
    F.function = &Kp1_gslfun;
    F.params = &params;
    gsl_root_fsolver_set(solver, &F, s_lo, s_hi);
    // // Alternative: Newton-Raphson method; sometimes fails with zero derivative
    // const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
    // gsl_root_fdfsolver *solver = gsl_root_fdfsolver_alloc(T);
    // gsl_function_fdf FDF;
    // FDF.f = &Kp1_gslfun;
    // FDF.df = &Kp2_gslfun;
    // FDF.fdf = &Kp1_Kp2_gslfun;
    // FDF.params = &params;
    // s = (s_lo + s_hi) / 2.0;
    // double s0;
    // gsl_root_fdfsolver_set(solver, &FDF, s);
    int status_solver, status_stop;
    int iter = 0;
    do {
        iter++;
        status_solver = gsl_root_fsolver_iterate(solver);
        s_lo = gsl_root_fsolver_x_lower(solver);
        s_hi = gsl_root_fsolver_x_upper(solver);
        status_stop = gsl_root_test_delta(s_hi, s_lo, epsabs, epsrel);
        // s0 = s;
        // status_solver = gsl_root_fdfsolver_iterate(solver);
        // s = gsl_root_fdfsolver_root(solver);
        // status_stop = gsl_root_test_delta(s, s0, epsabs, epsrel);
    } while(status_solver == GSL_SUCCESS && status_stop == GSL_CONTINUE && iter < maxiter);
    s = gsl_root_fsolver_root(solver);
    gsl_root_fsolver_free(solver);
    // gsl_root_fdfsolver_free(solver);
    if(status_solver) {
        std::string errmsg_solver = "problem in gsl_root_fsolver_iterate:\n  ";
        errmsg_solver += gsl_strerror(status_solver);
        if(stop_on_error) {
            Rcpp::stop(errmsg_solver);
        } else {
            Rcpp::warning(errmsg_solver);
        }
    }
    if(status_stop) {
        std::string errmsg_stop = "problem in gsl_root_test_delta():\n  ";
        errmsg_stop += gsl_strerror(status_stop);
        if(stop_on_error) {
            Rcpp::stop(errmsg_stop);
        } else {
            Rcpp::warning(errmsg_stop);
        }
    }
    return status_solver;
}

//' @describeIn qfrm_cpp
//'   \code{dqfm_butler()}
//'
// [[Rcpp::export]]
SEXP d_butler_Ed(const double quantile,
                 const Eigen::MatrixXd A, const Eigen::MatrixXd B,
                 const Eigen::ArrayXd mu, int order_spa, bool stop_on_error,
                 double tol_zero, double epsabs, double epsrel, int maxiter) {
    MatrixXd A_qB =  A - quantile * B;
    SelfAdjointEigenSolver<MatrixXd> eigA_qB(A_qB);
    const ArrayXd L = eigA_qB.eigenvalues();
    if((L >= -tol_zero).all() || (L <= tol_zero).all()) {
        return Rcpp::List::create(
            Rcpp::Named("value") = 0);
    }
    const MatrixXd U = eigA_qB.eigenvectors();
    const ArrayXd nu = U.transpose() * mu.matrix();
    const ArrayXd theta = nu.pow(2.0);
    const MatrixXd H = U.transpose() * B * U;
    double s;
    int status;
    status = butler_spa_root_find(s, L, theta, epsabs, epsrel, maxiter, stop_on_error);
    ArrayXd Xii_s = (1.0 - 2.0 * s * L).inverse();
    VectorXd Xiinu = Xii_s * nu;
    double J_s = J_fun(Xii_s, L, H, Xiinu);
    double Kp2_s = Kder_fun(Xii_s, L, theta, 2.0);
    double Mx_s = Mx_fun(s, L, theta, Xii_s);
    double value = Mx_s * J_s / sqrt(M_2PI * Kp2_s);
    if(order_spa > 1) {
        double Kp3_s = Kder_fun(Xii_s, L, theta, 3.0);
        double Kp4_s = Kder_fun(Xii_s, L, theta, 4.0);
        double Jp1_s = Jp1_fun(Xii_s, L, H, Xiinu);
        double Jp2_s = Jp2_fun(Xii_s, L, H, Xiinu);
        double k3h = Kp3_s / std::pow(Kp2_s, 1.5);
        double k4h = Kp4_s / std::pow(Kp2_s, 2.0);
        double cf = (k4h / 8.0 - 5.0 / 24.0 * std::pow(k3h, 2.0) +
                     Jp1_s * k3h / 2.0 / J_s / std::sqrt(Kp2_s) -
                     Jp2_s / 2.0 / J_s / Kp2_s);
        value *= 1.0 + cf;
    }
    return Rcpp::List::create(
        Rcpp::Named("value") = value);
}

//' @describeIn qfrm_cpp
//'   \code{pqfm_butler()}
//'
// [[Rcpp::export]]
SEXP p_butler_Ed(const double quantile,
                 const Eigen::MatrixXd A, const Eigen::MatrixXd B,
                 const Eigen::ArrayXd mu, int order_spa,
                 bool stop_on_error, double tol_zero,
                 double epsabs, double epsrel, int maxiter) {
    MatrixXd A_qB =  A - quantile * B;
    SelfAdjointEigenSolver<MatrixXd> eigA_qB(A_qB);
    const ArrayXd L = eigA_qB.eigenvalues();
    if((L >= -tol_zero).all()) {
        return Rcpp::List::create(
            Rcpp::Named("value") = 0);
    }
    if((L <= tol_zero).all()) {
        return Rcpp::List::create(
            Rcpp::Named("value") = 1);
    }
    const MatrixXd U = eigA_qB.eigenvectors();
    const ArrayXd nu = U.transpose() * mu.matrix();
    const ArrayXd theta = nu.pow(2.0);
    double value;
    double s;
    int status;
    status = butler_spa_root_find(s, L, theta, epsabs, epsrel, maxiter, stop_on_error);
    if(std::abs(s) * (1 - epsrel) < std::max(epsabs, tol_zero)) {
        ArrayXd Xii_0 = ArrayXd::Ones(L.size());
        double Kp2_0 = Kder_fun(Xii_0, L, theta, 2);
        double Kp3_0 = Kder_fun(Xii_0, L, theta, 3);
        value = 0.5 + Kp3_0 / 6.0 * M_1_SQRT_2PI / std::pow(Kp2_0, 1.5);
    } else {
        ArrayXd Xii_s = (1.0 - 2.0 * s * L).inverse();
        double K_s = Kx_fun(s, L, theta, Xii_s);
        double Kp2_s = Kder_fun(Xii_s, L, theta, 2.0);
        double w = std::copysign(std::sqrt(-2.0 * K_s), s);
        double u = s * std::sqrt(Kp2_s);
        double cf = 1.0 / w - 1.0 / u;
        if(order_spa > 1) {
            double Kp3_s = Kder_fun(Xii_s, L, theta, 3.0);
            double Kp4_s = Kder_fun(Xii_s, L, theta, 4.0);
            double k3h = Kp3_s / std::pow(Kp2_s, 1.5);
            double k4h = Kp4_s / std::pow(Kp2_s, 2.0);
            cf -= ((k4h / 8.0 - 5.0 / 24.0 * std::pow(k3h, 2.0)) / u -
                   std::pow(u, -3.0) - k3h / 2 / std::pow(u, 2.0) +
                   std::pow(w, -3.0));
        }
        value = R::pnorm(w, 0.0, 1.0, 1, 0) + R::dnorm(w, 0.0, 1.0, 0) * cf;
    }
    return Rcpp::List::create(
        Rcpp::Named("value") = value);
}
