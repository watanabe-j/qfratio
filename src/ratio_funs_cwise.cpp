#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

#include "dk_funs.h"
#include "dk_funs_cwise.h"
#include "hgs_funs.h"
#include "ratio_funs.h"

using Eigen::log;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, noncentral, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApIq_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::MatrixXd UA, const double b1,
                   const Eigen::VectorXd mu,
                   const double p = 1, const double q = 1, const int m = 100,
                   const double thr_margin = 100) {
    const int n = LA.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd mud = UA.transpose() * mu;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = h2_ij_vE(LAh, zeromat, mud, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1)
                             + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                   const double b1, const double b2,
                   const double p = 1, const double q = 1, const int m = 100,
                   const double thr_margin = 100)
{
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const double b1, const double b2,
                   const double p = 1, const double q = 1, const int m = 100, 
                   const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - b1 * A;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                   const double b1, const double b2, const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const int m = 100,
                   const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = h2_ij_vE(LAh, LBh, mu, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const double b1, const double b2, const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const int m = 100,
                   const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - b1 * A;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = h2_ij_mE(Ah, Bh, mu, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const double b2, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, bool error_bound = true,
                     const double thr_margin = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = htil3_pjk_vE(LA, LBh, zeromat, mu, m, int(p), lscf, thr_margin).row(p); // , nthreads).row(p);
    dks.resize(m + 1, m + 1);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2) + q * log(b2)
                              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / b2) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s = std::max(q, r);
        lscf.setZero();
        ArrayXXd dkstm = hhat3_pjk_vE(LAp, LBh, zeromat, mu, m, int(p), lscf, thr_margin).row(p); // , nthreads).row(p);
        dkstm.resize(m + 1, m + 1);
        dkstm /= exp(lscf - lscf.minCoeff());
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        set_cumsum(dkst, cumsum_dkst);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        ArrayXd LAp_Bb2 = LAp / LB / b2;
        double dp = dtil1_i_vE(LAp_Bb2, mub, int(p), lscfdp, thr_margin)(p);
        double lBdet = log(LB * b2).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s + 1, s + m + 1).lgamma() - lgamma(s) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * log(2) + q * log(b2) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf.minCoeff());

        return Rcpp::List::create(
            Rcpp::Named("ansseq")     = ansseq,
            Rcpp::Named("errseq")     = errseq,
            Rcpp::Named("twosided")   = twosided,
            Rcpp::Named("diminished") = diminished);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq")     = ansseq,
            Rcpp::Named("diminished") = diminished);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LA,
                     const Eigen::MatrixXd UA, const Eigen::ArrayXd LB,
                     const double b2, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, bool error_bound = true,
                     const double thr_margin = 100, int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = htil3_pjk_mE(A, Bh, zeromat, mu, m, int(p), lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2) + q * log(b2)
                              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / b2) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s = std::max(q, r);
        MatrixXd Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
        lscf.setZero();
        ArrayXXd dkstm = hhat3_pjk_mE(Ap, Bh, zeromat, mu, m, int(p), lscf, thr_margin, nthreads).row(p);
        dkstm.resize(m + 1, m + 1);
        dkstm /= exp(lscf - lscf.minCoeff());
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        set_cumsum(dkst, cumsum_dkst);
        MatrixXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        MatrixXd Ap_Bb2 = Bisqr * Ap * Bisqr / b2;
        double dp = dtil1_i_mE(Ap_Bb2, mub, int(p), lscfdp, thr_margin)(p);
        double lBdet = log(LB * b2).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s + 1, s + m + 1).lgamma() - lgamma(s) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * log(2) + q * log(b2) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf.minCoeff());

        return Rcpp::List::create(
            Rcpp::Named("ansseq")     = ansseq,
            Rcpp::Named("errseq")     = errseq,
            Rcpp::Named("twosided")   = twosided,
            Rcpp::Named("diminished") = diminished);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq")     = ansseq,
            Rcpp::Named("diminished") = diminished);
    }
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const double b1, const double b2,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const double b1, const double b2,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - b1 * A;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const double b1, const double b2, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = h3_ijk_vE(LAh, LBh, zeromat, mu, m, lscf, thr_margin); // , nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const double b1, const double b2, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - b1 * A;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = h3_ijk_mE(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cvEc(const Eigen::ArrayXd LB, const Eigen::ArrayXd LD,
                     const double b2, const double b3,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - b3 * LD;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d2_ij_vE(LBh, LDh, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cmEc(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                     const double b2, const double b3,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - b3 * D;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d2_ij_mE(Bh, Dh, m, lscf, thr_margin);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nvEc(const Eigen::ArrayXd LB, const Eigen::ArrayXd LD,
                     const double b2, const double b3, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - b3 * LD;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = h3_ijk_vE(zeromat, LBh, LDh, mu, m, lscf, thr_margin); // , nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r)  - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nmEc(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                     const double b2, const double b3, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - b3 * D;
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = h3_ijk_mE(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_cvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const Eigen::ArrayXd LD,
                     const double b2, const double b3,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - b3 * LD;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d3_pjk_vE(LA, LBh, LDh, m, int(p), lscf, thr_margin).row(p); // , nthreads).row(p);
    dks.resize(m + 1, m + 1);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_cmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const Eigen::MatrixXd D,
                     const double b2, const double b3,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - b3 * D;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = d3_pjk_mE(A, Bh, Dh, m, int(p), lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const Eigen::ArrayXd LD,
                     const double b2, const double b3,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - b3 * LD;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = htil3_pjk_vE(LA, LBh, LDh, mu, m, int(p), lscf, thr_margin).row(p); // , nthreads).row(p);
    dks.resize(m + 1, m + 1);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const Eigen::MatrixXd D,
                     const double b2, const double b3,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - b3 * D;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXXd dks = htil3_pjk_mE(A, Bh, Dh, mu, m, int(p), lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    ArrayXXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf.rowwise().reverse().matrix().diagonal().array() < 0) &&
                       dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).array()).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_cvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const Eigen::ArrayXd LD,
                     const double b1, const double b2, const double b3,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) { // , int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - b3 * LD;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = d3_ijk_vE(LAh, LBh, LDh, m, lscf, thr_margin); // , nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_cmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const Eigen::MatrixXd D,
                     const double b1, const double b2, const double b3,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - b1 * A;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - b3 * D;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = d3_ijk_mE(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const Eigen::ArrayXd LD,
                     const double b1, const double b2, const double b3,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - b1 * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - b2 * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - b3 * LD;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = h3_ijk_vE(LAh, LBh, LDh, mu, m, lscf, thr_margin); // , nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const Eigen::MatrixXd D,
                     const double b1, const double b2, const double b3,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const int m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const int n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - b1 * A;
    MatrixXd Bh = (ArrayXd::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - b3 * D;
    ArrayXXd lscf = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXXd dks = h3_ijk_mE(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = ((lscf.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).array() < 0) &&
                      dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).array()).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
