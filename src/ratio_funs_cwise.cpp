#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

#include "dk_funs.h"
#include "dk_funs_cwise.h"
#include "hgs_funs.h"

using Eigen::log;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Index;

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatXd;


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, noncentral, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApIq_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::MatrixXd UA, const double bA,
                   const Eigen::VectorXd mu,
                   const double p = 1, const double q = 1, const Eigen::Index m = 100,
                   const double thr_margin = 100, int nthreads = 1) {
    const Index n = LA.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd mud = UA.transpose() * mu;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = h2_ij_vEc(LAh, zeromat, mud, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA)
                             + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                   const double bA, const double bB,
                   const double p = 1, const double q = 1, const Eigen::Index m = 100,
                   const double thr_margin = 100, int nthreads = 1)
{
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d2_ij_vEc(LAh, LBh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const double bA, const double bB,
                   const double p = 1, const double q = 1, const Eigen::Index m = 100, 
                   const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d2_ij_mEc(Ah, Bh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                   const double bA, const double bB, const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const Eigen::Index m = 100,
                   const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = h2_ij_vEc(LAh, LBh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const double bA, const double bB, const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const Eigen::Index m = 100,
                   const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = h2_ij_mEc(Ah, Bh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const double bB, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, bool error_bound = true,
                     const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = htil3_pjk_vEc(LA, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
                              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s = std::max(q, r);
        lscf.setZero();
        ArrayXd dkstm = hhat3_pjk_vEc(LAp, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
        // dkstm.resize(m + 1, m + 1);
        dkstm /= exp(lscf - lscf.minCoeff());
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        set_cumsum(dkst, cumsum_dkst);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s + 1, s + m + 1).lgamma() - lgamma(s) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * M_LN2 + q * log(bB) + lgamma(p + 1);
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
                     const double bB, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, bool error_bound = true,
                     const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = htil3_pjk_mEc(A, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
                              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s = std::max(q, r);
        MatrixXd Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
        lscf.setZero();
        ArrayXd dkstm = hhat3_pjk_mEc(Ap, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
        // dkstm.resize(m + 1, m + 1);
        dkstm /= exp(lscf - lscf.minCoeff());
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        set_cumsum(dkst, cumsum_dkst);
        DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = dtil1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s + 1, s + m + 1).lgamma() - lgamma(s) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * M_LN2 + q * log(bB) + lgamma(p + 1);
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
                     const double bA, const double bB,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d2_ij_vEc(LAh, LBh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const double bA, const double bB,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d2_ij_mEc(Ah, Bh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nvEc(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                     const double bA, const double bB, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = h3_ijk_vEc(LAh, LBh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nmEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const double bA, const double bB, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = h3_ijk_mEc(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cvEc(const Eigen::ArrayXd LB, const Eigen::ArrayXd LD,
                     const double bB, const double bD,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d2_ij_vEc(LBh, LDh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cmEc(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                     const double bB, const double bD,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    // Here, DiagMat Bh is the 2nd par; r & q should be used accordingly in hgs_2dE
    ArrayXd dks = d2_ij_mEc(Dh, Bh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dEc(dks, r, q, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & vector, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nvEc(const Eigen::ArrayXd LB, const Eigen::ArrayXd LD,
                     const double bB, const double bD, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = h3_ijk_vEc(zeromat, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r)  - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & matrix, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nmEc(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                     const double bB, const double bD, const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = h3_ijk_mEc(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bB, const double bD,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d3_pjk_vEc(LA, LBh, LDh, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bB, const double bD,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = d3_pjk_mEc(A, Bh, Dh, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bB, const double bD,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = htil3_pjk_vEc(LA, LBh, LDh, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bB, const double bD,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks = htil3_pjk_mEc(A, Bh, Dh, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bA, const double bB, const double bD,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = d3_ijk_vEc(LAh, LBh, LDh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bA, const double bB, const double bD,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = d3_ijk_mEc(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bA, const double bB, const double bD,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = h3_ijk_vEc(LAh, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
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
                     const double bA, const double bB, const double bD,
                     const Eigen::ArrayXd mu,
                     const double p = 1, const double q = 1, const double r = 1,
                     const Eigen::Index m = 100, const double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks = h3_ijk_mEc(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
