#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

// // These are to use gsl in ApIq_int_nmE
// #include <RcppGSL.h>
// // [[Rcpp::depends(RcppGSL)]]
// #include <gsl/gsl_sf_hyperg.h>

#include "dk_funs.h"
#include "hgs_funs.h"

using Eigen::log;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Index;

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatXd;


//' @describeIn qfrm_cpp
//'   \code{qfm_Ap_int()}, central
//'
// [[Rcpp::export]]
SEXP Ap_int_cmE(const Eigen::MatrixXd A, const double p = 1, 
                const double thr_margin = 100) {
    ArrayXd lscf = ArrayXd::Zero(p + 1);
    double dp = d1_i_mE(A, p, lscf, thr_margin)(p);
    double ans = exp(p * M_LN2 + lgamma(p + 1) - lscf(p)) * dp;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfm_Ap_int()}, noncentral
//'
// [[Rcpp::export]]
SEXP Ap_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd mu,
                const double p = 1, const double thr_margin = 100) {
    ArrayXd lscf = ArrayXd::Zero(p + 1);
    double dp = dtil1_i_mE(A, mu, p, lscf, thr_margin)(p);
    double ans = exp(p * M_LN2 + lgamma(p + 1) - lscf(p)) * dp;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}


//' @describeIn qfrm_cpp
//'   \code{qfpm_ABpq_int()}, central & vector
//'
// [[Rcpp::export]]
SEXP ABpq_int_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const double p = 1, const double q = 1, const double thr_margin = 100) {
    ArrayXd lscf = ArrayXd::Zero(q + 1);
    double dpq = d2_pj_vE(LA, LB, q, p, lscf, thr_margin)(p, q);
    double ans = exp((p + q) * M_LN2 + lgamma(p + 1) + lgamma(q + 1) - lscf(q)) * dpq;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfpm_ABpq_int()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ABpq_int_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const double p = 1, const double q = 1, const double thr_margin = 100) {
    DiagMatXd B = LB.matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(q + 1);
    double dpq = d2_pj_mE(A, B, q, p, lscf, thr_margin)(p, q);
    double ans = exp((p + q) * M_LN2 + lgamma(p + 1) + lgamma(q + 1) - lscf(q)) * dpq;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfpm_ABpq_int()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ABpq_int_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1) {
    double dpq = dtil2_pq_vE(LA, LB, mu, p, q)(p, q);
    double ans = exp((p + q) * M_LN2 + lgamma(p + 1) + lgamma(q + 1)) * dpq;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfpm_ABpq_int()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ABpq_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1) {
    DiagMatXd B = LB.matrix().asDiagonal();
    double dpq = dtil2_pq_mE(A, B, mu, p, q)(p, q);
    double ans = exp((p + q) * M_LN2 + lgamma(p + 1) + lgamma(q + 1)) * dpq;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}


//' @describeIn qfrm_cpp
//'   \code{qfpm_ABDpqr_int()}, central & vector
//'
// [[Rcpp::export]]
SEXP ABDpqr_int_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd LD,
                    const double p = 1, const double q = 1, const double r = 1, 
                    const double thr_margin = 100) {
    ArrayXd lscf = ArrayXd::Zero(q + r + 1);
    double dpqr = d3_pjk_vE(LA, LB, LD, q + r, p, lscf, thr_margin, 1)(p, (q + r + 1) * r + q);
    double ans = exp((p + q + r) * M_LN2 + lgamma(p + 1) + lgamma(q + 1) + lgamma(r + 1) - lscf(q + r)) * dpqr;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfpm_ABDpqr_int()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ABDpqr_int_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const Eigen::MatrixXd D,
                    const double p = 1, const double q = 1, const double r = 1,
                    const double thr_margin = 100) {
    DiagMatXd B = LB.matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(q + r + 1);
    double dpqr = d3_pjk_mE(A, B, D, q + r, p, lscf, thr_margin, 1)(p, (q + r + 1) * r + q);
    double ans = exp((p + q + r) * M_LN2 + lgamma(p + 1) + lgamma(q + 1) + lgamma(r + 1) - lscf(q + r)) * dpqr;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfpm_ABDpqr_int()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ABDpqr_int_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd LD, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1) {
    double dpqr = dtil3_pqr_vE(LA, LB, LD, mu, p, q, r)(p, (q + 1) * (r + 1) - 1);
    double ans = exp((p + q + r) * M_LN2 + lgamma(p + 1) + lgamma(q + 1) + lgamma(r + 1)) * dpqr;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfpm_ABDpqr_int()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ABDpqr_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const Eigen::MatrixXd D, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1) {
    DiagMatXd B = LB.matrix().asDiagonal();
    double dpqr = dtil3_pqr_mE(A, B, D, mu, p, q, r)(p, (q + 1) * (r + 1) - 1);
    double ans = exp((p + q + r) * M_LN2 + lgamma(p + 1) + lgamma(q + 1) + lgamma(r + 1)) * dpqr;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_int()}, central
//'
// [[Rcpp::export]]
SEXP ApIq_int_cmE(const Eigen::MatrixXd A,
                  const double p = 1, const double q = 1, 
                  const double thr_margin = 100) {
    const Index n = A.rows();
    const double n_ = n;
    ArrayXd lscf = ArrayXd::Zero(p + 1);
    double dp = d1_i_mE(A, p, lscf, thr_margin)(p);
    double ans = exp((p - q) * M_LN2 + lgamma(p + 1) + lgamma(n_ / 2 + p - q) -
                     lgamma(n_ / 2 + p) - lscf(p)) * dp;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_int()}, noncentral
//'
// [[Rcpp::export]]
SEXP ApIq_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd mu,
                  const double p = 1, const double q = 1, 
                  const double thr_margin = 100) {
    const Index n = A.rows();
    const double n_ = n;
    ArrayXd aps = a1_pk_mE(A, mu, p, thr_margin).row(p);
    ArrayXd ls = ArrayXd::LinSpaced(p + 1, 0, p);
    double nsqnorm2 = -mu.matrix().squaredNorm() / 2;

    // Calling gsl::hyperg_1F1 from R. This is slow
    Rcpp::Function hyperg_1F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_1F1"];
    Rcpp::NumericVector hgv = hyperg_1F1(q, n_ / 2 + p + ls, nsqnorm2);
    Eigen::Map<ArrayXd> hgres(hgv.begin(), p + 1);

    // // Using the C++ library GSL with RcppGSL; this is ideal but
    // // not quite portable as the library need to be installed separately.
    // ArrayXd hgres(p + 1);
    // for(Index i = 0; i <= p; i++) {
    //     hgres(i) = gsl_sf_hyperg_1F1(q, n_ / 2 + p + ls(i), nsqnorm2);
    // }

    ArrayXd ansseq =
        exp((p - q) * M_LN2 + lgamma(p + 1)
            + (ls + n_ / 2 + p - q).lgamma() - M_LN2 * ls
            - (ls + 1).lgamma() - (ls + n_ / 2 + p).lgamma()) * hgres * aps;
    return Rcpp::List::create(Rcpp::Named("ansseq") = ansseq);
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, central
//'
// [[Rcpp::export]]
SEXP ApIq_npi_cvE(const Eigen::ArrayXd LA,
                    const double bA,
                    const double p = 1, const double q = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100) {
    const Index n = LA.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d1_i_vE(LAh, m, lscf, thr_margin);
    ArrayXd ansseq = hgs_1dE(dks, -p, n_ / 2, ((p - q) * M_LN2 - p * log(bA)
                             + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);

    if(error_bound) {
        ArrayXd cumsum_dkst(m + 1);
        ArrayXd signseq = get_sign_rfp1(-p, m + 1);
        dks /= exp(lscf - lscf(m));
        set_cumsum(dks, cumsum_dkst);
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, -p + 1, -p + 1 + m).lgamma() - lgamma(-p) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + 1, n_ / 2 + 1 + m).lgamma() +
            lgamma(n_ / 2 + p - q);
        lcoefe += (p - q) * M_LN2 - p * log(bA);
        ArrayXd errseq = exp(lcoefe - log(bA * LA).sum() / 2) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));
        errseq *= signseq;

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, noncentral
//'
// [[Rcpp::export]]
SEXP ApIq_npi_nvE(const Eigen::ArrayXd LA, const Eigen::MatrixXd UA, const double bA,
                  const Eigen::VectorXd mu,
                  const double p = 1, const double q = 1, const Eigen::Index m = 100,
                  const double thr_margin = 100, int nthreads = 1) {
    const Index n = LA.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd mud = UA.transpose() * mu;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h2_ij_vE(LAh, zeromat, mud, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA)
                             + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_int()}, central & vector
//'
// [[Rcpp::export]]
SEXP ApBq_int_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                  const double bB,
                  const double p = 1, const double q = 1,
                  const Eigen::Index m = 100, bool error_bound = true, const double thr_margin = 100) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_pj_vE(LA, LBh, m, p, lscf, thr_margin).row(p);
    ArrayXd ansseq = hgs_1dE(dks, q, n_ / 2 + p, ((p - q) * M_LN2 + q * log(bB)
                             + lgamma(p + 1) + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2 + p)), lscf);

    if(error_bound) {
        bool twosided = (LA < 0).any() && ((int(p) % 1) == 1);
        ArrayXd LAp = abs(LA);
        double deldif2 = 0;
        ArrayXd dkst(m + 1);
        if(twosided) {
            lscf.setZero();
            dkst = d2_pj_vE(LAp, LBh, m, p, lscf, thr_margin).row(p);
        } else {
            dkst = dks;
        }
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = d1_i_vE((LAp / LB / bB).eval(), p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, q + 1, q + m + 1).lgamma() - lgamma(q) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q);
        lcoefe += (p - q) * M_LN2 + q * log(bB) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq,
            Rcpp::Named("twosided") = twosided);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_int()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ApBq_int_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LA,
                  const Eigen::MatrixXd UA, const Eigen::ArrayXd LB,
                  const double bB,
                  const double p = 1, const double q = 1,
                  const Eigen::Index m = 100, bool error_bound = true, const double thr_margin = 100) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_pj_mE(A, Bh, m, p, lscf, thr_margin).row(p);
    ArrayXd ansseq = hgs_1dE(dks, q, n_ / 2 + p, ((p - q) * M_LN2 + q * log(bB)
                             + lgamma(p + 1) + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2 + p)), lscf);

    if(error_bound) {
        bool twosided = (LA < 0).any() && ((int(p) % 1) == 1);
        ArrayXd LAp = abs(LA);
        double deldif2 = 0;
        MatrixXd Ap(n, n);
        ArrayXd dkst(m + 1);
        if(twosided) {
            Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
            lscf.setZero();
            dkst = d2_pj_mE(Ap, Bh, m, p, lscf, thr_margin).row(p);
        } else {
            Ap = A;
            dkst = dks;
        }
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = d1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, q + 1, q + m + 1).lgamma() - lgamma(q) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q);
        lcoefe += (p - q) * M_LN2 + q * log(bB) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq,
            Rcpp::Named("twosided") = twosided);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_int()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ApBq_int_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                  const double bB, const Eigen::ArrayXd mu,
                  const double p = 1, const double q = 1,
                  const Eigen::Index m = 100, bool error_bound = true, 
                  const double thr_margin = 100) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = htil2_pj_vE(LA, LBh, mu, m, p, lscf, thr_margin).row(p);
    ArrayXd lscfp = lscf.row(p);
    ArrayXd ansseq = hgs_1dE(dks, q, n_ / 2 + p, ((p - q) * M_LN2 + q * log(bB)
                             + lgamma(p + 1) + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2 + p)), lscf);

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(2 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        lscf.setZero();
        ArrayXd dkst = hhat2_pj_vE(LAp, LBh, mu, m, p, lscf, thr_margin).row(p);
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, q + 1, q + m + 1).lgamma() - lgamma(q) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q);
        lcoefe += (p - q) * M_LN2 + q * log(bB) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq,
            Rcpp::Named("twosided") = twosided);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_int()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ApBq_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LA,
                  const Eigen::MatrixXd UA, const Eigen::ArrayXd LB,
                  const double bB, const Eigen::ArrayXd mu,
                  const double p = 1, const double q = 1,
                  const Eigen::Index m = 100, bool error_bound = true, 
                  const double thr_margin = 100) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = htil2_pj_mE(A, Bh, mu, m, p, lscf, thr_margin).row(p);
    ArrayXd lscfp = lscf.row(p);
    ArrayXd ansseq = hgs_1dE(dks, q, n_ / 2 + p, ((p - q) * M_LN2 + q * log(bB)
                             + lgamma(p + 1) + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2 + p)), lscf);

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(2 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        MatrixXd Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
        lscf.setZero();
        ArrayXd dkst = hhat2_pj_mE(Ap, Bh, mu, m, p, lscf, thr_margin).row(p);
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = dtil1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, q + 1, q + m + 1).lgamma() - lgamma(q) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q);
        lcoefe += (p - q) * M_LN2 + q * log(bB) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq,
            Rcpp::Named("twosided") = twosided);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}



//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & vector
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                  const double bA, const double bB,
                  const double p = 1, const double q = 1, const Eigen::Index m = 100, 
                  const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                  const double bA, const double bB,
                  const double p = 1, const double q = 1, const Eigen::Index m = 100,
                  const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                  const double bA, const double bB, const Eigen::ArrayXd mu,
                  const double p = 1, const double q = 1, const Eigen::Index m = 100,
                  const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h2_ij_vE(LAh, LBh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                  const double bA, const double bB, const Eigen::ArrayXd mu,
                  const double p = 1, const double q = 1, const Eigen::Index m = 100,
                  const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h2_ij_mE(Ah, Bh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, central & vector
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const double bB,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_pj_vE(LA, LBh, m, p, lscf, thr_margin).row(p);
    ArrayXd ansseq = hgs_1dE(dks, q, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);

    if(error_bound) {
        bool twosided = (LA < 0).any() && ((int(p) % 1) == 1);
        ArrayXd LAp = abs(LA);
        double deldif2 = 0;
        double s = q;
        ArrayXd dkst(m + 1);
        ArrayXd cumsum_dkst(m + 1);
        if(twosided) {
            lscf.setZero();
            dkst = d2_pj_vE(LAp, LBh, m, p, lscf, thr_margin).row(p);
        } else {
            dkst = dks;
        }
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = d1_i_vE((LAp / LB / bB).eval(), p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s + 1, s + m + 1).lgamma() - lgamma(s) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * M_LN2 + q * log(bB) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq,
            Rcpp::Named("twosided") = twosided);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LA,
                    const Eigen::MatrixXd UA, const Eigen::ArrayXd LB,
                    const double bB,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_pj_mE(A, Bh, m, p, lscf, thr_margin).row(p);
    ArrayXd ansseq = hgs_1dE(dks, q, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);

    if(error_bound) {
        bool twosided = (LA < 0).any() && ((int(p) % 1) == 1);
        ArrayXd LAp = abs(LA);
        double deldif2 = 0;
        double s = q;
        ArrayXd dkst(m + 1);
        ArrayXd cumsum_dkst(m + 1);
        MatrixXd Ap(n, n);
        if(twosided) {
            Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
            lscf.setZero();
            dkst = d2_pj_mE(Ap, Bh, m, p, lscf, thr_margin).row(p);
        } else {
            Ap = A;
            dkst = dks;
        }
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp = d1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), p, lscfdp, thr_margin)(p);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s + 1, s + m + 1).lgamma() - lgamma(s) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p + 1, n_ / 2 + p + m + 1).lgamma() +
            lgamma(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * M_LN2 + q * log(bB) + lgamma(p + 1);
        ArrayXd errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq,
            Rcpp::Named("errseq") = errseq,
            Rcpp::Named("twosided") = twosided);
    } else {
        return Rcpp::List::create(
            Rcpp::Named("ansseq") = ansseq);
    }
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const double bB, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, bool error_bound = true, 
                    const double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = htil3_pjk_vE(LA, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
                              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s = std::max(q, r);
        lscf.setZero();
        ArrayXd dkstm = hhat3_pjk_vE(LAp, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
        // dkstm.resize(m + 1, m + 1);
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
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
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

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
//'   \code{qfmrm_ApBIqr_int()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LA,
                    const Eigen::MatrixXd UA, const Eigen::ArrayXd LB,
                    const double bB, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = htil3_pjk_mE(A, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
                              + lgamma(p + 1) + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s = std::max(q, r);
        MatrixXd Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
        lscf.setZero();
        ArrayXd dkstm = hhat3_pjk_mE(Ap, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
        // dkstm.resize(m + 1, m + 1);
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
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
                         exp(lcoefe + log(cumsum_dkst) - lscf(m));

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
//'   \code{qfmrm_ApBIqr_npi()}, central & vector
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const double bA, const double bB,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const double bA, const double bB,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const double bA, const double bB, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h3_ijk_vE(LAh, LBh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const double bA, const double bB, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h3_ijk_mE(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & vector
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cvE(const Eigen::ArrayXd LB, const Eigen::ArrayXd LD,
                    const double bB, const double bD,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d2_ij_vE(LBh, LDh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & matrix
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cmE(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                    const double bB, const double bD,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    // Here, DiagMat Bh is the 2nd par; r & q should be used accordingly in hgs_2dE
    ArrayXd dks = d2_ij_mE(Dh, Bh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, r, q, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nvE(const Eigen::ArrayXd LB, const Eigen::ArrayXd LD,
                    const double bB, const double bD, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h3_ijk_vE(zeromat, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r)  - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nmE(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                    const double bB, const double bD, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    MatrixXd zeromat = MatrixXd::Zero(n, n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h3_ijk_mE(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, central & vector
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd LD,
                    const double bB, const double bD,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d3_pjk_vE(LA, LBh, LDh, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const Eigen::MatrixXd D,
                    const double bB, const double bD,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d3_pjk_mE(A, Bh, Dh, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
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
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = htil3_pjk_vE(LA, LBh, LDh, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
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
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = htil3_pjk_mE(A, Bh, Dh, mu, m, p, lscf, thr_margin, nthreads).row(p);
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgamma(p + 1)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2 + p)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, central & vector
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_cvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd LD,
                    const double bA, const double bB, const double bD,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d3_ijk_vE(LAh, LBh, LDh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, central & matrix
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_cmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
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
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d3_ijk_mE(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, noncentral & vector
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_nvE(const Eigen::ArrayXd LA, const Eigen::ArrayXd LB,
                    const Eigen::ArrayXd LD,
                    const double bA, const double bB, const double bD,
                    const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, const double thr_margin = 100,
                    int nthreads = 1) {
    const Index n = LB.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
    ArrayXd LDh = ArrayXd::Ones(n) - bD * LD;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h3_ijk_vE(LAh, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, noncentral & matrix
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_nmE(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
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
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h3_ijk_mE(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
