#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

#include "dk_funs.h"
#include "hgs_funs.h"

using std::lgammal;
using Eigen::log;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Index;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXl;
typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;


Eigen::Array<long double, Eigen::Dynamic, 1> LinSpaced_lgammal(const Eigen::Index n, const long double par) {
    long double n_ = n;
    ArrayXl ans(n);
    for(long double i = 0; i < n_; i++) {
        ans(i) = lgammal(par + i);
    }
    return ans;
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, noncentral
//'
// [[Rcpp::export]]
SEXP ApIq_npi_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                   const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> UA,
                   const long double bA,
                   const Eigen::Matrix<long double, Eigen::Dynamic, 1> mu,
                   const long double p = 1, const long double q = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 1) {
    const Index n = LA.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl zeromat = ArrayXl::Zero(n);
    ArrayXl mud = UA.transpose() * mu;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h2_ij_vE(LAh, zeromat, mud, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA)
                             + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                             // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                   const long double bA, const long double bB,
                   const long double p = 1, const long double q = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, central & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBq_npi_cmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                   const long double bA, const long double bB,
                   const long double p = 1, const long double q = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                   const long double bA, const long double bB,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                   const long double p = 1, const long double q = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h2_ij_vE(LAh, LBh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, noncentral & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBq_npi_nmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                   const long double bA, const long double bB,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                   const long double p = 1, const long double q = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h2_ij_mE(Ah, Bh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, bool error_bound = true,
                     const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_vE(LA, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXl LAp = abs(LA);
        ArrayXl mub = sqrt(3 / bB) * mu / LB.sqrt();
        long double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        long double s = std::max(q, r);
        lscf.setZero();
        ArrayXXl dkstm = hhat3_pjk_vE(LAp, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
        dkstm.resize(m + 1, m + 1);
        ArrayXl dkst = sum_counterdiagE(dkstm);
        ArrayXl cumsum_dkst(m + 1);
        // ArrayXl lscfp0 = lscf.row(p).head(m + 1);
        // dkst /= exp(lscfp0 - lscfp0(m));
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        ArrayXl lscfdp = ArrayXl::Zero(p + 1);
        long double dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        long double lBdet = log(LB * bB).sum();
        ArrayXl lcoefe =
            LinSpaced_lgammal(m + 1, s + 1) - lgammal(s) -
            LinSpaced_lgammal(m + 1, n_ / 2 + p + 1) +
            lgammal(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * M_LN2 + q * log(bB) + lgammal(p + 1);
        ArrayXl errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         // exp(lcoefe + log(cumsum_dkst)); // - lscf);
                         // exp(lcoefe + log(cumsum_dkst) - lscfp0(m));
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
//'   \code{qfmrm_ApBIqr_int()}, noncentral & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> UA,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, bool error_bound = true,
                     const long double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_mE(A, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXl LAp = abs(LA);
        ArrayXl mub = sqrt(3 / bB) * mu / LB.sqrt();
        long double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        long double s = std::max(q, r);
        MatrixXl Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
        lscf.setZero();
        ArrayXXl dkstm = hhat3_pjk_mE(Ap, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
        dkstm.resize(m + 1, m + 1);
        ArrayXl dkst = sum_counterdiagE(dkstm);
        ArrayXl cumsum_dkst(m + 1);
        // ArrayXl lscfp0 = lscf.row(p).head(m + 1);
        // dkst /= exp(lscfp0 - lscfp0(m));
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        MatrixXl Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXl lscfdp = ArrayXl::Zero(p + 1);
        long double dp = dtil1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        long double lBdet = log(LB * bB).sum();
        ArrayXl lcoefe =
            LinSpaced_lgammal(m + 1, s + 1) - lgammal(s) -
            LinSpaced_lgammal(m + 1, n_ / 2 + p + 1) +
            lgammal(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * M_LN2 + q * log(bB) + lgammal(p + 1);
        ArrayXl errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
                         // exp(lcoefe + log(cumsum_dkst)); // - lscf);
                         // exp(lcoefe + log(cumsum_dkst) - lscfp0(m));
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
//'   \code{qfmrm_ApBIqr_npi()}, central & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bA, const long double bB,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, central & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_cmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bA, const long double bB,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bA, const long double bB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_vE(LAh, LBh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_npi()}, noncentral & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_nmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bA, const long double bB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_mE(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & vector, long double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LD,
                     const long double bB, const long double bD,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                     int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - bD * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_vE(LBh, LDh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, central & matrix, long double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_cmEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bB, const long double bD,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                     int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_mE(Bh, Dh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & vector, long double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LD,
                     const long double bB, const long double bD,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - bD * LD;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_vE(zeromat, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r)  - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r)  - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, noncentral & matrix, long double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_nmEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bB, const long double bD, const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_mE(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, central & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_cvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA, const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LD,
                     const long double bB, const long double bD,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - bD * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_pjk_vE(LA, LBh, LDh, m, p, lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, central & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_cmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A, const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bB, const long double bD,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_pjk_mE(A, Bh, Dh, m, p, lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, noncentral & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LD,
                     const long double bB, const long double bD,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - bD * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_vE(LA, LBh, LDh, mu, m, p, lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, noncentral & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_nmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bB, const long double bD,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_mE(A, Bh, Dh, mu, m, p, lscf, thr_margin, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2
                              + q * log(bB) + r * log(bD) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() &&
                      dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, central & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_cvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LD,
                     const long double bA, const long double bB, const long double bD,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - bD * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_ijk_vE(LAh, LBh, LDh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, central & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_cmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bA, const long double bB, const long double bD,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_ijk_mE(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, noncentral & vector, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_nvEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA, const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                    const Eigen::Array<long double, Eigen::Dynamic, 1> LD,
                    const long double bA, const long double bB, const long double bD,
                    const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                    const long double p = 1, const long double q = 1, const long double r = 1,
                    const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 1) {
    const Index n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - bD * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_vE(LAh, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, noncentral & matrix, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_nmEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A, const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                    const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                    const long double bA, const long double bB, const long double bD,
                    const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                    const long double p = 1, const long double q = 1, const long double r = 1,
                    const Eigen::Index m = 100, const long double thr_margin = 100, int nthreads = 0) {
    const Index n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_mE(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    if((lscf < 0).any()) {
        for(Index i = 0; i <= m && !diminished; i++) {
            diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
