#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

#include "dk_funs.h"
#include "hgs_funs.h"

using std::lgammal;
using Eigen::log;
using Eigen::SelfAdjointEigenSolver;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXl;
typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;


// Eigen version of \code{sum_counterdiag()}
Eigen::Array<long double, Eigen::Dynamic, 1> sum_counterdiagE(const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& X) {
    const int n = X.rows();
    ArrayXl ans = ArrayXl::Zero(n);
    long double x;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j <= i; j++) {
            x = X(i - j, j);
            if(!std::isnan(x)) ans(i) += x;
        }
    }
    return ans;
}

// Eigen version of \code{sum_counterdiag3D()}
// X is a wide ArrayXXl, n * (n * n)
Eigen::Array<long double, Eigen::Dynamic, 1> sum_counterdiag3DE(const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& X) {
    const int n = X.rows();
    ArrayXl ans = ArrayXl::Zero(n);
    long double x;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j <= i; j++) {
            for(int k = 0; k <= (i - j); k++){
                x = X(i - j - k, j + k * n);
                if(!std::isnan(x)) ans(i) += x;
            }
        }
    }
    return ans;
}

Eigen::Array<long double, Eigen::Dynamic, 1> LinSpaced_lgammal(const int n, const long double par) {
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
                   const long double b1,
                   const Eigen::Matrix<long double, Eigen::Dynamic, 1> mu,
                   const long double p = 1, const long double q = 1,
                   const int m = 100, bool error_bound = false) {
    const int n = LA.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl zeromat = ArrayXl::Zero(n);
    ArrayXl mud = UA.transpose() * mu;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h2_ij_vE(LAh, zeromat, mud, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1)
                             + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                             // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                   const long double b1, const long double b2,
                   const long double p = 1, const long double q = 1,
                   const int m = 100, bool error_bound = false) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_vE(LAh, LBh, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                   const long double b1, const long double b2,
                   const long double p = 1, const long double q = 1,
                   const int m = 100, bool error_bound = false) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - b1 * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_mE(Ah, Bh,  m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                   const long double b1, const long double b2,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                   const long double p = 1, const long double q = 1,
                   const int m = 100, bool error_bound = false) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h2_ij_vE(LAh, LBh, mu, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                   const long double b1, const long double b2,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                   const long double p = 1, const long double q = 1,
                   const int m = 100, bool error_bound = false) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - b1 * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h2_ij_mE(Ah, Bh, mu, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b2,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, bool error_bound = true) { // , int nthreads = 1) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_vE(LA, LBh, zeromat, mu, m, int(p), lscf).row(p); // , nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2) + q * log(b2)
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXl LAp = abs(LA);
        ArrayXl mub = sqrt(3 / b2) * mu / LB.sqrt();
        long double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        long double s = std::max(q, r);
        lscf.setZero();
        ArrayXXl dkstm = hhat3_pjk_vE(LAp, LBh, zeromat, mu, m, int(p), lscf).row(p); // , nthreads).row(p);
        dkstm.resize(m + 1, m + 1);
        ArrayXl dkst = sum_counterdiagE(dkstm);
        ArrayXl cumsum_dkst(m + 1);
        // ArrayXl lscfp0 = lscf.row(p).head(m + 1);
        // dkst /= exp(lscfp0 - lscfp0(m));
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        ArrayXl lscfdp = ArrayXl::Zero(p + 1);
        ArrayXl LAp_Bb2 = LAp / LB / b2;
        long double dp = dtil1_i_vE(LAp_Bb2, mub, int(p), lscfdp)(p);
        long double lBdet = log(LB * b2).sum();
        ArrayXl lcoefe =
            LinSpaced_lgammal(m + 1, s + 1) - lgammal(s) -
            LinSpaced_lgammal(m + 1, n_ / 2 + p + 1) +
            lgammal(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * log(2) + q * log(b2) + lgammal(p + 1);
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
                     const long double b2,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, bool error_bound = true, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_mE(A, Bh, zeromat, mu, m, int(p), lscf, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2) + q * log(b2)
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(p + 1) + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXl LAp = abs(LA);
        ArrayXl mub = sqrt(3 / b2) * mu / LB.sqrt();
        long double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        long double s = std::max(q, r);
        MatrixXl Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
        lscf.setZero();
        ArrayXXl dkstm = hhat3_pjk_mE(Ap, Bh, zeromat, mu, m, int(p), lscf, nthreads).row(p);
        dkstm.resize(m + 1, m + 1);
        ArrayXl dkst = sum_counterdiagE(dkstm);
        ArrayXl cumsum_dkst(m + 1);
        // ArrayXl lscfp0 = lscf.row(p).head(m + 1);
        // dkst /= exp(lscfp0 - lscfp0(m));
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        MatrixXl Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
        ArrayXl lscfdp = ArrayXl::Zero(p + 1);
        MatrixXl Ap_Bb2 = Bisqr * Ap * Bisqr / b2;
        long double dp = dtil1_i_mE(Ap_Bb2, mub, int(p), lscfdp)(p);
        long double lBdet = log(LB * b2).sum();
        ArrayXl lcoefe =
            LinSpaced_lgammal(m + 1, s + 1) - lgammal(s) -
            LinSpaced_lgammal(m + 1, n_ / 2 + p + 1) +
            lgammal(n_ / 2 + p - q - r);
        lcoefe += (p - q - r) * log(2) + q * log(b2) + lgammal(p + 1);
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
                     const long double b1, const long double b2,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_vE(LAh, LBh, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b1, const long double b2,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - b1 * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_mE(Ah, Bh, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, -p, q, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b1, const long double b2,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_vE(LAh, LBh, zeromat, mu, m, lscf); // , nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                     const long double b1, const long double b2,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - b1 * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_mE(Ah, Bh, zeromat, mu, m, lscf, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) - p * log(b1) + q * log(b2)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - b3 * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_vE(LBh, LDh, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - b3 * D;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d2_ij_mE(Bh, Dh, m, lscf);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - b3 * LD;
    ArrayXl zeromat = ArrayXl::Zero(n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_vE(zeromat, LBh, LDh, mu, m, lscf); // , nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r)  - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r)  - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                     const long double b2, const long double b3, const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - b3 * D;
    MatrixXl zeromat = MatrixXl::Zero(n, n);
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_mE(zeromat, Bh, Dh, mu, m, lscf, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - b3 * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_pjk_vE(LA, LBh, LDh, m, int(p), lscf).row(p); // , nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - b3 * D;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_pjk_mE(A, Bh, Dh, m, int(p), lscf, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - b3 * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_vE(LA, LBh, LDh, mu, m, int(p), lscf).row(p); // , nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b2, const long double b3,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - b3 * D;
    // ArrayXXl lscf = ArrayXXl::Zero(p + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = htil3_pjk_mE(A, Bh, Dh, mu, m, int(p), lscf, nthreads).row(p);
    dks.resize(m + 1, m + 1);
    // ArrayXXl lscfp = lscf.row(p);
    // lscfp.resize(m + 1, m + 1);
    ArrayXXl ansmat = hgs_2dE(dks, q, r, n_ / 2 + p, ((p - q - r) * log(2)
                              + q * log(b2) + r * log(b3) + lgammal(p + 1)
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscfp);
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2 + p)));
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = dks.rowwise().reverse().matrix().diagonal().cwiseEqual(0).any();
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
                     const long double b1, const long double b2, const long double b3,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100) { // , int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - b3 * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_ijk_vE(LAh, LBh, LDh, m, lscf); // , nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                     const long double b1, const long double b2, const long double b3,
                     const long double p = 1, const long double q = 1, const long double r = 1,
                     const int m = 100, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - b1 * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - b3 * D;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = d3_ijk_mE(Ah, Bh, Dh, m, lscf, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                    const long double b1, const long double b2, const long double b3,
                    const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                    const long double p = 1, const long double q = 1, const long double r = 1,
                    const int m = 100) { // , int nthreads = 1) {
    const int n = LB.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - b1 * LA;
    ArrayXl LBh = ArrayXl::Ones(n) - b2 * LB;
    ArrayXl LDh = ArrayXl::Ones(n) - b3 * LD;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_vE(LAh, LBh, LDh, mu, m, lscf); // , nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
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
                    const long double b1, const long double b2, const long double b3,
                    const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                    const long double p = 1, const long double q = 1, const long double r = 1,
                    const int m = 100, int nthreads = 0) {
    const int n = LB.size();
    const long double n_ = n;
    MatrixXl Ah = MatrixXl::Identity(n, n) - b1 * A;
    MatrixXl Bh = (ArrayXl::Ones(n) - b2 * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - b3 * D;
    // ArrayXXl lscf = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXXl dks = h3_ijk_mE(Ah, Bh, Dh, mu, m, lscf, nthreads);
    ArrayXXl ansmat = hgs_3dE(dks, -p, q, r, n_ / 2, ((p - q - r) * log(2)
                              - p * log(b1) + q * log(b2) + r * log(b3)
                              + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)), lscf);
                              // + lgammal(n_ / 2 + p - q - r) - lgammal(n_ / 2)));
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = false;
    for(int i = 0; i <= m && !diminished; i++) {
        diminished = dks.block(0, i * (m + 1), m + 1, m + 1).rowwise().reverse().matrix().diagonal(i).cwiseEqual(0).any();
    }
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
