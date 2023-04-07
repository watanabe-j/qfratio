#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

#include "dk_funs.h"
#include "hgs_funs.h"

using std::lgammal;
using Eigen::log;
using Eigen::abs;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Index;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXl;
typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;
typedef Eigen::DiagonalMatrix<long double, Eigen::Dynamic> DiagMatXl;


Eigen::Array<long double, Eigen::Dynamic, 1> LinSpaced_lgammal(const Eigen::Index n, const long double par) {
    ArrayXl ans(n);
    for(Index i = 0; i < n; i++) {
        ans(i) = lgammal(par + (long double)(i));
    }
    return ans;
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, noncentral, long double
//'
// [[Rcpp::export]]
SEXP ApIq_npi_nEl(const Eigen::Array<long double, Eigen::Dynamic, 1> LA,
                  const long double bA,
                  const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                  const long double p_ = 1, const long double q_ = 1,
                  const Eigen::Index m = 100, const long double thr_margin = 100,
                  int nthreads = 1) {
    const Index n = LA.size();
    const long double n_ = n;
    ArrayXl LAh = ArrayXl::Ones(n) - bA * LA;
    ArrayXl zeromat = ArrayXl::Zero(n);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl dks = h2_ij_vE(LAh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXl ansmat = hgs_2dE(dks, -p_, q_, n_ / 2, ((p_ - q_) * M_LN2 - p_ * log(bA)
                             + lgammal(n_ / 2 + p_ - q_) - lgammal(n_ / 2)), lscf);
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_npi()}, long double
//'
// [[Rcpp::export]]
SEXP ApBq_npi_El(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                   const long double bA, const long double bB,
                 const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                 const long double p_ = 1, const long double q_ = 1, const Eigen::Index m = 100,
                 const long double thr_margin = 100, int nthreads = 0,
                 const long double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const long double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl dks((m + 1) * (m + 2) / 2);
    if(use_vec) {
        ArrayXl LAh = ArrayXl::Ones(n) - bA * A.diagonal().array();
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
        if(central) {
            dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h2_ij_vE(LAh, LBh, mu, m, lscf, thr_margin, nthreads);
}
    } else {
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
        if(central) {
            dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h2_ij_mE(Ah, Bh, mu, m, lscf, thr_margin, nthreads);
}
}
    ArrayXl ansmat = hgs_2dE(dks, -p_, q_, n_ / 2, ((p_ - q_) * M_LN2 - p_ * log(bA) + q_ * log(bB)
                              + lgammal(n_ / 2 + p_ - q_) - lgammal(n_ / 2)), lscf);
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nEl(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bB,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p_ = 1, const long double q_ = 1, const long double r_ = 1,
                     const Eigen::Index m = 100, bool error_bound = true,
                    const long double thr_margin = 100, int nthreads = 0,
                    const long double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index n = LB.size();
    const long double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    ArrayXl dks((m + 1) * (m + 2) / 2);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl LA(n);
    MatrixXl UA(n, n);
    ArrayXl LBh(n);
    DiagMatXl Bh(n);
    if(use_vec) {
        LA = A.diagonal().array();
        LBh = ArrayXl::Ones(n) - bB * LB;
        ArrayXl zeromat = ArrayXl::Zero(n);
        dks = htil3_pjk_vE(LA, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixXl> eigA(A);
        LA = eigA.eigenvalues();
        UA = eigA.eigenvectors();
        Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXl zeromat = MatrixXl::Zero(n, n);
        dks = htil3_pjk_mE(A, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    }
    // dks.resize(m + 1, m + 1);
    ArrayXl ansmat = hgs_2dE(dks, q_, r_, n_ / 2 + p_, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB)
                              + lgammal(p_ + 1) + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2 + p_)), lscf);
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXl LAp = abs(LA);
        ArrayXl mub = sqrt(3 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s_ = std::max(q_, r_);
        lscf.setZero();
        ArrayXl dkstm((m + 1) * (m + 2) / 2);
        ArrayXl lscfdp = ArrayXl::Zero(p + 1);
        double dp;
        if(use_vec) {
            ArrayXl zeromat = ArrayXl::Zero(n);
            dkstm = hhat3_pjk_vE(LAp, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
            dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
       } else {
            MatrixXl zeromat = MatrixXl::Zero(n, n);
            MatrixXl Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
            dkstm = hhat3_pjk_mE(Ap, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
            DiagMatXl Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
            dp = dtil1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        }
        // dkstm.resize(m + 1, m + 1);
        ArrayXl dkst = sum_counterdiagE(dkstm);
        ArrayXl cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        double lBdet = log(LB * bB).sum();
        ArrayXl lcoefe =
            LinSpaced_lgammal(m + 1, s_ + 1) - lgammal(s_) -
            LinSpaced_lgammal(m + 1, n_ / 2 + p_ + 1) +
            lgammal(n_ / 2 + p_ - q_ - r_);
        lcoefe += (p_ - q_ - r_) * M_LN2 + q_ * log(bB) + lgammal(p_ + 1);
        ArrayXl errseq = exp(lcoefe + (deldif2 + log(dp) - lscfdp(p) - lBdet / 2)) -
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
//'   \code{qfmrm_ApBIqr_npi()}, long double
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_El(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const long double bA, const long double bB,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p_ = 1, const long double q_ = 1, const long double r_ = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 0, const long double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const long double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl dks;
    ArrayXl ansseq(m + 1);
    if(central) {
        if(use_vec) {
            ArrayXl LAh = ArrayXl::Ones(n) - bA * A.diagonal().array();
            ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
            dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
        } else {
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
            dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
        }
    ArrayXl ansmat = hgs_2dE(dks, -p_, q_, n_ / 2, ((p_ - q_ - r_) * M_LN2 - p_ * log(bA) + q_ * log(bB)
                              + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2)), lscf);
        ansseq = sum_counterdiagE(ansmat);
    } else {
        if(use_vec) {
            ArrayXl LAh = ArrayXl::Ones(n) - bA * A.diagonal().array();
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
    ArrayXl zeromat = ArrayXl::Zero(n);
            dks = h3_ijk_vE(LAh, LBh, zeromat, mu, m, lscf, thr_margin, nthreads);
        } else {
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl zeromat = MatrixXl::Zero(n, n);
            dks = h3_ijk_mE(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
        }
    ArrayXl ansmat = hgs_3dE(dks, -p_, q_, r_, n_ / 2, ((p_ - q_ - r_) * M_LN2 - p_ * log(bA) + q_ * log(bB)
                              + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2)), lscf);
        ansseq = sum_counterdiag3DE(ansmat);
    }
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, long double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_El(const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bB, const long double bD,
                   const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p_ = 1, const long double q_ = 1, const long double r_ = 1,
                     const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 0, const long double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const long double n_ = n;
    bool use_vec = is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl dks;
    ArrayXl ansseq(m + 1);
    if(central) {
        if(use_vec) {
            ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
            ArrayXl LDh = ArrayXl::Ones(n) - bD * D.diagonal().array();
            dks = d2_ij_vE(LDh, LBh, m, lscf, thr_margin, nthreads);
        } else {
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    // Here, DiagMat Bh is the 2nd par; r & q should be used accordingly in hgs_2dE
            dks = d2_ij_mE(Dh, Bh, m, lscf, thr_margin, nthreads);
        }
    ArrayXl ansmat = hgs_2dE(dks, r_, q_, n_ / 2, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB) + r_ * log(bD)
                              + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2)), lscf);
        ansseq = sum_counterdiagE(ansmat);
    } else {
        if(use_vec) {
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
            ArrayXl LDh = ArrayXl::Ones(n) - bD * D.diagonal().array();
    ArrayXl zeromat = ArrayXl::Zero(n);
            dks = h3_ijk_vE(zeromat, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
        } else {
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
    MatrixXl zeromat = MatrixXl::Zero(n, n);
            dks = h3_ijk_mE(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
        }
    ArrayXl ansmat = hgs_3dE(dks, -p_, q_, r_, n_ / 2, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB) + r_ * log(bD)
                              + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2)), lscf);
        ansseq = sum_counterdiag3DE(ansmat);
    }
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_El(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A, const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bB, const long double bD,
                     const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                     const long double p_ = 1, const long double q_ = 1, const long double r_ = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 0, const long double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index n = LB.size();
    const long double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl dks((m + 1) * (m + 2) / 2);
    if(use_vec) {
        ArrayXl LA = A.diagonal().array();
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
        ArrayXl LDh = ArrayXl::Ones(n) - bD * D.diagonal().array();
        if(central) {
            dks = d3_pjk_vE(LA, LBh, LDh, m, p, lscf, thr_margin, nthreads).row(p);
        } else {
            dks = htil3_pjk_vE(LA, LBh, LDh, mu, m, p, lscf, thr_margin, nthreads).row(p);
}
    } else {
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
        if(central) {
            dks = d3_pjk_mE(A, Bh, Dh, m, p, lscf, thr_margin, nthreads).row(p);
        } else {
            dks = htil3_pjk_mE(A, Bh, Dh, mu, m, p, lscf, thr_margin, nthreads).row(p);
        }
    }
    // dks.resize(m + 1, m + 1);
    ArrayXl ansmat = hgs_2dE(dks, q_, r_, n_ / 2 + p_, ((p_ - q_ - r_) * M_LN2
                              + q_ * log(bB) + r_ * log(bD) + lgammal(p_ + 1)
                              + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2 + p_)), lscf);
    ArrayXl ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, long double
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_El(const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> A, const Eigen::Array<long double, Eigen::Dynamic, 1> LB,
                     const Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> D,
                     const long double bA, const long double bB, const long double bD,
                    const Eigen::Array<long double, Eigen::Dynamic, 1> mu,
                    const long double p_ = 1, const long double q_ = 1, const long double r_ = 1,
                   const Eigen::Index m = 100, const long double thr_margin = 100,
                   int nthreads = 0, const long double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const long double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXl lscf = ArrayXl::Zero(m + 1);
    ArrayXl dks((m + 1) * (m + 2) * (m + 3) / 6);
    if(use_vec) {
        ArrayXl LAh = ArrayXl::Ones(n) - bA * A.diagonal().array();
    ArrayXl LBh = ArrayXl::Ones(n) - bB * LB;
        ArrayXl LDh = ArrayXl::Ones(n) - bD * D.diagonal().array();
        if(central) {
            dks = d3_ijk_vE(LAh, LBh, LDh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h3_ijk_vE(LAh, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
}
    } else {
    MatrixXl Ah = MatrixXl::Identity(n, n) - bA * A;
    DiagMatXl Bh = (ArrayXl::Ones(n) - bB * LB).matrix().asDiagonal();
    MatrixXl Dh = MatrixXl::Identity(n, n) - bD * D;
        if(central) {
            dks = d3_ijk_mE(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h3_ijk_mE(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
        }
    }
    ArrayXl ansmat = hgs_3dE(dks, -p_, q_, r_, n_ / 2, ((p_ - q_ - r_) * M_LN2
                              - p_ * log(bA) + q_ * log(bB) + r_ * log(bD)
                              + lgammal(n_ / 2 + p_ - q_ - r_) - lgammal(n_ / 2)), lscf);
    ArrayXl ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
