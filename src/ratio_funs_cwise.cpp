#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

#include "dk_funs.h"
#include "dk_funs_cwise.h"
#include "hgs_funs.h"

using Eigen::log;
using Eigen::abs;
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
SEXP ApIq_npi_nEc(const Eigen::ArrayXd LA, const Eigen::MatrixXd UA, const double bA,
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
//'   \code{qfrm_ApBq_npi()}, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBq_npi_Ec(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                 const double bA, const double bB, const Eigen::ArrayXd mu,
                 const double p = 1, const double q = 1, const Eigen::Index m = 100,
                 const double thr_margin = 100, int nthreads = 0,
                 const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks((m + 1) * (m + 2) / 2);
    if(use_vec) {
        ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
        ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
        if(central) {
            dks = d2_ij_vEc(LAh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h2_ij_vEc(LAh, LBh, mu, m, lscf, thr_margin, nthreads);
        }
    } else {
        MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
        DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        if(central) {
            dks = d2_ij_mEc(Ah, Bh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h2_ij_mEc(Ah, Bh, mu, m, lscf, thr_margin, nthreads);
        }
    }
    ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q) * M_LN2 - p * log(bA) + q * log(bB)
                              + lgamma(n_ / 2 + p - q) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, noncentral, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nEc(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const double bB, const Eigen::ArrayXd mu,
                    const double p = 1, const double q = 1, const double r = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100, int nthreads = 0,
                    const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    ArrayXd dks((m + 1) * (m + 2) / 2);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd LA(n);
    MatrixXd UA(n, n);
    ArrayXd LBh(n);
    DiagMatXd Bh(n);
    if(use_vec) {
        LA = A.diagonal().array();
        LBh = ArrayXd::Ones(n) - bB * LB;
        ArrayXd zeromat = ArrayXd::Zero(n);
        dks = htil3_pjk_vEc(LA, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
        LA = eigA.eigenvalues();
        UA = eigA.eigenvectors();
        Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXd zeromat = MatrixXd::Zero(n, n);
        dks = htil3_pjk_mEc(A, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    }
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dEc(dks, q, r, n_ / 2 + p, ((p - q - r) * M_LN2 + q * log(bB)
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
        lscf.setZero();
        ArrayXd dkstm((m + 1) * (m + 2) / 2);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp;
        if(use_vec) {
            ArrayXd zeromat = ArrayXd::Zero(n);
            dkstm = hhat3_pjk_vEc(LAp, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
            dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        } else {
            MatrixXd zeromat = MatrixXd::Zero(n, n);
            MatrixXd Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
            dkstm = hhat3_pjk_mEc(Ap, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
            DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
            dp = dtil1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        }
        // dkstm.resize(m + 1, m + 1);
        ArrayXd dkst = sum_counterdiagE(dkstm);
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
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
//'   \code{qfmrm_ApBIqr_npi()}, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_Ec(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const double bA, const double bB, const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const double r = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf;
    ArrayXd dks;
    ArrayXd ansseq(m + 1);
    if(central) {
        lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
        if(use_vec) {
            ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            dks = d2_ij_vEc(LAh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            dks = d2_ij_mEc(Ah, Bh, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_2dEc(dks, -p, q, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                                + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiagE(ansmat);
    } else {
        lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
        if(use_vec) {
            ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            ArrayXd zeromat = ArrayXd::Zero(n);
            dks = h3_ijk_vEc(LAh, LBh, zeromat, mu, m, lscf, thr_margin, nthreads);
        } else {
            MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            MatrixXd zeromat = MatrixXd::Zero(n, n);
            dks = h3_ijk_mEc(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 - p * log(bA) + q * log(bB)
                                + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiag3DE(ansmat);
    }
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_Ec(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                   const double bB, const double bD, const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const double r = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf;
    ArrayXd dks;
    ArrayXd ansseq(m + 1);
    if(central) {
        lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
        if(use_vec) {
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
            dks = d2_ij_vEc(LDh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
            // Here, DiagMat Bh is the 2nd par; r & q should be used accordingly in hgs_2dE
            dks = d2_ij_mEc(Dh, Bh, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_2dEc(dks, r, q, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                                + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiagE(ansmat);
    } else {
        lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
        if(use_vec) {
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
            ArrayXd zeromat = ArrayXd::Zero(n);
            dks = h3_ijk_vEc(zeromat, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
        } else {
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
            MatrixXd zeromat = MatrixXd::Zero(n, n);
            dks = h3_ijk_mEc(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2 + q * log(bB) + r * log(bD)
                                + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiag3DE(ansmat);
    }
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_Ec(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const Eigen::MatrixXd D,
                   const double bB, const double bD,
                   const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const double r = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) / 2);
    ArrayXd dks((m + 1) * (m + 2) / 2);
    if(use_vec) {
        ArrayXd LA = A.diagonal().array();
        ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
        ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
        if(central) {
            dks = d3_pjk_vEc(LA, LBh, LDh, m, p, lscf, thr_margin, nthreads).row(p);
        } else {
            dks = htil3_pjk_vEc(LA, LBh, LDh, mu, m, p, lscf, thr_margin, nthreads).row(p);
        }
    } else {
        DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
        if(central) {
            dks = d3_pjk_mEc(A, Bh, Dh, m, p, lscf, thr_margin, nthreads).row(p);
        } else {
            dks = htil3_pjk_mEc(A, Bh, Dh, mu, m, p, lscf, thr_margin, nthreads).row(p);
        }
    }
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
//'   \code{qfmrm_ApBDqr_npi()}, coefficient-wise scaling
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_Ec(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const Eigen::MatrixXd D,
                   const double bA, const double bB, const double bD,
                   const Eigen::ArrayXd mu,
                   const double p = 1, const double q = 1, const double r = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    ArrayXd dks((m + 1) * (m + 2) * (m + 3) / 6);
    if(use_vec) {
        ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
        ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
        ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
        if(central) {
            dks = d3_ijk_vEc(LAh, LBh, LDh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h3_ijk_vEc(LAh, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
        }
    } else {
        MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
        DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
        if(central) {
            dks = d3_ijk_mEc(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h3_ijk_mEc(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
        }
    }
    ArrayXd ansmat = hgs_3dEc(dks, -p, q, r, n_ / 2, ((p - q - r) * M_LN2
                              - p * log(bA) + q * log(bB) + r * log(bD)
                              + lgamma(n_ / 2 + p - q - r) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = ((lscf < 0) && dks.cwiseEqual(0)).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
