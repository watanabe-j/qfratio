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
using Eigen::abs;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Index;

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatXd;



//' @describeIn qfrm_cpp
//'   \code{qfm_Ap_int()}
//'
// [[Rcpp::export]]
SEXP Ap_int_E(const Eigen::MatrixXd A, const Eigen::ArrayXd mu,
              const double p_ = 1, const double thr_margin = 100,
              const double tol_zero = 2.2e-14) {
    const Index p = p_;
    ArrayXd lscf = ArrayXd::Zero(p + 1);
    double dp;
    if(is_zero_E(mu, tol_zero)) {
        dp = d1_i_mE(A, p, lscf, thr_margin)(p);
    } else {
        dp = dtil1_i_mE(A, mu, p, lscf, thr_margin)(p);
    }
    double ans = exp(p_ * M_LN2 + lgamma(p_ + 1) - lscf(p)) * dp;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}


//' @describeIn qfrm_cpp
//'   \code{qfpm_ABpq_int()}
//'
// [[Rcpp::export]]
SEXP ABpq_int_E(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                const Eigen::ArrayXd mu,
                const double p_ = 1, const double q_ = 1,
                const double thr_margin = 100,
                const double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index q = q_;
    ArrayXd lscf = ArrayXd::Zero(q + 1);
    double dpq;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    if(use_vec) {
        ArrayXd LA = A.diagonal().array();
        if(central) {
            dpq = d2_pj_vE(LA, LB, q, p, lscf, thr_margin)(p, q);
        } else {
            dpq = dtil2_pq_vE(LA, LB, mu, p, q)(p, q);
        }
    } else {
        DiagMatXd B = LB.matrix().asDiagonal();
        if(central) {
            dpq = d2_pj_mE(A, B, q, p, lscf, thr_margin)(p, q);
        } else {
            dpq = dtil2_pq_mE(A, B, mu, p, q)(p, q);
        }
    }
    double ans = exp((p_ + q_) * M_LN2 + lgamma(p_ + 1) + lgamma(q_ + 1) - lscf(q)) * dpq;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}


//' @describeIn qfrm_cpp
//'   \code{qfpm_ABDpqr_int()}
//'
// [[Rcpp::export]]
SEXP ABDpqr_int_E(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                  const Eigen::MatrixXd D, const Eigen::ArrayXd mu,
                  const double p_ = 1, const double q_ = 1, const double r_ = 1,
                  const double thr_margin = 100,
                  const double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index q = q_;
    const Index r = r_;
    ArrayXd lscf = ArrayXd::Zero(q + r + 1);
    double dpqr;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    if(use_vec) {
        ArrayXd LA = A.diagonal().array();
        ArrayXd LD = D.diagonal().array();
        if(central) {
            dpqr = d3_pjk_vE(LA, LB, LD, q + r, p, lscf, thr_margin, 1)(p, (2 * q + r + 3) * r / 2 + q);
        } else {
            dpqr = dtil3_pqr_vE(LA, LB, LD, mu, p, q, r)(p, (q + 1) * (r + 1) - 1);
        }
    } else {
        DiagMatXd B = LB.matrix().asDiagonal();
        if(central) {
            dpqr = d3_pjk_mE(A, B, D, q + r, p, lscf, thr_margin, 1)(p, (2 * q + r + 3) * r / 2 + q);
        } else {
            dpqr = dtil3_pqr_mE(A, B, D, mu, p, q, r)(p, (q + 1) * (r + 1) - 1);
        }
    }
    double ans = exp((p_ + q_ + r_) * M_LN2 + lgamma(p_ + 1) + lgamma(q_ + 1) + lgamma(r_ + 1) - lscf(q + r)) * dpqr;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_int()}, central
//'
// [[Rcpp::export]]
SEXP ApIq_int_cE(const Eigen::MatrixXd A,
                 const double p_ = 1, const double q_ = 1, 
                 const double thr_margin = 100) {
    const Index p = p_;
    const Index n = A.rows();
    const double n_ = n;
    ArrayXd lscf = ArrayXd::Zero(p + 1);
    double dp = d1_i_mE(A, p, lscf, thr_margin)(p);
    double ans = exp((p_ - q_) * M_LN2 + lgamma(p_ + 1) + lgamma(n_ / 2 + p_ - q_) -
                     lgamma(n_ / 2 + p_) - lscf(p)) * dp;
    return Rcpp::List::create(Rcpp::Named("ans") = ans);
}

//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_int()}, noncentral
//'
// [[Rcpp::export]]
SEXP ApIq_int_nE(const Eigen::MatrixXd A, const Eigen::ArrayXd mu,
                 const double p_ = 1, const double q_ = 1, 
                 const double thr_margin = 100) {
    const Index p = p_;
    const Index n = A.rows();
    const double n_ = n;
    ArrayXd aps = a1_pk_mE(A, mu, p, thr_margin).row(p);
    ArrayXd ls = ArrayXd::LinSpaced(p + 1, 0, p_);
    double nsqnorm2 = -mu.matrix().squaredNorm() / 2;

    // Calling gsl::hyperg_1F1 from R. This is slow
    Rcpp::Function hyperg_1F1 = Rcpp::Environment::namespace_env("gsl")["hyperg_1F1"];
    Rcpp::NumericVector hgv = hyperg_1F1(q_, n_ / 2 + p_ + ls, nsqnorm2);
    Eigen::Map<ArrayXd> hgres(hgv.begin(), p + 1);

    // // Using the C++ library GSL with RcppGSL; this is ideal but
    // // not quite portable as the library need to be installed separately.
    // ArrayXd hgres(p + 1);
    // for(Index i = 0; i <= p; i++) {
    //     hgres(i) = gsl_sf_hyperg_1F1(q, n_ / 2 + p + ls(i), nsqnorm2);
    // }

    ArrayXd ansseq =
        exp((p_ - q_) * M_LN2 + lgamma(p_ + 1)
            + (ls + n_ / 2 + p_ - q_).lgamma() - M_LN2 * ls
            - (ls + 1).lgamma() - (ls + n_ / 2 + p_).lgamma()) * hgres * aps;
    return Rcpp::List::create(Rcpp::Named("ansseq") = ansseq);
}


//' @describeIn qfrm_cpp
//'   \code{qfrm_ApIq_npi()}, central
//'
// [[Rcpp::export]]
SEXP ApIq_npi_cE(const Eigen::ArrayXd LA,
                 const double bA,
                 const double p_ = 1, const double q_ = 1,
                 const Eigen::Index m = 100, bool error_bound = true,
                 const double thr_margin = 100) {
    const Index n = LA.size();
    const double n_ = n;
    const double m_ = m;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = d1_i_vE(LAh, m, lscf, thr_margin);
    ArrayXd ansseq = hgs_1dE(dks, -p_, n_ / 2, ((p_ - q_) * M_LN2 - p_ * log(bA)
                             + lgamma(n_ / 2 + p_ - q_) - lgamma(n_ / 2)), lscf);

    if(error_bound) {
        ArrayXd cumsum_dkst(m + 1);
        ArrayXd signseq = get_sign_rfp1(-p_, m + 1);
        dks /= exp(lscf - lscf(m));
        set_cumsum(dks, cumsum_dkst);
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, -p_ + 1, -p_ + 1 + m_).lgamma() - lgamma(-p_) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + 1, n_ / 2 + 1 + m_).lgamma() +
            lgamma(n_ / 2 + p_ - q_);
        lcoefe += (p_ - q_) * M_LN2 - p_ * log(bA);
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
//'   \code{qfrm_ApIq_npi()}, noncentral, double
//'
// [[Rcpp::export]]
SEXP ApIq_npi_nEd(const Eigen::ArrayXd LA, const double bA,
                  const Eigen::ArrayXd mu,
                  const double p_ = 1, const double q_ = 1, const Eigen::Index m = 100,
                  const double thr_margin = 100, int nthreads = 1) {
    const Index n = LA.size();
    const double n_ = n;
    ArrayXd LAh = ArrayXd::Ones(n) - bA * LA;
    ArrayXd zeromat = ArrayXd::Zero(n);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks = h2_ij_vE(LAh, zeromat, mu, m, lscf, thr_margin, nthreads);
    ArrayXd ansmat = hgs_2dE(dks, -p_, q_, n_ / 2, ((p_ - q_) * M_LN2 - p_ * log(bA)
                             + lgamma(n_ / 2 + p_ - q_) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfrm_ApBq_int()}
//'
// [[Rcpp::export]]
SEXP ApBq_int_E(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                const double bB, const Eigen::ArrayXd mu,
                const double p_ = 1, const double q_ = 1,
                const Eigen::Index m = 100, bool error_bound = true, 
                const double thr_margin = 100,
                const double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index n = LB.size();
    const double n_ = n;
    const double m_ = m;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd dks(m + 1);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd LA(n);
    MatrixXd UA(n, n);
    ArrayXd LBh(n);
    DiagMatXd Bh(n);
    if(use_vec) {
        LA = A.diagonal().array();
        LBh = ArrayXd::Ones(n) - bB * LB;
        if(central) {
            dks = d2_pj_vE(LA, LBh, m, p, lscf, thr_margin).row(p);
        } else {
            dks = htil2_pj_vE(LA, LBh, mu, m, p, lscf, thr_margin).row(p);
        }
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
        LA = eigA.eigenvalues();
        UA = eigA.eigenvectors();
        Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        if(central) {
            dks = d2_pj_mE(A, Bh, m, p, lscf, thr_margin).row(p);
        } else {
            dks = htil2_pj_mE(A, Bh, mu, m, p, lscf, thr_margin).row(p);
        }
    }
    ArrayXd ansseq = hgs_1dE(dks, q_, n_ / 2 + p_, ((p_ - q_) * M_LN2 + q_ * log(bB)
                             + lgamma(p_ + 1) + lgamma(n_ / 2 + p_ - q_) - lgamma(n_ / 2 + p_)), lscf);

    if(error_bound) {
        bool twosided = !central || ((LA < 0).any() && ((p % 1) == 1));
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(2 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        MatrixXd Ap(n, n);
        ArrayXd dkst(m + 1);
        if(twosided) {
            lscf.setZero();
            if(use_vec) {
                if(central) {
                    dkst = d2_pj_vE(LAp, LBh, m, p, lscf, thr_margin).row(p);
                } else {
                    dkst = hhat2_pj_vE(LAp, LBh, mu, m, p, lscf, thr_margin).row(p);
                }
            } else {
                Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
                if(central) {
                    dkst = d2_pj_mE(Ap, Bh, m, p, lscf, thr_margin).row(p);
                } else {
                    dkst = hhat2_pj_mE(Ap, Bh, mu, m, p, lscf, thr_margin).row(p);
                }
            }
        } else {
            Ap = A;
            dkst = dks;
        }
        double dp;
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        if(use_vec) {
            if(central) {
                dp = d1_i_vE((LAp / LB / bB).eval(), p, lscfdp, thr_margin)(p);
            } else {
                dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
            }
        } else {
            DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
            if(central) {
                dp = d1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), p, lscfdp, thr_margin)(p);
            } else {
                dp = dtil1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), mub, p, lscfdp, thr_margin)(p);
            }
        }
        ArrayXd cumsum_dkst(m + 1);
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, q_ + 1, q_ + m_ + 1).lgamma() - lgamma(q_) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p_ + 1, n_ / 2 + p_ + m_ + 1).lgamma() +
            lgamma(n_ / 2 + p_ - q_);
        lcoefe += (p_ - q_) * M_LN2 + q_ * log(bB) + lgamma(p_ + 1);
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
//'   \code{qfrm_ApBq_npi()}, double
//'
// [[Rcpp::export]]
SEXP ApBq_npi_Ed(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                 const double bA, const double bB, const Eigen::ArrayXd mu,
                 const double p_ = 1, const double q_ = 1, const Eigen::Index m = 100,
                 const double thr_margin = 100, int nthreads = 0,
                 const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks((m + 1) * (m + 2) / 2);
    if(use_vec) {
        ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
        ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
        if(central) {
            dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h2_ij_vE(LAh, LBh, mu, m, lscf, thr_margin, nthreads);
        }
    } else {
        MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
        DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        if(central) {
            dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h2_ij_mE(Ah, Bh, mu, m, lscf, thr_margin, nthreads);
        }
    }
    ArrayXd ansmat = hgs_2dE(dks, -p_, q_, n_ / 2, ((p_ - q_) * M_LN2 - p_ * log(bA) + q_ * log(bB)
                              + lgamma(n_ / 2 + p_ - q_) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}



//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBIqr_int()}, central
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_cEd(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const double bB,
                    const double p_ = 1, const double q_ = 1, const double r_ = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100,
                    const double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index n = LB.size();
    const double n_ = n;
    const double m_ = m;
    bool use_vec = is_diag_E(A, tol_zero);
    ArrayXd dks(m + 1);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd LA(n);
    MatrixXd UA(n, n);
    ArrayXd LBh(n);
    DiagMatXd Bh(n);
    if(use_vec) {
        LA = A.diagonal().array();
        LBh = ArrayXd::Ones(n) - bB * LB;
        dks = d2_pj_vE(LA, LBh, m, p, lscf, thr_margin).row(p);
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
        LA = eigA.eigenvalues();
        UA = eigA.eigenvectors();
        Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        dks = d2_pj_mE(A, Bh, m, p, lscf, thr_margin).row(p);
    }
    ArrayXd ansseq = hgs_1dE(dks, q_, n_ / 2 + p_, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB)
              + lgamma(p_ + 1) + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2 + p_)), lscf);

    if(error_bound) {
        bool twosided = (LA < 0).any() && ((p % 1) == 1);
        ArrayXd LAp = abs(LA);
        double deldif2 = 0;
        double s_ = q_;
        ArrayXd dkst(m + 1);
        ArrayXd cumsum_dkst(m + 1);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp;
        if(use_vec) {
            if(twosided) {
                lscf.setZero();
                dkst = d2_pj_vE(LAp, LBh, m, p, lscf, thr_margin).row(p);
            } else {
                dkst = dks;
            }
            dp = d1_i_vE((LAp / LB / bB).eval(), p, lscfdp, thr_margin)(p);
        } else {
            MatrixXd Ap(n, n);
            if(twosided) {
                Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
                lscf.setZero();
                dkst = d2_pj_mE(Ap, Bh, m, p, lscf, thr_margin).row(p);
            } else {
                Ap = A;
                dkst = dks;
            }
            DiagMatXd Bisqr = LB.sqrt().matrix().asDiagonal().inverse();
            dp = d1_i_mE((Bisqr * Ap * Bisqr / bB).eval(), p, lscfdp, thr_margin)(p);
        }
        dkst /= exp(lscf - lscf(m));
        set_cumsum(dkst, cumsum_dkst);
        double lBdet = log(LB * bB).sum();
        ArrayXd lcoefe =
            ArrayXd::LinSpaced(m + 1, s_ + 1, s_ + m_ + 1).lgamma() - lgamma(s_) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p_ + 1, n_ / 2 + p_ + m_ + 1).lgamma() +
            lgamma(n_ / 2 + p_ - q_ - r_);
        lcoefe += (p_ - q_ - r_) * M_LN2 + q_ * log(bB) + lgamma(p_ + 1);
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
//'   \code{qfmrm_ApBIqr_int()}, noncentral, double
//'
// [[Rcpp::export]]
SEXP ApBIqr_int_nEd(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                    const double bB, const Eigen::ArrayXd mu,
                    const double p_ = 1, const double q_ = 1, const double r_ = 1,
                    const Eigen::Index m = 100, bool error_bound = true,
                    const double thr_margin = 100, int nthreads = 0,
                    const double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index n = LB.size();
    const double n_ = n;
    const double m_ = m;
    bool use_vec = is_diag_E(A, tol_zero);
    ArrayXd dks((m + 1) * (m + 2) / 2);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd LA(n);
    MatrixXd UA(n, n);
    ArrayXd LBh(n);
    DiagMatXd Bh(n);
    if(use_vec) {
        LA = A.diagonal().array();
        LBh = ArrayXd::Ones(n) - bB * LB;
        ArrayXd zeromat = ArrayXd::Zero(n);
        dks = htil3_pjk_vE(LA, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    } else {
        Eigen::SelfAdjointEigenSolver<MatrixXd> eigA(A);
        LA = eigA.eigenvalues();
        UA = eigA.eigenvectors();
        Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXd zeromat = MatrixXd::Zero(n, n);
        dks = htil3_pjk_mE(A, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
    }
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q_, r_, n_ / 2 + p_, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB)
                              + lgamma(p_ + 1) + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2 + p_)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();

    if(error_bound) {
        bool twosided = true;
        ArrayXd LAp = abs(LA);
        ArrayXd mub = sqrt(3 / bB) * mu / LB.sqrt();
        double deldif2 = (mub.matrix().squaredNorm() - mu.matrix().squaredNorm()) / 2;
        double s_ = std::max(q_, r_);
        lscf.setZero();
        ArrayXd dkstm((m + 1) * (m + 2) / 2);
        ArrayXd lscfdp = ArrayXd::Zero(p + 1);
        double dp;
        if(use_vec) {
            ArrayXd zeromat = ArrayXd::Zero(n);
            dkstm = hhat3_pjk_vE(LAp, LBh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
            dp = dtil1_i_vE((LAp / LB / bB).eval(), mub, p, lscfdp, thr_margin)(p);
        } else {
            MatrixXd zeromat = MatrixXd::Zero(n, n);
            MatrixXd Ap = UA * LAp.matrix().asDiagonal() * UA.transpose();
            dkstm = hhat3_pjk_mE(Ap, Bh, zeromat, mu, m, p, lscf, thr_margin, nthreads).row(p);
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
            ArrayXd::LinSpaced(m + 1, s_ + 1, s_ + m_ + 1).lgamma() - lgamma(s_) -
            ArrayXd::LinSpaced(m + 1, n_ / 2 + p_ + 1, n_ / 2 + p_ + m_ + 1).lgamma() +
            lgamma(n_ / 2 + p_ - q_ - r_);
        lcoefe += (p_ - q_ - r_) * M_LN2 + q_ * log(bB) + lgamma(p_ + 1);
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
//'   \code{qfmrm_ApBIqr_npi()}, double
//'
// [[Rcpp::export]]
SEXP ApBIqr_npi_Ed(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const double bA, const double bB, const Eigen::ArrayXd mu,
                   const double p_ = 1, const double q_ = 1, const double r_ = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks;
    ArrayXd ansseq(m + 1);
    if(central) {
        if(use_vec) {
            ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            dks = d2_ij_vE(LAh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            dks = d2_ij_mE(Ah, Bh, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_2dE(dks, -p_, q_, n_ / 2, ((p_ - q_ - r_) * M_LN2 - p_ * log(bA) + q_ * log(bB)
                                + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiagE(ansmat);
    } else {
        if(use_vec) {
            ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            ArrayXd zeromat = ArrayXd::Zero(n);
            dks = h3_ijk_vE(LAh, LBh, zeromat, mu, m, lscf, thr_margin, nthreads);
        } else {
            MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            MatrixXd zeromat = MatrixXd::Zero(n, n);
            dks = h3_ijk_mE(Ah, Bh, zeromat, mu, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_3dE(dks, -p_, q_, r_, n_ / 2, ((p_ - q_ - r_) * M_LN2 - p_ * log(bA) + q_ * log(bB)
                                + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiag3DE(ansmat);
    }
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_IpBDqr_gen()}, double
//'
// [[Rcpp::export]]
SEXP IpBDqr_gen_Ed(const Eigen::ArrayXd LB, const Eigen::MatrixXd D,
                   const double bB, const double bD, const Eigen::ArrayXd mu,
                   const double p_ = 1, const double q_ = 1, const double r_ = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks;
    ArrayXd ansseq(m + 1);
    if(central) {
        if(use_vec) {
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
            dks = d2_ij_vE(LDh, LBh, m, lscf, thr_margin, nthreads);
        } else {
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
            // Here, DiagMat Bh is the 2nd par; r & q should be used accordingly in hgs_2dE
            dks = d2_ij_mE(Dh, Bh, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_2dE(dks, r_, q_, n_ / 2, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB) + r_ * log(bD)
                                + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiagE(ansmat);
    } else {
        if(use_vec) {
            ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
            ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
            ArrayXd zeromat = ArrayXd::Zero(n);
            dks = h3_ijk_vE(zeromat, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
        } else {
            DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
            MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
            MatrixXd zeromat = MatrixXd::Zero(n, n);
            dks = h3_ijk_mE(zeromat, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
        }
        ArrayXd ansmat = hgs_3dE(dks, -p_, q_, r_, n_ / 2, ((p_ - q_ - r_) * M_LN2 + q_ * log(bB) + r_ * log(bD)
                                + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2)), lscf);
        ansseq = sum_counterdiag3DE(ansmat);
    }
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_int()}, double
//'
// [[Rcpp::export]]
SEXP ApBDqr_int_Ed(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const Eigen::MatrixXd D,
                   const double bB, const double bD,
                   const Eigen::ArrayXd mu,
                   const double p_ = 1, const double q_ = 1, const double r_ = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index p = p_;
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks((m + 1) * (m + 2) / 2);
    if(use_vec) {
        ArrayXd LA = A.diagonal().array();
        ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
        ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
        if(central) {
            dks = d3_pjk_vE(LA, LBh, LDh, m, p, lscf, thr_margin, nthreads).row(p);
        } else {
            dks = htil3_pjk_vE(LA, LBh, LDh, mu, m, p, lscf, thr_margin, nthreads).row(p);
        }
    } else {
        DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
        if(central) {
            dks = d3_pjk_mE(A, Bh, Dh, m, p, lscf, thr_margin, nthreads).row(p);
        } else {
            dks = htil3_pjk_mE(A, Bh, Dh, mu, m, p, lscf, thr_margin, nthreads).row(p);
        }
    }
    // dks.resize(m + 1, m + 1);
    ArrayXd ansmat = hgs_2dE(dks, q_, r_, n_ / 2 + p_, ((p_ - q_ - r_) * M_LN2
                              + q_ * log(bB) + r_ * log(bD) + lgamma(p_ + 1)
                              + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2 + p_)), lscf);
    ArrayXd ansseq = sum_counterdiagE(ansmat);
    // ansseq /= exp(lscf);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}


//' @describeIn qfrm_cpp
//'   \code{qfmrm_ApBDqr_npi()}, double
//'
// [[Rcpp::export]]
SEXP ApBDqr_npi_Ed(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                   const Eigen::MatrixXd D,
                   const double bA, const double bB, const double bD,
                   const Eigen::ArrayXd mu,
                   const double p_ = 1, const double q_ = 1, const double r_ = 1,
                   const Eigen::Index m = 100, const double thr_margin = 100,
                   int nthreads = 0, const double tol_zero = 2.2e-14) {
    const Index n = LB.size();
    const double n_ = n;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    ArrayXd lscf = ArrayXd::Zero(m + 1);
    ArrayXd dks((m + 1) * (m + 2) * (m + 3) / 6);
    if(use_vec) {
        ArrayXd LAh = ArrayXd::Ones(n) - bA * A.diagonal().array();
        ArrayXd LBh = ArrayXd::Ones(n) - bB * LB;
        ArrayXd LDh = ArrayXd::Ones(n) - bD * D.diagonal().array();
        if(central) {
            dks = d3_ijk_vE(LAh, LBh, LDh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h3_ijk_vE(LAh, LBh, LDh, mu, m, lscf, thr_margin, nthreads);
        }
    } else {
        MatrixXd Ah = MatrixXd::Identity(n, n) - bA * A;
        DiagMatXd Bh = (ArrayXd::Ones(n) - bB * LB).matrix().asDiagonal();
        MatrixXd Dh = MatrixXd::Identity(n, n) - bD * D;
        if(central) {
            dks = d3_ijk_mE(Ah, Bh, Dh, m, lscf, thr_margin, nthreads);
        } else {
            dks = h3_ijk_mE(Ah, Bh, Dh, mu, m, lscf, thr_margin, nthreads);
        }
    }
    ArrayXd ansmat = hgs_3dE(dks, -p_, q_, r_, n_ / 2, ((p_ - q_ - r_) * M_LN2
                              - p_ * log(bA) + q_ * log(bB) + r_ * log(bD)
                              + lgamma(n_ / 2 + p_ - q_ - r_) - lgamma(n_ / 2)), lscf);
    ArrayXd ansseq = sum_counterdiag3DE(ansmat);
    bool diminished = (lscf < 0).any() && dks.cwiseEqual(0).any();
    return Rcpp::List::create(
        Rcpp::Named("ansseq")     = ansseq,
        Rcpp::Named("diminished") = diminished);
}
