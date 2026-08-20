#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>

// These are to use gsl
#include "gsl/integration/gsl_integration.h"
#include "gsl/err/gsl_errno.h"

#include "dk_funs.h"
#include "hgs_funs.h"

using Eigen::exp;
using Eigen::log;
using Eigen::abs;
using Eigen::ArrayXi;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::Index;

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagMatXd;


struct bao_tB_params_m {
    const Eigen::MatrixXd *A;
    const Eigen::ArrayXd *LB;
    const Eigen::ArrayXd *mu;
    const double *p_;
    const double *q_;
    const Eigen::Index *p_i;
    const double *epsabs;
    const double *epsrel;
    const int *limit;
};

double bao_tB_fun_int_c_m(double u, void *p)
{
    struct bao_tB_params_m *params = (struct bao_tB_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    ArrayXd lscfdp = ArrayXd::Zero((*p_i) + 1);
    double d_til = d1_i_mE(Ar, *p_i, lscfdp, 100.0)(*p_i);

    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() * d_til;

    return out;
}

double bao_tB_fun_int_n_m(double u, void *p)
{
    struct bao_tB_params_m *params = (struct bao_tB_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *mu = (params->mu);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    ArrayXd mu_til = Delta * (*mu);
    ArrayXd lscfdp = ArrayXd::Zero((*p_i) + 1);
    double d_til = dtil1_i_mE(Ar, mu_til, *p_i, lscfdp, 100.0)(*p_i);

    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() *
                 exp(mu_til.matrix().squaredNorm() / 2.0) * d_til;

    return out;
}

struct bao_mr_params_m {
    const Eigen::MatrixXd *A;
    const Eigen::ArrayXd *LB;
    const Eigen::MatrixXd *D;
    const Eigen::ArrayXd *mu;
    const double *p_;
    const double *q_;
    const double *r_;
    const Eigen::Index *p_i;
    const double *epsabs;
    const double *epsrel;
    const int *limit;
};

double bao_mr_fun_int_c_m(double t, void *p)
{
    struct bao_mr_params_m *params = (struct bao_mr_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const MatrixXd *D = (params->D);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);
    
    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Dr = DeltaD * (*D) * DeltaD;
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    SelfAdjointEigenSolver<MatrixXd> eigDr(Dr);
    ArrayXd LDr = eigDr.eigenvalues();
    MatrixXd HDr = eigDr.eigenvectors();
    MatrixXd Ar_HDr = HDr.transpose() * Ar * HDr;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_m params_pass;
    params_pass.A = &Ar_HDr;
    params_pass.LB = &LDr;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    gsl_function F;
    F.function = &bao_tB_fun_int_c_m;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_mr_fun_int_n_m(double t, void *p)
{
    struct bao_mr_params_m *params = (struct bao_mr_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const MatrixXd *D = (params->D);
    const ArrayXd *mu = (params->mu);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Dr = DeltaD * (*D) * DeltaD;
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    SelfAdjointEigenSolver<MatrixXd> eigDr(Dr);
    ArrayXd LDr = eigDr.eigenvalues();
    MatrixXd HDr = eigDr.eigenvectors();
    MatrixXd Ar_HDr = HDr.transpose() * Ar * HDr;
    ArrayXd mu_til = HDr.transpose() * (Delta * (*mu)).matrix();

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_m params_pass;
    params_pass.A = &Ar_HDr;
    params_pass.LB = &LDr;
    params_pass.mu = &mu_til;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    gsl_function F;
    F.function = &bao_tB_fun_int_n_m;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}



struct bao_tB_params_v {
    const Eigen::ArrayXd *LA;
    const Eigen::ArrayXd *LB;
    const Eigen::ArrayXd *mu;
    const double *p_;
    const double *q_;
    const Eigen::Index *p_i;
    const double *epsabs;
    const double *epsrel;
    const int *limit;
};

double bao_tB_fun_int_c_v(double u, void *p)
{
    struct bao_tB_params_v *params = (struct bao_tB_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LAr = (*LA) * Delta2;

    ArrayXd lscfdp = ArrayXd::Zero((*p_i) + 1);
    double d_til = d1_i_vE(LAr, *p_i, lscfdp, 100.0)(*p_i);

    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() * d_til;

    return out;
}

double bao_tB_fun_int_n_v(double u, void *p)
{
    struct bao_tB_params_v *params = (struct bao_tB_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *mu = (params->mu);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LAr = (*LA) * Delta2;
    ArrayXd mu_til = Delta * (*mu);
    ArrayXd lscfdp = ArrayXd::Zero((*p_i) + 1);
    double d_til = dtil1_i_vE(LAr, mu_til, *p_i, lscfdp, 100.0)(*p_i);

    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() *
                 exp(mu_til.matrix().squaredNorm() / 2.0) * d_til;

    return out;
}

struct bao_mr_params_v {
    const Eigen::ArrayXd *LA;
    const Eigen::ArrayXd *LB;
    const Eigen::ArrayXd *LD;
    const Eigen::ArrayXd *mu;
    const double *p_;
    const double *q_;
    const double *r_;
    const Eigen::Index *p_i;
    const double *epsabs;
    const double *epsrel;
    const int *limit;
};

double bao_mr_fun_int_c_v(double t, void *p)
{
    struct bao_mr_params_v *params = (struct bao_mr_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *LD = (params->LD);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LDr = (*LD) * Delta2;
    ArrayXd LAr = (*LA) * Delta2;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_v params_pass;
    params_pass.LA = &LAr;
    params_pass.LB = &LDr;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    params_pass.epsabs = epsabs;
    params_pass.epsrel = epsrel;
    params_pass.limit = limit;
    gsl_function F;
    F.function = &bao_tB_fun_int_c_v;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_mr_fun_int_n_v(double t, void *p)
{
    struct bao_mr_params_v *params = (struct bao_mr_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *LD = (params->LD);
    const ArrayXd *mu = (params->mu);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LDr = (*LD) * Delta2;
    ArrayXd LAr = (*LA) * Delta2;
    ArrayXd mu_til = (Delta * (*mu));

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_v params_pass;
    params_pass.LA = &LAr;
    params_pass.LB = &LDr;
    params_pass.mu = &mu_til;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    params_pass.epsabs = epsabs;
    params_pass.epsrel = epsrel;
    params_pass.limit = limit;
    gsl_function F;
    F.function = &bao_tB_fun_int_n_v;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}



struct bao_sA_params {
    const Eigen::ArrayXd *LA;
    const Eigen::ArrayXd *mu;
    const double *p_;
    const Index *p_i;
};

double bao_sA_fun_npi_c(double s, void *p)
{
    struct bao_sA_params *params = (struct bao_sA_params *)p;
    const ArrayXd *LA = (params->LA);
    const double *p_ = (params->p_);
    const Index *p_i = (params->p_i);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * s * (*LA));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd R = (*LA) * Delta2;
    ArrayXd lscfdp = ArrayXd::Zero((*p_i) + 1);
    double d_til = d1_i_vE(R, *p_i, lscfdp, 100.0)(*p_i);

    double out = std::pow(s, (*p_) - 1.0) * Delta.prod() * d_til;

    return out;
}

double bao_sA_fun_npi_n(double s, void *p)
{
    struct bao_sA_params *params = (struct bao_sA_params *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *mu = (params->mu);
    const double *p_ = (params->p_);
    const Index *p_i = (params->p_i);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * s * (*LA));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd R = (*LA) * Delta2;
    ArrayXd mu_til = Delta * (*mu);
    ArrayXd lscfdp = ArrayXd::Zero((*p_i) + 1);
    double d_til = dtil1_i_vE(R, mu_til, *p_i, lscfdp, 100.0)(*p_i);

    double out = std::pow(s, (*p_) - 1.0) * Delta.prod() *
                 exp(mu_til.matrix().squaredNorm() / 2.0) * d_til;

    return out;
}

double bao_tB_fun_npi_c_m(double u, void *p)
{
    struct bao_tB_params_m *params = (struct bao_tB_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    SelfAdjointEigenSolver<MatrixXd> eigAr(Ar, Eigen::EigenvaluesOnly);
    ArrayXd LAr = eigAr.eigenvalues();

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_sA_params params_pass;
    params_pass.LA = &LAr;
    params_pass.p_ = p_;
    params_pass.p_i = p_i;
    gsl_function F;
    F.function = &bao_sA_fun_npi_c;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_tB_fun_npi_n_m(double u, void *p)
{
    struct bao_tB_params_m *params = (struct bao_tB_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *mu = (params->mu);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    SelfAdjointEigenSolver<MatrixXd> eigAr(Ar);
    ArrayXd LAr = eigAr.eigenvalues();
    MatrixXd HAr = eigAr.eigenvectors();
    ArrayXd mu_til = HAr.transpose() * (Delta * (*mu)).matrix();

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_sA_params params_pass;
    params_pass.LA = &LAr;
    params_pass.mu = &mu_til;
    params_pass.p_ = p_;
    params_pass.p_i = p_i;
    gsl_function F;
    F.function = &bao_sA_fun_npi_n;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_mr_fun_npi_c_m(double t, void *p)
{
    struct bao_mr_params_m *params = (struct bao_mr_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const MatrixXd *D = (params->D);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Dr = DeltaD * (*D) * DeltaD;
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    SelfAdjointEigenSolver<MatrixXd> eigDr(Dr);
    ArrayXd LDr = eigDr.eigenvalues();
    MatrixXd HDr = eigDr.eigenvectors();
    MatrixXd Ar_HDr = HDr.transpose() * Ar * HDr;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_m params_pass;
    params_pass.A = &Ar_HDr;
    params_pass.LB = &LDr;
    params_pass.p_ = p_;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    params_pass.epsabs = epsabs;
    params_pass.epsrel = epsrel;
    params_pass.limit = limit;
    gsl_function F;
    F.function = &bao_tB_fun_npi_c_m;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_mr_fun_npi_n_m(double t, void *p)
{
    struct bao_mr_params_m *params = (struct bao_mr_params_m *)p;
    const MatrixXd *A = (params->A);
    const ArrayXd *LB = (params->LB);
    const MatrixXd *D = (params->D);
    const ArrayXd *mu = (params->mu);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    DiagMatXd DeltaD = Delta.matrix().asDiagonal();
    MatrixXd Dr = DeltaD * (*D) * DeltaD;
    MatrixXd Ar = DeltaD * (*A) * DeltaD;
    SelfAdjointEigenSolver<MatrixXd> eigDr(Dr);
    ArrayXd LDr = eigDr.eigenvalues();
    MatrixXd HDr = eigDr.eigenvectors();
    MatrixXd Ar_HDr = HDr.transpose() * Ar * HDr;
    ArrayXd mu_til = HDr.transpose() * (Delta * (*mu)).matrix();

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_m params_pass;
    params_pass.A = &Ar_HDr;
    params_pass.LB = &LDr;
    params_pass.mu = &mu_til;
    params_pass.p_ = p_;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    params_pass.epsabs = epsabs;
    params_pass.epsrel = epsrel;
    params_pass.limit = limit;
    gsl_function F;
    F.function = &bao_tB_fun_npi_n_m;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}


double bao_tB_fun_npi_c_v(double u, void *p)
{
    struct bao_tB_params_v *params = (struct bao_tB_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LAr = (*LA) * Delta2;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_sA_params params_pass;
    params_pass.LA = &LAr;
    params_pass.p_ = p_;
    params_pass.p_i = p_i;
    gsl_function F;
    F.function = &bao_sA_fun_npi_c;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_tB_fun_npi_n_v(double u, void *p)
{
    struct bao_tB_params_v *params = (struct bao_tB_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *mu = (params->mu);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * u * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LAr = (*LA) * Delta2;
    ArrayXd mu_til = Delta * (*mu);

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_sA_params params_pass;
    params_pass.LA = &LAr;
    params_pass.mu = &mu_til;
    params_pass.p_ = p_;
    params_pass.p_i = p_i;
    gsl_function F;
    F.function = &bao_sA_fun_npi_n;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(u, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_mr_fun_npi_c_v(double t, void *p)
{
    struct bao_mr_params_v *params = (struct bao_mr_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *LD = (params->LD);
    const double *p_ = (params->p_);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const Index *p_i = (params->p_i);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LDr = (*LD) * Delta2;
    ArrayXd LAr = (*LA) * Delta2;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_v params_pass;
    params_pass.LA = &LAr;
    params_pass.LB = &LDr;
    params_pass.p_ = p_;
    params_pass.q_ = r_;
    params_pass.p_i = p_i;
    params_pass.epsabs = epsabs;
    params_pass.epsrel = epsrel;
    params_pass.limit = limit;
    gsl_function F;
    F.function = &bao_tB_fun_npi_c_v;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}

double bao_mr_fun_npi_n_v(double t, void *p)
{
    struct bao_mr_params_v *params = (struct bao_mr_params_v *)p;
    const ArrayXd *LA = (params->LA);
    const ArrayXd *LB = (params->LB);
    const ArrayXd *LD = (params->LD);
    const ArrayXd *mu = (params->mu);
    const double *p_ = (params->p_);
    const Index *p_i = (params->p_i);
    const double *q_ = (params->q_);
    const double *r_ = (params->r_);
    const double *epsabs = (params->epsabs);
    const double *epsrel = (params->epsrel);
    const int *limit = (params->limit);

    ArrayXd Delta2 = 1.0 / (1.0 + 2.0 * t * (*LB));
    ArrayXd Delta = Delta2.sqrt();
    ArrayXd LDr = (*LD) * Delta2;
    ArrayXd LAr = (*LA) * Delta2;
    ArrayXd mu_til = (Delta * (*mu));

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(*limit);
    double value, error;
    int status;
    struct bao_tB_params_v params_pass;
    params_pass.LA = &LAr;
    params_pass.LB = &LDr;
    params_pass.mu = &mu_til;
    params_pass.p_ = p_;
    params_pass.p_i = p_i;
    params_pass.q_ = r_;
    params_pass.epsabs = epsabs;
    params_pass.epsrel = epsrel;
    params_pass.limit = limit;
    gsl_function F;
    F.function = &bao_tB_fun_npi_n_v;
    F.params = &params_pass;
    status = gsl_integration_qagiu(&F, 0, *epsabs, *epsrel, *limit, w,
                                   &value, &error);
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        // if (stop_on_error)
        //     Rcpp::stop(errmsg);
        // else
            Rcpp::warning(errmsg);
    }
    double out = std::pow(t, (*q_) - 1.0) * Delta.prod() * value;

    return out;
}



//' @describeIn qfrm_cpp
//'   \code{qfmrm_integ_int()}, double
//'
// [[Rcpp::export]]
SEXP integ_mr_int_Ed(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const Eigen::MatrixXd D, const Eigen::ArrayXd mu,
                     const double p_, const double q_, const double r_,
                     bool stop_on_error, const double tol_zero,
                     double epsabs, double epsrel, int limit)
{
    const Index p_i = p_;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    double cons = exp(-mu.matrix().squaredNorm() / 2.0 + p_ * M_LN2
                      + lgamma(p_ + 1.0) - lgamma(q_) - lgamma(r_));
    epsabs /= cons;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);
    double value, error;
    int status;
    gsl_function F;
    if (use_vec) {
        ArrayXd LA = A.diagonal();
        ArrayXd LD = D.diagonal();
        struct bao_mr_params_v params;
        params.LA = &LA;
        params.LB = &LB;
        params.LD = &LD;
        params.mu = &mu;
        params.p_ = &p_;
        params.q_ = &q_;
        params.r_ = &r_;
        params.p_i = &p_i;
        params.epsabs = &epsabs;
        params.epsrel = &epsrel;
        params.limit = &limit;
        F.params = &params;
        if (central) {
            F.function = &bao_mr_fun_int_c_v;
        }
        else {
            F.function = &bao_mr_fun_int_n_v;
        }
        // Function call must be within the same scope with params;
        // otherwise segfault occurs
        status = gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, w,
                                       &value, &error);
    }
    else {
        struct bao_mr_params_m params;
        params.A = &A;
        params.LB = &LB;
        params.D = &D;
        params.mu = &mu;
        params.p_ = &p_;
        params.q_ = &q_;
        params.r_ = &r_;
        params.p_i = &p_i;
        params.epsabs = &epsabs;
        params.epsrel = &epsrel;
        params.limit = &limit;
        F.params = &params;
        if (central) {
            F.function = &bao_mr_fun_int_c_m;
        }
        else {
            F.function = &bao_mr_fun_int_n_m;
        }
        status = gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, w,
                                       &value, &error);
    }
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        if (stop_on_error)
            Rcpp::stop(errmsg);
        else
            Rcpp::warning(errmsg);
    }

    value *= cons;
    error *= cons;

    return Rcpp::List::create(
        Rcpp::Named("value")     = value,
        Rcpp::Named("abs.error") = error);
}

//' @describeIn qfrm_cpp
//'   \code{qfmrm_integ_npi()}, double
//'
// [[Rcpp::export]]
SEXP integ_mr_npi_Ed(const Eigen::MatrixXd A, const Eigen::ArrayXd LB,
                     const Eigen::MatrixXd D, const Eigen::ArrayXd mu,
                     const double p_, const double q_, const double r_,
                     bool stop_on_error, const double tol_zero,
                     double epsabs, double epsrel, int limit)
{
    const double p_c = ceil(p_);
    const double p_r = p_c - p_;
    const Index p_i = p_c;
    bool use_vec = is_diag_E(A, tol_zero) && is_diag_E(D, tol_zero);
    bool central = is_zero_E(mu, tol_zero);
    double cons = exp(-mu.matrix().squaredNorm() / 2.0 + p_c * M_LN2
                       + lgamma(p_c + 1.0)
                       - lgamma(p_r) - lgamma(q_) - lgamma(r_));
    epsabs /= cons;

    gsl_set_error_handler_off();
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(limit);
    double value, error;
    int status;
    gsl_function F;
    if (use_vec) {
        ArrayXd LA = A.diagonal();
        ArrayXd LD = D.diagonal();
        struct bao_mr_params_v params;
        params.LA = &LA;
        params.LB = &LB;
        params.LD = &LD;
        params.mu = &mu;
        params.p_ = &p_r;
        params.q_ = &q_;
        params.r_ = &r_;
        params.p_i = &p_i;
        params.epsabs = &epsabs;
        params.epsrel = &epsrel;
        params.limit = &limit;
        F.params = &params;
        if (central) {
            F.function = &bao_mr_fun_npi_c_v;
        }
        else {
            F.function = &bao_mr_fun_npi_n_v;
        }
        status = gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, w,
                                       &value, &error);
    }
    else {
        struct bao_mr_params_m params;
        params.A = &A;
        params.LB = &LB;
        params.D = &D;
        params.mu = &mu;
        params.p_ = &p_r;
        params.q_ = &q_;
        params.r_ = &r_;
        params.p_i = &p_i;
        params.epsabs = &epsabs;
        params.epsrel = &epsrel;
        params.limit = &limit;
        F.params = &params;
        if (central) {
            F.function = &bao_mr_fun_npi_c_m;
        }
        else {
            F.function = &bao_mr_fun_npi_n_m;
        }
        status = gsl_integration_qagiu(&F, 0, epsabs, epsrel, limit, w,
                                       &value, &error);
    }
    gsl_integration_workspace_free(w);
    if (status) {
        std::string errmsg = "problem in gsl_integration_qagiu():\n  ";
        errmsg += gsl_strerror(status);
        if (stop_on_error)
            Rcpp::stop(errmsg);
        else
            Rcpp::warning(errmsg);
    }

    value *= cons;
    error *= cons;

    return Rcpp::List::create(
        Rcpp::Named("value")     = value,
        Rcpp::Named("abs.error") = error);
}
