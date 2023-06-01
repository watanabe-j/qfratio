#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::ArrayXd;
// using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
// using Eigen::LDLT;


// Eigen function to obtain matrix square root via Cholesky decomposition
// with pivoting (accommodates singular matrices)
inline Eigen::MatrixXd cholpiv_UE(const Eigen::MatrixXd& X) {
    Eigen::LDLT<MatrixXd> ldltX(X);
    MatrixXd D = ldltX.vectorD();
    MatrixXd U = ldltX.matrixU();
    return D.array().sqrt().matrix().asDiagonal() * U * ldltX.transpositionsP().transpose();
}

// Eigen function to obtain multivariate normal variables
// with specified mean vector and covariance
// // [[Rcpp::export]]
Eigen::MatrixXd rmvnE(const int N, const Eigen::VectorXd& mu,
                      const Eigen::MatrixXd& Sigma) {
    const int nv = Sigma.rows();

    // LDLT method (Cholesky with pivoting)
    MatrixXd cSigma = cholpiv_UE(Sigma);

    // // LLT method (Cholesky w/o pivoting) does not work for singular ones
    // MatrixXd cSigma = Sigma.llt().matrixU();

    // // Matrix square root using eigendecomposition
    // MatrixXd cSigma = matsqrtE(Sigma);

    // Rcpp::RNGScope scope; // Unnecessary when used via Rcpp attributes only
    Rcpp::NumericVector x = Rcpp::rnorm(N * nv, 0, 1);
    Eigen::Map<MatrixXd> X(x.begin(), N, nv);
    X *= cSigma;
    X.rowwise() += mu.transpose();
    return X;
}

//' @describeIn qfrm_cpp
//'   \code{rqfp()}
//'
// [[Rcpp::export]]
Eigen::ArrayXd rqfpE(const int nit,
                     const Eigen::MatrixXd A, const Eigen::MatrixXd B, const Eigen::MatrixXd D,
                     const double p_, const double q_, const double r_,
                     const Eigen::VectorXd mu, const Eigen::MatrixXd Sigma) {
    MatrixXd X = rmvnE(nit, mu, Sigma);
    ArrayXd qfAp(nit);
    ArrayXd qfBq(nit);
    ArrayXd qfDr(nit);
    if(p_ == 0) {
        qfAp.setOnes();
    } else {
        qfAp = (X * A * X.transpose()).diagonal().array().pow(p_);
    }
    if(q_ == 0) {
        qfBq.setOnes();
    } else {
        qfBq = (X * B * X.transpose()).diagonal().array().pow(q_);
    }
    if(r_ == 0) {
        qfDr.setOnes();
    } else {
        qfDr = (X * D * X.transpose()).diagonal().array().pow(r_);
    }
    ArrayXd ans = qfAp * qfBq * qfDr;
    return ans;
}
