#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

using Eigen::ArrayXd;
// using Eigen::ArrayXXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
// using Eigen::LDLT;


inline Eigen::MatrixXd cholpiv_UE(const Eigen::MatrixXd& X) {
    Eigen::LDLT<MatrixXd> ldltX(X);
    MatrixXd D = ldltX.vectorD();
    MatrixXd U = ldltX.matrixU();
    return D.array().sqrt().matrix().asDiagonal() * U * ldltX.transpositionsP();
}

// inline Eigen::MatrixXd cholpiv_LE(const Eigen::MatrixXd& X) {
//     Eigen::LDLT<MatrixXd> ldltX(X);
//     MatrixXd D = ldltX.vectorD();
//     MatrixXd L = ldltX.matrixL();
//     return ldltX.transpositionsP().transpose() * L * D.array().sqrt().matrix().asDiagonal();
// }

// inline Eigen::MatrixXd matsqrtE(const Eigen::MatrixXd& X) {
//     Eigen::SelfAdjointEigenSolver<MatrixXd> eigX(X);
//     return eigX.operatorSqrt();
// }

// [[Rcpp::export]]
Eigen::MatrixXd rmvnE(const int N, const Eigen::VectorXd& mu,
                      const Eigen::MatrixXd& Sigma) {
    const int nv = Sigma.rows();

    // LDLT method (Cholesky with pivoting)
    MatrixXd cSigma = cholpiv_UE(Sigma);

    // // LLT method (Cholesky w/o pivoting) does not work for singular ones
    // MatrixXd cSigma = Sigma.llt().matrixU();

    // // Matrix square root using eigendecomposition
    // MatrixXd cSigma = matsqrtE(Sigma);

    Rcpp::NumericVector x = Rcpp::rnorm(N * nv, 0, 1);
    Eigen::Map<MatrixXd> X(x.begin(), N, nv);
    X *= cSigma;
    X.rowwise() += mu.transpose();
    return X;
}

// [[Rcpp::export]]
Eigen::ArrayXd rqfrE(const int nit,
                     const Eigen::MatrixXd A, const Eigen::MatrixXd B,
                     const double p, const double q,
                     const Eigen::VectorXd mu, const Eigen::MatrixXd Sigma) {
    MatrixXd X = rmvnE(nit, mu, Sigma);
    // MatrixXd cAL = cholpiv_LE(A);
    // MatrixXd cBL = cholpiv_LE(B);
    // ArrayXd Num = (X * cAL).rowwise().squaredNorm().array().pow(p);
    // ArrayXd Den = (X * cBL).rowwise().squaredNorm().array().pow(q);
    ArrayXd Num = (X * A * X.transpose()).diagonal().array().pow(p);
    ArrayXd Den = (X * B * X.transpose()).diagonal().array().pow(q);
    ArrayXd ans = Num / Den;
    return ans;
}

// [[Rcpp::export]]
Eigen::ArrayXd rqfmrE(const int nit,
                     const Eigen::MatrixXd A, const Eigen::MatrixXd B, const Eigen::MatrixXd D,
                     const double p, const double q, const double r,
                     const Eigen::VectorXd mu, const Eigen::MatrixXd Sigma) {
    MatrixXd X = rmvnE(nit, mu, Sigma);
    ArrayXd Num = (X * A * X.transpose()).diagonal().array().pow(p);
    ArrayXd Den1 = (X * B * X.transpose()).diagonal().array().pow(q);
    ArrayXd Den2 = (X * D * X.transpose()).diagonal().array().pow(r);
    ArrayXd ans = Num / Den1 / Den2;
    return ans;
}

// [[Rcpp::export]]
Eigen::ArrayXd rqfpE(const int nit,
                     const Eigen::MatrixXd A, const Eigen::MatrixXd B, const Eigen::MatrixXd D,
                     const double p, const double q, const double r,
                     const Eigen::VectorXd mu, const Eigen::MatrixXd Sigma) {
    MatrixXd X = rmvnE(nit, mu, Sigma);
    ArrayXd Num1 = (X * A * X.transpose()).diagonal().array().pow(p);
    ArrayXd Num2 = (X * B * X.transpose()).diagonal().array().pow(q);
    ArrayXd Num3 = (X * D * X.transpose()).diagonal().array().pow(r);
    ArrayXd ans = Num1 * Num2 * Num3;
    return ans;
}
