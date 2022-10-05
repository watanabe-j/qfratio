#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <unsupported/Eigen/SpecialFunctions>
#include <cmath>
#include "hgs_funs.h"

using Eigen::log;
using Eigen::abs;
using Eigen::sign;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;
using Eigen::Matrix;
using Eigen::Array;
using Eigen::DenseBase;
using Eigen::MatrixBase;
using Eigen::ArrayBase;
using Eigen::Dynamic;

typedef Eigen::Array<long double, Eigen::Dynamic, 1> ArrayXl;
typedef Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic> ArrayXXl;
typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;

// Eigen function template to put cumulative sum of Data into Out
template <typename Derived>
void set_cumsum(const Eigen::DenseBase<Derived>& Data, Eigen::DenseBase<Derived>& Out) {
    std::partial_sum(Data.derived().data(), Data.derived().data() + Data.size(), Out.derived().data());
}
// Template instantiations, so that these be used in other places via header
template void set_cumsum<ArrayXd>(const Eigen::DenseBase<ArrayXd>& Data, Eigen::DenseBase<ArrayXd>& Out);
template void set_cumsum<ArrayXl>(const Eigen::DenseBase<ArrayXl>& Data, Eigen::DenseBase<ArrayXl>& Out);

// // Eigen function to calculate a series of log (of absolute value of)
// // rising factorial (a)_k from parameter a and length n (k = 0, 1, ... n - 1).
// // When a is a negative integer, the result will be NaN from the (a+1)-th term
// // (as it should be)
// template <typename T>
// Eigen::Array<T, Dynamic, 1> get_lrf(const T a, const int n) {
//     typedef Eigen::Array<T, Dynamic, 1> ArrayXx;
//     ArrayXx ans(n);
//     if((a < 0) && (T(int(a)) == a)) { // If a is negative integer
//         ArrayXx data = ArrayXx::LinSpaced(n, -a + 1, -a - n + 2);
//         data(0) = 1;
//         data = data.log();
//         set_cumsum(data, ans);
//     } else {
//         ans = ArrayXx::LinSpaced(n, a, a + n - 1).lgamma();
//         ans -= lgamma(a);
//     }
//     return ans;
// }
// template ArrayXd get_lrf(const double a, const int n);
// template ArrayXl get_lrf(const long double a, const int n);

// Overloaded functions get_lrf for double and long double
// Non-template versions are required because Eigen does not have
// lgamma for a long double Array.
Eigen::ArrayXd get_lrf(const double a, const int n) {
    ArrayXd ans(n);
    if((a < 0) && (double(int(a)) == a)) { // If a is negative integer
        ArrayXd data = ArrayXd::LinSpaced(n, -a + 1, -a - n + 2);
        data(0) = 1;
        data = data.log();
        set_cumsum(data, ans);
    } else {
        ans = ArrayXd::LinSpaced(n, a, a + n - 1).lgamma();
        ans -= lgamma(a);
    }
    return ans;
}

Eigen::Array<long double, Eigen::Dynamic, 1> get_lrf(const long double a, const int n) {
    ArrayXl ans(n);
    if((a < 0) && ((long double)(int(a)) == a)) { // If a is negative integer
        ArrayXl data = ArrayXl::LinSpaced(n, -a + 1, -a - n + 2);
        data(0) = 1;
        data = data.log();
        set_cumsum(data, ans);
    } else {
        for(int i = 0; i < n; i++) ans[i] = std::lgammal(a + i) - std::lgammal(a);
    }
    return ans;
}


// Eigen function template to obtain the signs of a series of rising factorial (a)_k
// from parameter a and length n (k = 0, 1, ... n - 1)
// And template instantiations for double and long double
template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1>  get_sign_rf(const T a, const int n) {
    typedef Eigen::Array<T, Dynamic, 1> ArrayXx;
    ArrayXx ans(n);
    ArrayXx Signs = sign(ArrayXx::LinSpaced(n, a - 1, a + n - 2));
    Signs(0) = 1;
    std::partial_sum(Signs.data(), Signs.data() + n, ans.data(), std::multiplies<T>());
    return ans;
}
template ArrayXd get_sign_rf(const double a, const int n);
template ArrayXl get_sign_rf(const long double a, const int n);

// Eigen function template to obtain the signs of a series of rising factorial (a+1)_k
// And template instantiations for double and long double
template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1> get_sign_rfp1(const T a, const int n) {
    typedef Eigen::Array<T, Dynamic, 1> ArrayXx;
    ArrayXx ans(n);
    ArrayXx Signs = sign(ArrayXx::LinSpaced(n, a, a + n - 1));
    std::partial_sum(Signs.data(), Signs.data() + n, ans.data(), std::multiplies<T>());
    return ans;
}
template ArrayXd get_sign_rfp1(const double a, const int n);
// template ArrayXl get_sign_rfp1(const long double a, const int n);

// // This function makes the following n * (n * n) Matrix
// // from a Vector {0, 1, ... n}:
// // 0 0 .. 0 1 1 .. 1 ... n n .. n
// // 0 0 .. 0 1 1 .. 1 ... n n .. n
// // . . .. . . . .. . ... . . .. .
// // 0 . .. 0 1 1 .. 1 ... n n .. n
// template <typename Derived>
// Eigen::Matrix<typename Derived::Scalar, Dynamic, Dynamic> ResizeFor3d(const Eigen::MatrixBase<Derived>& X) {
//     const int n = X.rows();
//     Matrix<typename Derived::Scalar, Dynamic, Dynamic> mX = X.transpose().colwise().replicate(n);
//     mX.resize(n * n, 1);
//     return mX.transpose().colwise().replicate(n);
// }

// Eigen version of hgs_1d()
//
// Takes lscf as a separate parameter
//
// // [[Rcpp::export]]
Eigen::ArrayXd hgs_1dE(const Eigen::ArrayXd& dks,
                       const double a1, const double b, const double lconst,
                       const Eigen::ArrayXd& lscf) {
    const int m = dks.size() - 1;
    ArrayXd Alnum = get_lrf(a1, m + 1);
    ArrayXd Alden = get_lrf(b, m + 1);
    ArrayXd ansseq(m + 1);
    ArrayXd Asgns = get_sign_rf(a1, m + 1);
    ansseq = exp(Alnum - Alden + log(abs(dks)) + lconst - lscf);
    ansseq *= Asgns * sign(dks);
    return ansseq;
}


// // Eigen template version of hgs_2d()
// //
// // // [[Rcpp::export]]
// template <typename Derived>
// Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
// hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
//         const typename Derived::Scalar a1, const typename Derived::Scalar a2,
//         const double b, const typename Derived::Scalar lconst,
//         const Eigen::ArrayBase<Derived>& lscf) {
//     const int m = dks.rows() - 1;
//     typedef Matrix<typename Derived::Scalar, Dynamic, 1> VectorXx;
//     typedef Array<typename Derived::Scalar, Dynamic, 1> ArrayXx;
//     typedef Array<typename Derived::Scalar, Dynamic, Dynamic> ArrayXXx;
//     // For denominator Array, calculated as double and then casted
//     VectorXd seq0m = VectorXd::LinSpaced(m + 1, 0, m);
//     ArrayXx Alnumi = get_lrf(a1, m + 1);
//     ArrayXx Alnumj = get_lrf(a2, m + 1);
//     ArrayXXd Aldend(m + 1, m + 1);
//     Aldend = lgamma(b + (seq0m.rowwise().replicate(m + 1) +
//                            seq0m.transpose().colwise().replicate(m + 1)).array());
//     Aldend -= lgamma(b);
//     ArrayXXx Alden = Aldend.cast<typename Derived::Scalar>();
//     ArrayXXx ansmat = ArrayXXx::Zero(m + 1, m + 1);
//     ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
//     ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
//     ansmat.colwise() += Alnumi;
//     ansmat.rowwise() += Alnumj.transpose();
//     ansmat -= Alden;
//     ansmat += log(abs(dks)) + lconst;
//     ansmat -= lscf;
//     ansmat = exp(ansmat);
//     ansmat.colwise() *= Asgnsi;
//     ansmat.rowwise() *= Asgnsj.transpose();
//     ansmat *= sign(dks);
//     return ansmat;
// }
// // template instantiations for double and long double cases
// template ArrayXXd hgs_2dE(const ArrayBase<ArrayXXd>& dks,
//                           const double a1, const double a2,
//                           const double b, const double lconst,
//                           const ArrayBase<ArrayXXd>& lscf);
// template ArrayXXl hgs_2dE(const ArrayBase<ArrayXXl>& dks,
//                           const long double a1, const long double a2,
//                           const double b, const long double lconst,
//                           const ArrayBase<ArrayXXl>& lscf);

// Eigen version of hgs_2d(); double
//
// // [[Rcpp::export]]
Eigen::ArrayXXd
hgs_2dE(const Eigen::ArrayXXd& dks,
        const double a1, const double a2, const double b,
        // const double lconst, const Eigen::ArrayXXd& lscf) {
        const double lconst, const Eigen::ArrayXd& lscf) {
    const int m = dks.rows() - 1;
    ArrayXd Alnumi = get_lrf(a1, m + 1);
    ArrayXd Alnumj = get_lrf(a2, m + 1);
    // VectorXd seq0m = VectorXd::LinSpaced(m + 1, 0, m);
    // ArrayXXd Alden(m + 1, m + 1);
    // Alden = lgamma(b + (seq0m.rowwise().replicate(m + 1) +
    //                        seq0m.transpose().colwise().replicate(m + 1)).array());
    // Alden -= lgamma(b);
    ArrayXXd ansmat = ArrayXXd::Zero(m + 1, m + 1);
    ArrayXd Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXd Asgnsj = get_sign_rf(a2, m + 1);
    ansmat.colwise() += Alnumi;
    ansmat.rowwise() += Alnumj.transpose();
    // ansmat -= Alden;
    for(int k = 0; k <= m; k++) {
        double lden = std::lgamma(b + k) - std::lgamma(b);
        for(int i = 0; i <= k; i++) ansmat(i, k - i) -= lden;
    }
    ansmat += log(abs(dks)) + lconst;
    // ansmat -= lscf;
    for(int k = 0; k <= m; k++) {
        for(int i = 0; i <= k; i++) ansmat(i, k - i) -= lscf(k);
    }
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    ansmat.rowwise() *= Asgnsj.transpose();
    ansmat *= sign(dks);
    return ansmat;
}

// Eigen version of hgs_2d(); long double
//
// // [[Rcpp::export]]
Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>
hgs_2dE(const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& dks,
        const long double a1, const long double a2, const long double b,
        const long double lconst,
        // const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& lscf) {
        const Eigen::Array<long double, Eigen::Dynamic, 1>& lscf) {
    const int m = dks.rows() - 1;
    ArrayXl Alnumi = get_lrf(a1, m + 1);
    ArrayXl Alnumj = get_lrf(a2, m + 1);
    ArrayXXl ansmat = ArrayXXl::Zero(m + 1, m + 1);
    ArrayXl Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXl Asgnsj = get_sign_rf(a2, m + 1);
    ansmat.colwise() += Alnumi;
    ansmat.rowwise() += Alnumj.transpose();
    for(int k = 0; k <= m; k++) {
        long double lden = std::lgammal(b + k) - std::lgammal(b);
        for(int i = 0; i <= k; i++) ansmat(i, k - i) -= lden;
    }
    ansmat += log(abs(dks)) + lconst;
    // ansmat -= lscf;
    for(int k = 0; k <= m; k++) {
        for(int i = 0; i <= k; i++) ansmat(i, k - i) -= lscf(k);
    }
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    ansmat.rowwise() *= Asgnsj.transpose();
    ansmat *= sign(dks);
    return ansmat;
}


// // Eigen temlate version of hgs_3d()
// //
// // // [[Rcpp::export]]
// template <typename Derived>
// Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
// hgs_3dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1,
//         const typename Derived::Scalar a2, const typename Derived::Scalar a3,
//         const double b, const typename Derived::Scalar lconst,
//         const Eigen::ArrayBase<Derived>& lscf) {
//     const int m = dks.rows() - 1;
//     typedef Matrix<typename Derived::Scalar, Dynamic, 1> VectorXx;
//     typedef Array<typename Derived::Scalar, Dynamic, 1> ArrayXx;
//     typedef Array<typename Derived::Scalar, Dynamic, Dynamic> ArrayXXx;
//     VectorXd seq0m = VectorXd::LinSpaced(m + 1, 0, m);
//     ArrayXx Alnumi = get_lrf(a1, m + 1);
//     ArrayXx Alnumj = get_lrf(a2, m + 1);
//     ArrayXx Alnumk = get_lrf(a3, m + 1);
//     ArrayXXd Aldend(m + 1, (m + 1) * (m + 1));
//     Aldend = (b + (seq0m.rowwise().replicate(m + 1) +
//                       seq0m.transpose().colwise().replicate(m + 1)).rowwise().replicate(m + 1).array() +
//              ResizeFor3d(seq0m).array()).lgamma();
//     Aldend -= lgamma(b);
//     ArrayXXx Alden = Aldend.cast<typename Derived::Scalar>();
//     ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
//     ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
//     ArrayXx Asgnsk = get_sign_rf(a3, m + 1);
//     ArrayXXx ansmat = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
//     ansmat.colwise() += Alnumi;
//     for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() += Alnumj.transpose();
//     for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) += Alnumk(k);
//     ansmat -= Alden;
//     ansmat += log(abs(dks)) + lconst;
//     ansmat -= lscf;
//     ansmat = exp(ansmat);
//     ansmat.colwise() *= Asgnsi;
//     for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() *= Asgnsj.transpose();
//     for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) *= Asgnsk(k);
//     ansmat *= sign(dks);
//     return ansmat;
// }
// template ArrayXXd hgs_3dE(const ArrayBase<ArrayXXd>& dks, const double a1,
//                           const double a2, const double a3,
//                           const double b, const double lconst,
//                           const ArrayBase<ArrayXXd>& lscf);
// template ArrayXXl hgs_3dE(const ArrayBase<ArrayXXl>& dks, const long double a1,
//                           const long double a2, const long double a3,
//                           const double b, const long double lconst,
//                           const ArrayBase<ArrayXXl>& lscf);

// Eigen version of hgs_3d(), double
//
// // [[Rcpp::export]]
Eigen::ArrayXXd
hgs_3dE(const Eigen::ArrayXXd& dks,
        const double a1, const double a2, const double a3,
        // const double b, const double lconst, const Eigen::ArrayXXd& lscf) {
        const double b, const double lconst, const Eigen::ArrayXd& lscf) {
    const int m = dks.rows() - 1;
    ArrayXd Alnumi = get_lrf(a1, m + 1);
    ArrayXd Alnumj = get_lrf(a2, m + 1);
    ArrayXd Alnumk = get_lrf(a3, m + 1);
    // VectorXd seq0m = VectorXd::LinSpaced(m + 1, 0, m);
    // ArrayXXd Alden(m + 1, (m + 1) * (m + 1));
    // Alden = (b + (seq0m.rowwise().replicate(m + 1) +
    //                   seq0m.transpose().colwise().replicate(m + 1)).rowwise().replicate(m + 1).array() +
    //          ResizeFor3d(seq0m).array()).lgamma();
    // Alden -= lgamma(b);
    ArrayXd Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXd Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXd Asgnsk = get_sign_rf(a3, m + 1);
    ArrayXXd ansmat = ArrayXXd::Zero(m + 1, (m + 1) * (m + 1));
    ansmat.colwise() += Alnumi;
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() += Alnumj.transpose();
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) += Alnumk(k);
    // ansmat -= Alden;
    for(int k = 0; k <= m; k++) {
        double lden = std::lgamma(b + k) - std::lgamma(b);
        for(int i = 0; i <= k; i++) {
            for(int j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lden;
            }
        }
    }
    ansmat += log(abs(dks)) + lconst;
    // ansmat -= lscf;
    for(int k = 0; k <= m; k++) {
        for(int i = 0; i <= k; i++) {
            for(int j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lscf(k);
            }
        }
    }
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() *= Asgnsj.transpose();
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) *= Asgnsk(k);
    ansmat *= sign(dks);
    return ansmat;
}

// Eigen version of hgs_3d(), long double
//
// // [[Rcpp::export]]
Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>
hgs_3dE(const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& dks,
        const long double a1, const long double a2, const long double a3,
        const long double b, const long double lconst,
        // const Eigen::Array<long double, Eigen::Dynamic, Eigen::Dynamic>& lscf) {
        const Eigen::Array<long double, Eigen::Dynamic, 1>& lscf) {
    const int m = dks.rows() - 1;
    ArrayXl Alnumi = get_lrf(a1, m + 1);
    ArrayXl Alnumj = get_lrf(a2, m + 1);
    ArrayXl Alnumk = get_lrf(a3, m + 1);
    ArrayXl Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXl Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXl Asgnsk = get_sign_rf(a3, m + 1);
    ArrayXXl ansmat = ArrayXXl::Zero(m + 1, (m + 1) * (m + 1));
    ansmat.colwise() += Alnumi;
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() += Alnumj.transpose();
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) += Alnumk(k);
    for(int k = 0; k <= m; k++) {
        double lden = std::lgammal(b + k) - std::lgammal(b);
        for(int i = 0; i <= k; i++) {
            for(int j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lden;
            }
        }
    }
    ansmat += log(abs(dks)) + lconst;
    // ansmat -= lscf;
    for(int k = 0; k <= m; k++) {
        for(int i = 0; i <= k; i++) {
            for(int j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lscf(k);
            }
        }
    }
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() *= Asgnsj.transpose();
    for(int k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) *= Asgnsk(k);
    ansmat *= sign(dks);
    return ansmat;
}
