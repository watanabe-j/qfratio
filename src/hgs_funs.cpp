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
using Eigen::Index;

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
// Eigen::Array<T, Dynamic, 1> get_lrf(const T a, const Eigen::Index n) {
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
// template ArrayXd get_lrf(const double a, const Index n);
// template ArrayXl get_lrf(const long double a, const Index n);

// Overloaded functions get_lrf for double and long double
// Non-template versions are required because Eigen does not have
// lgamma for a long double Array.
Eigen::ArrayXd get_lrf(const double a, const Eigen::Index n) {
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

Eigen::Array<long double, Eigen::Dynamic, 1> get_lrf(const long double a, const Eigen::Index n) {
    ArrayXl ans(n);
    if((a < 0) && ((long double)(int(a)) == a)) { // If a is negative integer
        ArrayXl data = ArrayXl::LinSpaced(n, -a + 1, -a - n + 2);
        data(0) = 1;
        data = data.log();
        set_cumsum(data, ans);
    } else {
        for(Index i = 0; i < n; i++) ans[i] = std::lgammal(a + i) - std::lgammal(a);
    }
    return ans;
}


// Eigen function template to obtain the signs of a series of rising factorial (a)_k
// from parameter a and length n (k = 0, 1, ... n - 1)
// And template instantiations for double and long double
template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1>  get_sign_rf(const T a, const Eigen::Index n) {
    typedef Eigen::Array<T, Dynamic, 1> ArrayXx;
    ArrayXx ans(n);
    ArrayXx Signs = sign(ArrayXx::LinSpaced(n, a - 1, a + n - 2));
    Signs(0) = 1;
    std::partial_sum(Signs.data(), Signs.data() + n, ans.data(), std::multiplies<T>());
    return ans;
}
template ArrayXd get_sign_rf(const double a, const Index n);
template ArrayXl get_sign_rf(const long double a, const Index n);

// Eigen function template to obtain the signs of a series of rising factorial (a+1)_k
// And template instantiations for double and long double
template <typename T>
Eigen::Array<T, Eigen::Dynamic, 1> get_sign_rfp1(const T a, const Eigen::Index n) {
    typedef Eigen::Array<T, Dynamic, 1> ArrayXx;
    ArrayXx ans(n);
    ArrayXx Signs = sign(ArrayXx::LinSpaced(n, a, a + n - 1));
    std::partial_sum(Signs.data(), Signs.data() + n, ans.data(), std::multiplies<T>());
    return ans;
}
template ArrayXd get_sign_rfp1(const double a, const Index n);
// template ArrayXl get_sign_rfp1(const long double a, const Index n);


// Eigen version of hgs_1d()
//
// Takes lscf as a separate parameter
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_1dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1, 
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf) {
    const Index m = dks.size() - 1;
    typedef Array<typename Derived::Scalar, Dynamic, 1> ArrayXx;
    ArrayXx Alnum = get_lrf(a1, m + 1);
    ArrayXx Alden = get_lrf(b, m + 1);
    ArrayXx ansseq(m + 1);
    ArrayXx Asgns = get_sign_rf(a1, m + 1);
    ansseq = exp(Alnum - Alden + log(abs(dks)) + lconst - lscf);
    ansseq *= Asgns * sign(dks);
    return ansseq;
}
template ArrayXd hgs_1dE(const ArrayBase<ArrayXd>& dks, const double a1, 
                         const double b, const double lconst,
                         const ArrayBase<ArrayXd>& lscf);


// Eigen template version of hgs_2d()
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf) {
    const Index m = dks.rows() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXXx ansmat = ArrayXXx::Zero(m + 1, m + 1);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    ansmat.colwise() += Alnumi;
    ansmat.rowwise() += Alnumj.transpose();
    for(Index k = 0; k <= m; k++) {
        double lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index i = 0; i <= k; i++) ansmat(i, k - i) -= lden;
    }
    ansmat += log(abs(dks)) + lconst;
    for(Index k = 0; k <= m; k++) {
        for(Index i = 0; i <= k; i++) ansmat(i, k - i) -= lscf(k);
    }
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    ansmat.rowwise() *= Asgnsj.transpose();
    ansmat *= sign(dks);
    return ansmat;
}
// template instantiations for double and long double cases
template ArrayXXd hgs_2dE(const ArrayBase<ArrayXXd>& dks,
                          const double a1, const double a2,
                          const double b, const double lconst,
                          const Array<double, Dynamic, 1>& lscf);
template ArrayXXl hgs_2dE(const ArrayBase<ArrayXXl>& dks,
                          const long double a1, const long double a2,
                          const long double b, const long double lconst,
                          const Array<long double, Dynamic, 1>& lscf);

// Eigen template version of hgs_2d(); coefficient-wise lscf
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf) {
    const Index m = dks.rows() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXXx ansmat = ArrayXXx::Zero(m + 1, m + 1);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    ansmat.colwise() += Alnumi;
    ansmat.rowwise() += Alnumj.transpose();
    for(Index k = 0; k <= m; k++) {
        double lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index i = 0; i <= k; i++) ansmat(i, k - i) -= lden;
    }
    ansmat += log(abs(dks)) + lconst;
    ansmat -= lscf;
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    ansmat.rowwise() *= Asgnsj.transpose();
    ansmat *= sign(dks);
    return ansmat;
}
template ArrayXXd hgs_2dE(const ArrayBase<ArrayXXd> &dks,
                          const double a1, const double a2,
                          const double b, const double lconst,
                          const ArrayBase<ArrayXXd> &lscf);


// Eigen temlate version of hgs_3d()
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_3dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1,
        const typename Derived::Scalar a2, const typename Derived::Scalar a3,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf) {
    const Index m = dks.rows() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXx Alnumk = get_lrf(a3, m + 1);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXx Asgnsk = get_sign_rf(a3, m + 1);
    ArrayXXx ansmat = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
    ansmat.colwise() += Alnumi;
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() += Alnumj.transpose();
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) += Alnumk(k);
    for(Index k = 0; k <= m; k++) {
        double lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index i = 0; i <= k; i++) {
            for(Index j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lden;
            }
        }
    }
    ansmat += log(abs(dks)) + lconst;
    for(Index k = 0; k <= m; k++) {
        for(Index i = 0; i <= k; i++) {
            for(Index j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lscf(k);
            }
        }
    }
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() *= Asgnsj.transpose();
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) *= Asgnsk(k);
    ansmat *= sign(dks);
    return ansmat;
}
template ArrayXXd hgs_3dE(const ArrayBase<ArrayXXd>& dks, const double a1,
                          const double a2, const double a3,
                          const double b, const double lconst,
                          const Array<double, Dynamic, 1>& lscf);
template ArrayXXl hgs_3dE(const ArrayBase<ArrayXXl>& dks, const long double a1,
                          const long double a2, const long double a3,
                          const long double b, const long double lconst,
                          const Array<long double, Dynamic, 1>& lscf);

// Eigen template version of hgs_3d(), coefficient-wise lscf
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>
hgs_3dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar a3, const typename Derived::Scalar b,
        const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf) {
    const Index m = dks.rows() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    typedef Array<Scalar, Dynamic, Dynamic> ArrayXXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXx Alnumk = get_lrf(a3, m + 1);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXx Asgnsk = get_sign_rf(a3, m + 1);
    ArrayXXx ansmat = ArrayXXx::Zero(m + 1, (m + 1) * (m + 1));
    ansmat.colwise() += Alnumi;
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() += Alnumj.transpose();
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) += Alnumk(k);
    for(Index k = 0; k <= m; k++) {
        double lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index i = 0; i <= k; i++) {
            for(Index j = 0; j <= k - i; j++) {
                ansmat(i, j + (k - i - j) * (m + 1)) -= lden;
            }
        }
    }
    ansmat += log(abs(dks)) + lconst;
    ansmat -= lscf;
    ansmat = exp(ansmat);
    ansmat.colwise() *= Asgnsi;
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1).rowwise() *= Asgnsj.transpose();
    for(Index k = 0; k <= m; k++) ansmat.block(0, k * (m + 1), m + 1, m + 1) *= Asgnsk(k);
    ansmat *= sign(dks);
    return ansmat;
}
template ArrayXXd hgs_3dE(const ArrayBase<ArrayXXd>& dks, const double a1, 
                          const double a2, const double a3, const double b,
                          const double lconst, const ArrayBase<ArrayXXd>& lscf);


// Eigen template version of \code{sum_counterdiag()}
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiagE(const Eigen::ArrayBase<Derived>& X) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = X.rows();
    ArrayXx ans = ArrayXx::Zero(n);
    Scalar x;
    for(Index i = 0; i < n; i++) {
        for(Index j = 0; j <= i; j++) {
            x = X(i - j, j);
            if(!std::isnan(x)) ans(i) += x;
        }
    }
    return ans;
}
template ArrayXd sum_counterdiagE(const ArrayBase<ArrayXXd>& X);
template ArrayXl sum_counterdiagE(const ArrayBase<ArrayXXl>& X);

// Eigen template version of \code{sum_counterdiag3D()}
// X is a wide ArrayXXl, n * (n * n)
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiag3DE(const Eigen::ArrayBase<Derived>& X) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = X.rows();
    ArrayXx ans = ArrayXx::Zero(n);
    Scalar x;
    for(Index i = 0; i < n; i++) {
        for(Index j = 0; j <= i; j++) {
            for(Index k = 0; k <= (i - j); k++){
                x = X(i - j - k, j + k * n);
                if(!std::isnan(x)) ans(i) += x;
            }
        }
    }
    return ans;
}
template ArrayXd sum_counterdiag3DE(const ArrayBase<ArrayXXd>& X);
template ArrayXl sum_counterdiag3DE(const ArrayBase<ArrayXXl>& X);
