#include "config.h"
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

template <typename Derived>
bool is_zero_E(const Eigen::ArrayBase<Derived>& X, const typename Derived::Scalar tol) {
    return (X.abs() <= tol).all();
}
template bool is_zero_E<ArrayXd>(const Eigen::ArrayBase<ArrayXd>& X, const double tol);
template bool is_zero_E<ArrayXl>(const Eigen::ArrayBase<ArrayXl>& X, const long double tol);

template <typename Derived>
bool is_diag_E(const Eigen::MatrixBase<Derived>& X, const typename Derived::Scalar tol) {
    Matrix<typename Derived::Scalar, Dynamic, Dynamic> Xa = X;
    Xa.diagonal().setZero();
    return is_zero_E(Xa.array(), tol);
}
template bool is_diag_E(const Eigen::MatrixBase<MatrixXd>& X, const double tol);
template bool is_diag_E(const Eigen::MatrixBase<MatrixXl>& X, const long double tol);


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
// Specialized for ULT object (see config.h)
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_2dE(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf) {
    const Index m = dks.ULT_getM() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXx ansmat = ArrayXx::Zero((m + 1) * (m + 2) / 2);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULTcol(k, m + 1) += Alnumi.head(m + 1 - k) + Alnumj(k);
        Scalar lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) -= lden;
    }
    ansmat += log(abs(dks)) + lconst;
    for(Index k = 0; k <= m; k++) {
        for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) -= lscf(k);
    }
    ansmat = exp(ansmat);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULTcol(k, m + 1) *= Asgnsi.head(m + 1 - k) * Asgnsj(k);
    }
    ansmat *= sign(dks);
    return ansmat;
}
// template instantiations for double and long double cases
template ArrayXd hgs_2dE(const ArrayBase<ArrayXd>& dks,
                         const double a1, const double a2,
                         const double b, const double lconst,
                         const Array<double, Dynamic, 1>& lscf);
template ArrayXl hgs_2dE(const ArrayBase<ArrayXl>& dks,
                         const long double a1, const long double a2,
                         const long double b, const long double lconst,
                         const Array<long double, Dynamic, 1>& lscf);

// Eigen template version of hgs_2d(); coefficient-wise lscf
// Specialized for ULT object (see config.h)
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_2dEc(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf) {
    const Index m = dks.ULT_getM() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXx ansmat = ArrayXx::Zero((m + 1) * (m + 2) / 2);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULTcol(k, m + 1) += Alnumi.head(m + 1 - k) + Alnumj(k);
        Scalar lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index i = 0; i <= k; i++) ansmat.ULTat(i, k - i, m + 1) -= lden;
    }
    ansmat += log(abs(dks)) + lconst;
    ansmat -= lscf;
    ansmat = exp(ansmat);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULTcol(k, m + 1) *= Asgnsi.head(m + 1 - k) * Asgnsj(k);
    }
    ansmat *= sign(dks);
    return ansmat;
}
template ArrayXd hgs_2dEc(const ArrayBase<ArrayXd> &dks,
                          const double a1, const double a2,
                          const double b, const double lconst,
                          const ArrayBase<ArrayXd> &lscf);


// Eigen temlate version of hgs_3d()
// Specialized for ULC object (see config.h)
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_3dE(const Eigen::ArrayBase<Derived>& dks, const typename Derived::Scalar a1,
        const typename Derived::Scalar a2, const typename Derived::Scalar a3,
        const typename Derived::Scalar b, const typename Derived::Scalar lconst,
        const Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>& lscf) {
    const Index m = dks.ULC_getM() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXx Alnumk = get_lrf(a3, m + 1);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXx Asgnsk = get_sign_rf(a3, m + 1);
    ArrayXx ansmat = ArrayXx::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULCslice(k, m + 1) += Alnumk(k);
        Scalar lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index j = 0; j <= k; j++) {
            ansmat.ULCcol(j, k - j, m + 1) += Alnumi.head(m + 1 - k) + Alnumj(j);
            for(Index i = 0; i <= k - j; i++) {
                ansmat.ULCat(i, j, k - i - j, m + 1) -= lden;
            }
        }
    }
    ansmat += log(abs(dks)) + lconst;
    for(Index k = 0; k <= m; k++) {
        for(Index i = 0; i <= k; i++) {
            for(Index j = 0; j <= k - i; j++) {
                ansmat.ULCat(i, j, k - i - j, m + 1) -= lscf(k);
            }
        }
    }
    ansmat = exp(ansmat);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULCslice(k, m + 1) *= Asgnsk(k);
        for(Index j = 0; j <= k; j++) {
            ansmat.ULCcol(j, k - j, m + 1) *= Asgnsi.head(m + 1 - k) * Asgnsj(j);
        }
    }
    ansmat *= sign(dks);
    return ansmat;
}
template ArrayXd hgs_3dE(const ArrayBase<ArrayXd>& dks, const double a1,
                         const double a2, const double a3,
                         const double b, const double lconst,
                         const Array<double, Dynamic, 1>& lscf);
template ArrayXl hgs_3dE(const ArrayBase<ArrayXl>& dks, const long double a1,
                         const long double a2, const long double a3,
                         const long double b, const long double lconst,
                         const Array<long double, Dynamic, 1>& lscf);

// Eigen template version of hgs_3d(), coefficient-wise lscf
// Specialized for ULC object (see config.h)
//
// // [[Rcpp::export]]
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
hgs_3dEc(const Eigen::ArrayBase<Derived>& dks,
        const typename Derived::Scalar a1, const typename Derived::Scalar a2,
        const typename Derived::Scalar a3, const typename Derived::Scalar b,
        const typename Derived::Scalar lconst,
        const Eigen::ArrayBase<Derived>& lscf) {
    const Index m = dks.ULC_getM() - 1;
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    ArrayXx Alnumi = get_lrf(a1, m + 1);
    ArrayXx Alnumj = get_lrf(a2, m + 1);
    ArrayXx Alnumk = get_lrf(a3, m + 1);
    ArrayXx Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXx Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXx Asgnsk = get_sign_rf(a3, m + 1);
    ArrayXx ansmat = ArrayXx::Zero((m + 1) * (m + 2) * (m + 3) / 6);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULCslice(k, m + 1) += Alnumk(k);
        Scalar lden = std::lgamma(b + k) - std::lgamma(b);
        for(Index j = 0; j <= k; j++) {
            ansmat.ULCcol(j, k - j, m + 1) += Alnumi.head(m + 1 - k) + Alnumj(j);
            for(Index i = 0; i <= k - j; i++) {
                ansmat.ULCat(i, j, k - i - j, m + 1) -= lden;
            }
        }
    }
    ansmat += log(abs(dks)) + lconst;
    ansmat -= lscf;
    ansmat = exp(ansmat);
    for(Index k = 0; k <= m; k++) {
        ansmat.ULCslice(k, m + 1) *= Asgnsk(k);
        for(Index j = 0; j <= k; j++) {
            ansmat.ULCcol(j, k - j, m + 1) *= Asgnsi.head(m + 1 - k) * Asgnsj(j);
        }
    }
    ansmat *= sign(dks);
    return ansmat;
}
template ArrayXd hgs_3dEc(const ArrayBase<ArrayXd>& dks, const double a1, 
                          const double a2, const double a3, const double b,
                          const double lconst, const ArrayBase<ArrayXd>& lscf);


// Eigen template version of \code{sum_counterdiag()}
// Specialized for ULT object (see config.h)
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiagE(const Eigen::ArrayBase<Derived>& X) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = (sqrt(8 * X.size() + 2) - 1) / 2;
    ArrayXx ans = ArrayXx::Zero(n);
    Scalar x;
    for(Index i = 0; i < n; i++) {
        for(Index j = 0; j <= i; j++) {
            x = X.ULTat(i - j, j, n);
            if(!std::isnan(x)) ans(i) += x;
        }
    }
    return ans;
}
template ArrayXd sum_counterdiagE(const ArrayBase<ArrayXd>& X);
template ArrayXl sum_counterdiagE(const ArrayBase<ArrayXl>& X);

// Eigen template version of \code{sum_counterdiag3D()}
// Specialized for ULC object (see config.h)
template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
sum_counterdiag3DE(const Eigen::ArrayBase<Derived>& X) {
    typedef typename Derived::Scalar Scalar;
    typedef Array<Scalar, Dynamic, 1> ArrayXx;
    const Index n = X.ULC_getM();
    ArrayXx ans = ArrayXx::Zero(n);
    Scalar x;
    for(Index i = 0; i < n; i++) {
        for(Index j = 0; j <= i; j++) {
            for(Index k = 0; k <= (i - j); k++){
                x = X.ULCat(i - j - k, j, k, n);
                if(!std::isnan(x)) ans(i) += x;
            }
        }
    }
    return ans;
}
template ArrayXd sum_counterdiag3DE(const ArrayBase<ArrayXd>& X);
template ArrayXl sum_counterdiag3DE(const ArrayBase<ArrayXl>& X);
