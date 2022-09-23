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


void set_cumsum(const Eigen::ArrayXd& Data, Eigen::ArrayXd& Out) {
    std::partial_sum(Data.data(), Data.data() + Data.size(), Out.data());
}

// void set_cumprod_sign(const double a, const int n, Eigen::ArrayXd& Out) {
//     ArrayXd Data = sign(ArrayXd::LinSpaced(n, a, a + n - 1));
//     std::partial_sum(Data.data(), Data.data() + n, Out.data(), std::multiplies<double>());
// }

// void set_cumprod_sign_rf(const double a, const int n, Eigen::ArrayXd& Out) {
//     ArrayXd Data = sign(ArrayXd::LinSpaced(n, a - 1, a + n - 2));
//     Data(0) = 1;
//     std::partial_sum(Data.data(), Data.data() + n, Out.data(), std::multiplies<double>());
// }

// template<typename T> inline T sgn(T x) {
//     return (x > T(0)) - (x < T(0));
// }

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

Eigen::ArrayXd get_sign_rf(const double a, const int n) {
    ArrayXd ans(n);
    ArrayXd Signs = sign(ArrayXd::LinSpaced(n, a - 1, a + n - 2));
    Signs(0) = 1;
    std::partial_sum(Signs.data(), Signs.data() + n, ans.data(), std::multiplies<double>());
    return ans;
}

Eigen::ArrayXd get_sign_rfp1(const double a, const int n) {
    ArrayXd ans(n);
    ArrayXd Signs = sign(ArrayXd::LinSpaced(n, a, a + n - 1));
    std::partial_sum(Signs.data(), Signs.data() + n, ans.data(), std::multiplies<double>());
    return ans;
}

// This function makes the following n * (n * n) Matrix
// from a Vector {0, 1, ... n}:
// 0 0 .. 0 1 1 .. 1 ... n n .. n
// 0 0 .. 0 1 1 .. 1 ... n n .. n
// . . .. . . . .. . ... . . .. .
// 0 . .. 0 1 1 .. 1 ... n n .. n
Eigen::MatrixXd ResizeFor3d(const Eigen::VectorXd& X) {
    const int n = X.rows();
    MatrixXd mX = X.transpose().colwise().replicate(n);
    mX.resize(n * n, 1);
    return mX.transpose().colwise().replicate(n);
}

// C++ version of hgs_1d
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
    // ArrayXd Asgns(m + 1);
    // set_cumprod_sign_rf(a1, m + 1, Asgns);
    // Asgns(0) = 1;
    // for(int i = 1; i <= m; i ++) {
    //     Asgns(i) = Asgns(i - 1) * sgn(a1 + double(i) - 1);
    // }
    ansseq = exp(Alnum - Alden + log(abs(dks)) + lconst - lscf);
    ansseq *= Asgns * sign(dks);
    return ansseq;
}


// C++ version of hgs_2d
//
// // [[Rcpp::export]]
Eigen::ArrayXXd hgs_2dE(const Eigen::ArrayXXd& dks,
                        const double a1, const double a2,
                        const double b, const double lconst,
                        const Eigen::ArrayXXd& lscf) {
    const int m = dks.rows() - 1;
    VectorXd seq0m = VectorXd::LinSpaced(m + 1, 0, m);
    VectorXd Alnumi = get_lrf(a1, m + 1);
    VectorXd Alnumj = get_lrf(a2, m + 1);
    ArrayXXd Alden(m + 1, m + 1);
    Alden = lgamma(b + (seq0m.rowwise().replicate(m + 1) +
                           seq0m.transpose().colwise().replicate(m + 1)).array());
    Alden -= lgamma(b);
    ArrayXXd ansmat(m + 1, m + 1);
    ArrayXd Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXd Asgnsj = get_sign_rf(a2, m + 1);
    // ArrayXd Asgnsi(m + 1);
    // ArrayXd Asgnsj(m + 1);
    // set_cumprod_sign_rf(a1, m + 1, Asgnsi);
    // set_cumprod_sign_rf(a2, m + 1, Asgnsj);
    // Asgnsi(0) = 1;
    // Asgnsj(0) = 1;
    // for(int i = 1; i <= m; i++) {
    //     Asgnsi(i) = Asgnsi(i - 1) * sgn(a1 + double(i) - 1);
    //     Asgnsj(i) = Asgnsj(i - 1) * sgn(a2 + double(i) - 1);
    // }
    ansmat = exp((Alnumi.rowwise().replicate(m + 1) +
                  Alnumj.transpose().colwise().replicate(m + 1)).array()
                 - Alden + log(abs(dks)) + lconst - lscf);
    ansmat *= (Asgnsi.matrix() * Asgnsj.matrix().transpose()).array() * sign(dks);
    return ansmat;
}


// C++ version of hgs_3d
//
// // [[Rcpp::export]]
Eigen::ArrayXXd hgs_3dE(const Eigen::ArrayXXd& dks,
                        const double a1, const double a2, const double a3,
                        const double b, const double lconst,
                        const Eigen::ArrayXXd& lscf) {
    const int m = dks.rows() - 1;
    VectorXd seq0m = VectorXd::LinSpaced(m + 1, 0, m);
    VectorXd Alnumi = get_lrf(a1, m + 1);
    VectorXd Alnumj = get_lrf(a2, m + 1);
    VectorXd Alnumk = get_lrf(a3, m + 1);
    ArrayXXd Alden(m + 1, (m + 1) * (m + 1));
    Alden = (b + (seq0m.rowwise().replicate(m + 1) +
                      seq0m.transpose().colwise().replicate(m + 1)).rowwise().replicate(m + 1).array() +
             ResizeFor3d(seq0m).array()).lgamma();
    Alden -= lgamma(b);
    ArrayXXd ansmat(m + 1, (m + 1) * (m + 1));
    ArrayXd Asgnsi = get_sign_rf(a1, m + 1);
    ArrayXd Asgnsj = get_sign_rf(a2, m + 1);
    ArrayXd Asgnsk = get_sign_rf(a3, m + 1);
    // ArrayXd Asgnsi(m + 1);
    // ArrayXd Asgnsj(m + 1);
    // ArrayXd Asgnsk(m + 1);
    // set_cumprod_sign_rf(a1, m + 1, Asgnsi);
    // set_cumprod_sign_rf(a2, m + 1, Asgnsj);
    // set_cumprod_sign_rf(a3, m + 1, Asgnsk);
    // Asgnsi(0) = 1;
    // Asgnsj(0) = 1;
    // Asgnsk(0) = 1;
    // for(int i = 1; i <= m; i++) {
    //     Asgnsi(i) = Asgnsi(i - 1) * sgn(a1 + double(i) - 1);
    //     Asgnsj(i) = Asgnsj(i - 1) * sgn(a2 + double(i) - 1);
    //     Asgnsk(i) = Asgnsk(i - 1) * sgn(a3 + double(i) - 1);
    // }
    ansmat = exp( (Alnumi.rowwise().replicate(m + 1) +
                   Alnumj.transpose().colwise().replicate(m + 1) ).rowwise().replicate(m + 1).array() +
                   ResizeFor3d(Alnumk).array()
                 - Alden + log(abs(dks)) + lconst - lscf);
    ansmat *= ( (Asgnsi.matrix() * Asgnsj.matrix().transpose()).rowwise().replicate(m + 1) ).array() *
              ResizeFor3d(Asgnsk.matrix()).array() * sign(dks);
    return ansmat;
}
