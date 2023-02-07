// This header is to add user-defined methods in Eigen::ArrayBase via config.h
// These concern indexing of coefficients and "columns" in a packed
// Upper-Left Triangular (ULT) object, which stores the upper-left elements of
// a column-major M x M array in an ArrayXx:
//    0 ..j . M cols
//  0 x x x x
//  . x x x o
//  i x x o o
//  . x o o o
//  M rows    M * (M + 1) / 2 coefs
// Second, Upper-Left Stacked (ULS) objects assumes ULT stacked row-wise
// in an ArrayXXx (used in, e.g., d3_pjk_mEc)


// #include <RcppEigen.h>
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen; // automatically assumed

// Get M for ULT from size() of ArrayXx
inline Index ULT_getM() const {
    return (std::sqrt(8 * this->size() + 1) - 1) / 2;
}
inline Index ULT_getM() {
    return (std::sqrt(8 * this->size() + 1) - 1) / 2;
}

// (i,j)-th coef in ULT
inline Scalar ULTat(Index i, Index j) const {
    Index M = ULT_getM();
    return this->operator()(i + j * (2 * M - j + 1) / 2);
}
inline Scalar& ULTat(Index i, Index j) {
    Index M = ULT_getM();
    return this->operator()(i + j * (2 * M - j + 1) / 2);
}
// Same with known M for performance
inline Scalar ULTat(Index i, Index j, Index M) const {
    return this->operator()(i + j * (2 * M - j + 1) / 2);
}
inline Scalar& ULTat(Index i, Index j, Index M) {
    return this->operator()(i + j * (2 * M - j + 1) / 2);
}

// j-th "column" in ULT
inline VectorBlock<Derived> ULTcol(Index j) const {
    Index M = ULT_getM();
    return this->segment(j * (2 * M - j + 1) / 2, M - j);
}
inline VectorBlock<Derived> ULTcol(Index j) {
    Index M = ULT_getM();
    return this->segment(j * (2 * M - j + 1) / 2, M - j);
}
inline VectorBlock<Derived> ULTcol(Index j, Index M) const {
    return this->segment(j * (2 * M - j + 1) / 2, M - j);
}
inline VectorBlock<Derived> ULTcol(Index j, Index M) {
    return this->segment(j * (2 * M - j + 1) / 2, M - j);
}


// Get M for ULS from cols() of ArrayXXx
inline Index ULS_getM() const {
    return (std::sqrt(8 * this->cols() + 1) - 1) / 2;
}
inline Index ULS_getM() {
    return (std::sqrt(8 * this->cols() + 1) - 1) / 2;
}

// (i,j,k)-th coef in ULS
inline Scalar ULSat(Index i, Index j, Index k) const {
    Index M = ULS_getM();
    return this->operator()(i, j + k * (2 * M - k + 1) / 2);
}
inline Scalar& ULSat(Index i, Index j, Index k) {
    Index M = ULS_getM();
    return this->operator()(i, j + k * (2 * M - k + 1) / 2);
}
inline Scalar ULSat(Index i, Index j, Index k, Index M) const {
    return this->operator()(i, j + k * (2 * M - k + 1) / 2);
}
inline Scalar& ULSat(Index i, Index j, Index k, Index M) {
    return this->operator()(i, j + k * (2 * M - k + 1) / 2);
}

// (j,k)-th "column" in ULS
inline Block<Derived, Dynamic, 1, true> ULScol(Index j, Index k) const {
    Index M = ULS_getM();
    return this->col(j + k * (2 * M - k + 1) / 2);
}
inline Block<Derived, Dynamic, 1, true> ULScol(Index j, Index k) {
    Index M = ULS_getM();
    return this->col(j + k * (2 * M - k + 1) / 2);
}
inline Block<Derived, Dynamic, 1, true> ULScol(Index j, Index k, Index M) const {
    return this->col(j + k * (2 * M - k + 1) / 2);
}
inline Block<Derived, Dynamic, 1, true> ULScol(Index j, Index k, Index M) {
    return this->col(j + k * (2 * M - k + 1) / 2);
}

// // segment in j-th "column" in ULT
// inline VectorBlock<Derived> ULTcolseg(Index j, Index start, Index size) const {
//     Index M = ULT_getM();
//     return this->segment(start + j * (2 * M - j + 1) / 2, size);
// }
// inline VectorBlock<Derived> ULTcolseg(Index j, Index start, Index size) {
//     Index M = ULT_getM();
//     return this->segment(start + j * (2 * M - j + 1) / 2, size);
// }
// inline VectorBlock<Derived> ULTcolseg(Index j, Index start, Index size, Index M) const {
//     return this->segment(start + j * (2 * M - j + 1) / 2, size);
// }
// inline VectorBlock<Derived> ULTcolseg(Index j, Index start, Index size, Index M) {
//     return this->segment(start + j * (2 * M - j + 1) / 2, size);
// }

