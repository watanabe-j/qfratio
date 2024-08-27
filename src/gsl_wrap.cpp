#include <Rcpp.h>

// These are to use gsl
#include <RcppGSL.h>
// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>

using Rcpp::IntegerVector;
using Rcpp::IntegerMatrix;
using Rcpp::NumericVector;
using Rcpp::NumericMatrix;
using std::size_t;


//' @describeIn gsl_wrap
//'   wrapper of \code{gsl_hyperg_1F1_e()}, looping along \code{bvec}
//'
// [[Rcpp::export]]
SEXP hyperg_1F1_vec_b(const double a, const Rcpp::NumericVector bvec,
                      const double x) {
    gsl_sf_result hgtmp;
    gsl_set_error_handler_off();
    size_t n = bvec.size();
    NumericVector hgres(n);
    NumericVector hgerr(n);
    IntegerVector hgstatus(n);
    for(size_t i = 0; i < n; i++) {
        hgstatus(i) = gsl_sf_hyperg_1F1_e(a, bvec(i), x, &hgtmp);
        hgres(i) = hgtmp.val;
        hgerr(i) = hgtmp.err;
    }
    return Rcpp::List::create(
        Rcpp::Named("val")    = hgres,
        Rcpp::Named("err")    = hgerr,
        Rcpp::Named("status") = hgstatus);
}

//' @describeIn gsl_wrap
//'   wrapper of \code{gsl_hyperg_2F1_e()}, looping along \code{Amat} and
//'   recycling \code{cvec}
//'
// [[Rcpp::export]]
SEXP hyperg_2F1_mat_a_vec_c(const Rcpp::NumericMatrix Amat, const double b,
                            const Rcpp::NumericVector cvec, const double x) {
    gsl_sf_result hgtmp;
    gsl_set_error_handler_off();
    size_t n = cvec.size();
    NumericMatrix hgres(n, n);
    NumericMatrix hgerr(n, n);
    IntegerMatrix hgstatus(n, n);
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < n - i; j++) {
            hgstatus(j, i) = gsl_sf_hyperg_2F1_e(Amat(j, i), b, cvec(j), x,
                                                 &hgtmp);
            hgres(j, i) = hgtmp.val;
            hgerr(j, i) = hgtmp.err;
        }
    }
    return Rcpp::List::create(
        Rcpp::Named("val")    = hgres,
        Rcpp::Named("err")    = hgerr,
        Rcpp::Named("status") = hgstatus);
}
