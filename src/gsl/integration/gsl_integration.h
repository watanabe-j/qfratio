// This file is taken from GSL version 2.8 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.
// - 2023 Junya Watanabe

/* integration/gsl_integration.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_INTEGRATION_H__
#define __GSL_INTEGRATION_H__
#include <stdlib.h>
#include "../gsl_math.h" // edited for qfratio

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Workspace for adaptive integrators */

typedef struct
  {
    size_t limit;
    size_t size;
    size_t nrmax;
    size_t i;
    size_t maximum_level;
    double *alist;
    double *blist;
    double *rlist;
    double *elist;
    size_t *order;
    size_t *level;
  }
gsl_integration_workspace;

gsl_integration_workspace *
  gsl_integration_workspace_alloc (const size_t n);

void
  gsl_integration_workspace_free (gsl_integration_workspace * w);


/* Workspace for QAWS integrator */

typedef struct
{
  double alpha;
  double beta;
  int mu;
  int nu;
  double ri[25];
  double rj[25];
  double rg[25];
  double rh[25];
}
gsl_integration_qaws_table;

gsl_integration_qaws_table * 
gsl_integration_qaws_table_alloc (double alpha, double beta, int mu, int nu);

int
gsl_integration_qaws_table_set (gsl_integration_qaws_table * t,
                                double alpha, double beta, int mu, int nu);

void
gsl_integration_qaws_table_free (gsl_integration_qaws_table * t);

// /* Workspace for QAWO integrator */ // edited for qfratio

// enum gsl_integration_qawo_enum { GSL_INTEG_COSINE, GSL_INTEG_SINE }; // edited for qfratio

// typedef struct // edited for qfratio
// { // edited for qfratio
//   size_t n; // edited for qfratio
//   double omega; // edited for qfratio
//   double L; // edited for qfratio
//   double par; // edited for qfratio
//   enum gsl_integration_qawo_enum sine; // edited for qfratio
//   double *chebmo; // edited for qfratio
// } // edited for qfratio
// gsl_integration_qawo_table; // edited for qfratio

// gsl_integration_qawo_table *  // edited for qfratio
// gsl_integration_qawo_table_alloc (double omega, double L,  // edited for qfratio
//                                   enum gsl_integration_qawo_enum sine, // edited for qfratio
//                                   size_t n); // edited for qfratio

// int // edited for qfratio
// gsl_integration_qawo_table_set (gsl_integration_qawo_table * t, // edited for qfratio
//                                 double omega, double L, // edited for qfratio
//                                 enum gsl_integration_qawo_enum sine); // edited for qfratio

// int // edited for qfratio
// gsl_integration_qawo_table_set_length (gsl_integration_qawo_table * t, // edited for qfratio
//                                        double L); // edited for qfratio

// void // edited for qfratio
// gsl_integration_qawo_table_free (gsl_integration_qawo_table * t); // edited for qfratio


/* Definition of an integration rule */

typedef void gsl_integration_rule (const gsl_function * f,
                                   double a, double b,
                                   double *result, double *abserr,
                                   double *defabs, double *resabs);

void gsl_integration_qk15 (const gsl_function * f, double a, double b,
                           double *result, double *abserr,
                           double *resabs, double *resasc);

// void gsl_integration_qk21 (const gsl_function * f, double a, double b, // edited for qfratio
//                            double *result, double *abserr, // edited for qfratio
//                            double *resabs, double *resasc); // edited for qfratio

// void gsl_integration_qk31 (const gsl_function * f, double a, double b, // edited for qfratio
//                            double *result, double *abserr, // edited for qfratio
//                            double *resabs, double *resasc); // edited for qfratio

// void gsl_integration_qk41 (const gsl_function * f, double a, double b, // edited for qfratio
//                            double *result, double *abserr, // edited for qfratio
//                            double *resabs, double *resasc); // edited for qfratio

// void gsl_integration_qk51 (const gsl_function * f, double a, double b, // edited for qfratio
//                            double *result, double *abserr, // edited for qfratio
//                            double *resabs, double *resasc); // edited for qfratio

// void gsl_integration_qk61 (const gsl_function * f, double a, double b, // edited for qfratio
//                            double *result, double *abserr, // edited for qfratio
//                            double *resabs, double *resasc); // edited for qfratio

// void gsl_integration_qcheb (gsl_function * f, double a, double b,  // edited for qfratio
//                             double *cheb12, double *cheb24); // edited for qfratio

/* The low-level integration rules in QUADPACK are identified by small
   integers (1-6). We'll use symbolic constants to refer to them.  */

enum
  {
    GSL_INTEG_GAUSS15 = 1,      /* 15 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS21 = 2,      /* 21 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS31 = 3,      /* 31 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS41 = 4,      /* 41 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS51 = 5,      /* 51 point Gauss-Kronrod rule */
    GSL_INTEG_GAUSS61 = 6       /* 61 point Gauss-Kronrod rule */
  };

void 
gsl_integration_qk (const int n, const double xgk[], 
                    const double wg[], const double wgk[],
                    double fv1[], double fv2[],
                    const gsl_function *f, double a, double b,
                    double * result, double * abserr, 
                    double * resabs, double * resasc);


// int gsl_integration_qng (const gsl_function * f, // edited for qfratio
//                          double a, double b, // edited for qfratio
//                          double epsabs, double epsrel, // edited for qfratio
//                          double *result, double *abserr, // edited for qfratio
//                          size_t * neval); // edited for qfratio

// int gsl_integration_qag (const gsl_function * f, // edited for qfratio
//                          double a, double b, // edited for qfratio
//                          double epsabs, double epsrel, size_t limit, // edited for qfratio
//                          int key, // edited for qfratio
//                          gsl_integration_workspace * workspace, // edited for qfratio
//                          double *result, double *abserr); // edited for qfratio

// int gsl_integration_qagi (gsl_function * f, // edited for qfratio
//                           double epsabs, double epsrel, size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           double *result, double *abserr); // edited for qfratio

int gsl_integration_qagiu (gsl_function * f,
                           double a,
                           double epsabs, double epsrel, size_t limit,
                           gsl_integration_workspace * workspace,
                           double *result, double *abserr);

// int gsl_integration_qagil (gsl_function * f, // edited for qfratio
//                            double b, // edited for qfratio
//                            double epsabs, double epsrel, size_t limit, // edited for qfratio
//                            gsl_integration_workspace * workspace, // edited for qfratio
//                            double *result, double *abserr); // edited for qfratio


// int gsl_integration_qags (const gsl_function * f, // edited for qfratio
//                           double a, double b, // edited for qfratio
//                           double epsabs, double epsrel, size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           double *result, double *abserr); // edited for qfratio

// int gsl_integration_qagp (const gsl_function * f, // edited for qfratio
//                           double *pts, size_t npts, // edited for qfratio
//                           double epsabs, double epsrel, size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           double *result, double *abserr); // edited for qfratio

// int gsl_integration_qawc (gsl_function *f, // edited for qfratio
//                           const double a, const double b, const double c, // edited for qfratio
//                           const double epsabs, const double epsrel, const size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           double * result, double * abserr); // edited for qfratio

// int gsl_integration_qaws (gsl_function * f, // edited for qfratio
//                           const double a, const double b, // edited for qfratio
//                           gsl_integration_qaws_table * t, // edited for qfratio
//                           const double epsabs, const double epsrel, // edited for qfratio
//                           const size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           double *result, double *abserr); // edited for qfratio

// int gsl_integration_qawo (gsl_function * f, // edited for qfratio
//                           const double a, // edited for qfratio
//                           const double epsabs, const double epsrel, // edited for qfratio
//                           const size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           gsl_integration_qawo_table * wf, // edited for qfratio
//                           double *result, double *abserr); // edited for qfratio

// int gsl_integration_qawf (gsl_function * f, // edited for qfratio
//                           const double a, // edited for qfratio
//                           const double epsabs, // edited for qfratio
//                           const size_t limit, // edited for qfratio
//                           gsl_integration_workspace * workspace, // edited for qfratio
//                           gsl_integration_workspace * cycle_workspace, // edited for qfratio
//                           gsl_integration_qawo_table * wf, // edited for qfratio
//                           double *result, double *abserr); // edited for qfratio

// /* Workspace for fixed-order Gauss-Legendre integration */ // edited for qfratio

// typedef struct // edited for qfratio
//   { // edited for qfratio
//     size_t n;         /* number of points */ // edited for qfratio
//     double *x;        /* Gauss abscissae/points */ // edited for qfratio
//     double *w;        /* Gauss weights for each abscissae */ // edited for qfratio
//     int precomputed;  /* high precision abscissae/weights precomputed? */ // edited for qfratio
//   } // edited for qfratio
// gsl_integration_glfixed_table; // edited for qfratio


// gsl_integration_glfixed_table * gsl_integration_glfixed_table_alloc (size_t n); // edited for qfratio

// void gsl_integration_glfixed_table_free (gsl_integration_glfixed_table * t); // edited for qfratio

// /* Routine for fixed-order Gauss-Legendre integration */ // edited for qfratio

// double gsl_integration_glfixed (const gsl_function *f, // edited for qfratio
//                                 double a, // edited for qfratio
//                                 double b, // edited for qfratio
//                                 const gsl_integration_glfixed_table * t); // edited for qfratio

// /* Routine to retrieve the i-th Gauss-Legendre point and weight from t */ // edited for qfratio

// int gsl_integration_glfixed_point (double a, // edited for qfratio
//                                    double b, // edited for qfratio
//                                    size_t i, // edited for qfratio
//                                    double *xi, // edited for qfratio
//                                    double *wi, // edited for qfratio
//                                    const gsl_integration_glfixed_table * t); // edited for qfratio


// /* Cquad integration - Pedro Gonnet */ // edited for qfratio

// /* Data of a single interval */ // edited for qfratio
// typedef struct // edited for qfratio
// { // edited for qfratio
//   double a, b; // edited for qfratio
//   double c[64]; // edited for qfratio
//   double fx[33]; // edited for qfratio
//   double igral, err; // edited for qfratio
//   int depth, rdepth, ndiv; // edited for qfratio
// } gsl_integration_cquad_ival; // edited for qfratio


// /* The workspace is just a collection of intervals */ // edited for qfratio
// typedef struct // edited for qfratio
// { // edited for qfratio
//   size_t size; // edited for qfratio
//   gsl_integration_cquad_ival *ivals; // edited for qfratio
//   size_t *heap; // edited for qfratio
// } gsl_integration_cquad_workspace; // edited for qfratio

// gsl_integration_cquad_workspace * // edited for qfratio
// gsl_integration_cquad_workspace_alloc (const size_t n); // edited for qfratio

// void // edited for qfratio
// gsl_integration_cquad_workspace_free (gsl_integration_cquad_workspace * w); // edited for qfratio

// int // edited for qfratio
// gsl_integration_cquad (const gsl_function * f, double a, double b, // edited for qfratio
// 		                   double epsabs, double epsrel, // edited for qfratio
// 		                   gsl_integration_cquad_workspace * ws, // edited for qfratio
// 		                   double *result, double *abserr, size_t * nevals); // edited for qfratio

// /* Romberg integration workspace and routines */ // edited for qfratio

// typedef struct // edited for qfratio
// { // edited for qfratio
//   size_t n;       /* maximum number of steps */ // edited for qfratio
//   double *work1;  /* workspace for a row of R matrix, size n */ // edited for qfratio
//   double *work2;  /* workspace for a row of R matrix, size n */ // edited for qfratio
// } gsl_integration_romberg_workspace; // edited for qfratio

// gsl_integration_romberg_workspace *gsl_integration_romberg_alloc(const size_t n); // edited for qfratio
// void gsl_integration_romberg_free(gsl_integration_romberg_workspace * w); // edited for qfratio
// int gsl_integration_romberg(const gsl_function * f, const double a, const double b, // edited for qfratio
//                             const double epsabs, const double epsrel, double * result, // edited for qfratio
//                             size_t * neval, gsl_integration_romberg_workspace * w); // edited for qfratio

// /* IQPACK related structures and routines */ // edited for qfratio

// typedef struct // edited for qfratio
// { // edited for qfratio
//   double alpha; // edited for qfratio
//   double beta; // edited for qfratio
//   double a; // edited for qfratio
//   double b; // edited for qfratio
//   double zemu; // edited for qfratio
//   double shft; // edited for qfratio
//   double slp; // edited for qfratio
//   double al; // edited for qfratio
//   double be; // edited for qfratio
// } gsl_integration_fixed_params; // edited for qfratio

// typedef struct // edited for qfratio
// { // edited for qfratio
//   int (*check)(const size_t n, const gsl_integration_fixed_params * params); // edited for qfratio
//   int (*init)(const size_t n, double * diag, double * subdiag, gsl_integration_fixed_params * params); // edited for qfratio
// } gsl_integration_fixed_type; // edited for qfratio

// typedef struct // edited for qfratio
// { // edited for qfratio
//   size_t n;        /* number of nodes/weights */ // edited for qfratio
//   double *weights; /* quadrature weights */ // edited for qfratio
//   double *x;       /* quadrature nodes */ // edited for qfratio
//   double *diag;    /* diagonal of Jacobi matrix */ // edited for qfratio
//   double *subdiag; /* subdiagonal of Jacobi matrix */ // edited for qfratio
//   const gsl_integration_fixed_type * type; // edited for qfratio
// } gsl_integration_fixed_workspace; // edited for qfratio

// /* IQPACK integral types */ // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_legendre; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_chebyshev; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_gegenbauer; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_jacobi; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_laguerre; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_hermite; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_exponential; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_rational; // edited for qfratio
// GSL_VAR const gsl_integration_fixed_type * gsl_integration_fixed_chebyshev2; // edited for qfratio

// gsl_integration_fixed_workspace * // edited for qfratio
// gsl_integration_fixed_alloc(const gsl_integration_fixed_type * type, const size_t n, // edited for qfratio
//                             const double a, const double b, const double alpha, const double beta); // edited for qfratio

// void gsl_integration_fixed_free(gsl_integration_fixed_workspace * w); // edited for qfratio

// size_t gsl_integration_fixed_n(const gsl_integration_fixed_workspace * w); // edited for qfratio

// double *gsl_integration_fixed_nodes(const gsl_integration_fixed_workspace * w); // edited for qfratio

// double *gsl_integration_fixed_weights(const gsl_integration_fixed_workspace * w); // edited for qfratio

// int gsl_integration_fixed(const gsl_function * func, double * result, // edited for qfratio
//                           const gsl_integration_fixed_workspace * w); // edited for qfratio

// /* Lebedev quadrature */ // edited for qfratio

// typedef struct // edited for qfratio
// { // edited for qfratio
//   size_t n;        /* number of nodes/weights */ // edited for qfratio
//   double *weights; /* quadrature weights */ // edited for qfratio
//   double *x;       /* x quadrature nodes */ // edited for qfratio
//   double *y;       /* y quadrature nodes */ // edited for qfratio
//   double *z;       /* z quadrature nodes */ // edited for qfratio
//   double *theta;   /* theta quadrature nodes */ // edited for qfratio
//   double *phi;     /* phi quadrature nodes */ // edited for qfratio
// } gsl_integration_lebedev_workspace; // edited for qfratio

// gsl_integration_lebedev_workspace * gsl_integration_lebedev_alloc(const size_t n); // edited for qfratio

// void gsl_integration_lebedev_free(gsl_integration_lebedev_workspace * w); // edited for qfratio

// size_t gsl_integration_lebedev_n(const gsl_integration_lebedev_workspace * w); // edited for qfratio

__END_DECLS

#endif /* __GSL_INTEGRATION_H__ */
