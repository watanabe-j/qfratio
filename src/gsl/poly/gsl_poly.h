// This file is taken from GSL version 2.7.1 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.
// - 2023 Junya Watanabe

/* poly/gsl_poly.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Brian Gough
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

#ifndef __GSL_POLY_H__
#define __GSL_POLY_H__

#include <stdlib.h>
#include "../gsl_inline.h" // edited for qfratio
// #include <gsl/gsl_complex.h> // edited for qfratio

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


/* Evaluate polynomial
 *
 * c[0] + c[1] x + c[2] x^2 + ... + c[len-1] x^(len-1)
 *
 * exceptions: none
 */

/* real polynomial, real x */
INLINE_DECL double gsl_poly_eval(const double c[], const int len, const double x);

// /* real polynomial, complex x */ // edited for qfratio
// INLINE_DECL gsl_complex gsl_poly_complex_eval (const double c [], const int len, const gsl_complex z); // edited for qfratio

// /* complex polynomial, complex x */ // edited for qfratio
// INLINE_DECL gsl_complex gsl_complex_poly_complex_eval (const gsl_complex c [], const int len, const gsl_complex z); // edited for qfratio

// int gsl_poly_eval_derivs(const double c[], const size_t lenc, const double x, double res[], const size_t lenres); // edited for qfratio

#ifdef HAVE_INLINE
INLINE_FUN
double 
gsl_poly_eval(const double c[], const int len, const double x)
{
  int i;
  double ans = c[len-1];
  for(i=len-1; i>0; i--) ans = c[i-1] + x * ans;
  return ans;
}

// INLINE_FUN // edited for qfratio
// gsl_complex // edited for qfratio
// gsl_poly_complex_eval(const double c[], const int len, const gsl_complex z) // edited for qfratio
// { // edited for qfratio
//   int i; // edited for qfratio
//   gsl_complex ans; // edited for qfratio
//   GSL_SET_COMPLEX (&ans, c[len-1], 0.0); // edited for qfratio
//   for(i=len-1; i>0; i--) { // edited for qfratio
//     /* The following three lines are equivalent to // edited for qfratio
//        ans = gsl_complex_add_real (gsl_complex_mul (z, ans), c[i-1]);  // edited for qfratio
//        but faster */ // edited for qfratio
//     double tmp = c[i-1] + GSL_REAL (z) * GSL_REAL (ans) - GSL_IMAG (z) * GSL_IMAG (ans); // edited for qfratio
//     GSL_SET_IMAG (&ans, GSL_IMAG (z) * GSL_REAL (ans) + GSL_REAL (z) * GSL_IMAG (ans)); // edited for qfratio
//     GSL_SET_REAL (&ans, tmp); // edited for qfratio
//   }  // edited for qfratio
//   return ans; // edited for qfratio
// } // edited for qfratio

// INLINE_FUN // edited for qfratio
// gsl_complex // edited for qfratio
// gsl_complex_poly_complex_eval(const gsl_complex c[], const int len, const gsl_complex z) // edited for qfratio
// { // edited for qfratio
//   int i; // edited for qfratio
//   gsl_complex ans = c[len-1]; // edited for qfratio
//   for(i=len-1; i>0; i--) { // edited for qfratio
//     /* The following three lines are equivalent to // edited for qfratio
//        ans = gsl_complex_add (c[i-1], gsl_complex_mul (x, ans)); // edited for qfratio
//        but faster */ // edited for qfratio
//     double tmp = GSL_REAL (c[i-1]) + GSL_REAL (z) * GSL_REAL (ans) - GSL_IMAG (z) * GSL_IMAG (ans); // edited for qfratio
//     GSL_SET_IMAG (&ans, GSL_IMAG (c[i-1]) + GSL_IMAG (z) * GSL_REAL (ans) + GSL_REAL (z) * GSL_IMAG (ans)); // edited for qfratio
//     GSL_SET_REAL (&ans, tmp); // edited for qfratio
//   } // edited for qfratio
//   return ans; // edited for qfratio
// } // edited for qfratio
#endif /* HAVE_INLINE */

// /* Work with divided-difference polynomials, Abramowitz & Stegun 25.2.26 */ // edited for qfratio

// int // edited for qfratio
// gsl_poly_dd_init (double dd[], const double x[], const double y[], // edited for qfratio
//                   size_t size); // edited for qfratio

// INLINE_DECL double // edited for qfratio
// gsl_poly_dd_eval (const double dd[], const double xa[], const size_t size, const double x); // edited for qfratio

// #ifdef HAVE_INLINE // edited for qfratio
// INLINE_FUN // edited for qfratio
// double  // edited for qfratio
// gsl_poly_dd_eval(const double dd[], const double xa[], const size_t size, const double x) // edited for qfratio
// { // edited for qfratio
//   size_t i; // edited for qfratio
//   double y = dd[size - 1]; // edited for qfratio
//   for (i = size - 1; i--;) y = dd[i] + (x - xa[i]) * y; // edited for qfratio
//   return y; // edited for qfratio
// } // edited for qfratio
// #endif /* HAVE_INLINE */ // edited for qfratio


// int // edited for qfratio
// gsl_poly_dd_taylor (double c[], double xp, // edited for qfratio
//                     const double dd[], const double x[], size_t size, // edited for qfratio
//                     double w[]); // edited for qfratio

// int // edited for qfratio
// gsl_poly_dd_hermite_init (double dd[], double z[], const double xa[], const double ya[], // edited for qfratio
//                           const double dya[], const size_t size); // edited for qfratio

// /* Solve for real or complex roots of the standard quadratic equation, // edited for qfratio
//  * returning the number of real roots. // edited for qfratio
//  * // edited for qfratio
//  * Roots are returned ordered. // edited for qfratio
//  */ // edited for qfratio
// int gsl_poly_solve_quadratic (double a, double b, double c,  // edited for qfratio
//                               double * x0, double * x1); // edited for qfratio

// int  // edited for qfratio
// gsl_poly_complex_solve_quadratic (double a, double b, double c,  // edited for qfratio
//                                   gsl_complex * z0, gsl_complex * z1); // edited for qfratio


// /* Solve for real roots of the cubic equation // edited for qfratio
//  * x^3 + a x^2 + b x + c = 0, returning the // edited for qfratio
//  * number of real roots. // edited for qfratio
//  * // edited for qfratio
//  * Roots are returned ordered. // edited for qfratio
//  */ // edited for qfratio
// int gsl_poly_solve_cubic (double a, double b, double c,  // edited for qfratio
//                           double * x0, double * x1, double * x2); // edited for qfratio

// int  // edited for qfratio
// gsl_poly_complex_solve_cubic (double a, double b, double c,  // edited for qfratio
//                               gsl_complex * z0, gsl_complex * z1,  // edited for qfratio
//                               gsl_complex * z2); // edited for qfratio


// /* Solve for the complex roots of a general real polynomial */ // edited for qfratio

// typedef struct  // edited for qfratio
// {  // edited for qfratio
//   size_t nc ; // edited for qfratio
//   double * matrix ;  // edited for qfratio
// }  // edited for qfratio
// gsl_poly_complex_workspace ; // edited for qfratio

// gsl_poly_complex_workspace * gsl_poly_complex_workspace_alloc (size_t n); // edited for qfratio
// void gsl_poly_complex_workspace_free (gsl_poly_complex_workspace * w); // edited for qfratio

// int // edited for qfratio
// gsl_poly_complex_solve (const double * a, size_t n,  // edited for qfratio
//                         gsl_poly_complex_workspace * w, // edited for qfratio
//                         gsl_complex_packed_ptr z); // edited for qfratio

__END_DECLS

#endif /* __GSL_POLY_H__ */
