// This file is taken from GSL version 2.7.1 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.
// - 2023 Junya Watanabe

/* roots/convergence.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough
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

#include <config.h>
#include "../gsl_math.h" // edited for qfratio
#include "../err/gsl_errno.h" // edited for qfratio
#include "gsl_roots.h" // edited for qfratio

// int // edited for qfratio
// gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel) // edited for qfratio
// { // edited for qfratio
//   const double abs_lower = fabs(x_lower) ; // edited for qfratio
//   const double abs_upper = fabs(x_upper) ; // edited for qfratio

//   double min_abs, tolerance; // edited for qfratio

//   if (epsrel < 0.0) // edited for qfratio
//     GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL); // edited for qfratio
//    // edited for qfratio
//   if (epsabs < 0.0) // edited for qfratio
//     GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL); // edited for qfratio

//   if (x_lower > x_upper) // edited for qfratio
//     GSL_ERROR ("lower bound larger than upper bound", GSL_EINVAL); // edited for qfratio

//   if ((x_lower > 0.0 && x_upper > 0.0) || (x_lower < 0.0 && x_upper < 0.0))  // edited for qfratio
//     { // edited for qfratio
//       min_abs = GSL_MIN_DBL(abs_lower, abs_upper) ; // edited for qfratio
//     } // edited for qfratio
//   else // edited for qfratio
//     { // edited for qfratio
//       min_abs = 0; // edited for qfratio
//     } // edited for qfratio

//   tolerance = epsabs + epsrel * min_abs  ; // edited for qfratio
//    // edited for qfratio
//   if (fabs(x_upper - x_lower) < tolerance) // edited for qfratio
//     return GSL_SUCCESS; // edited for qfratio
//    // edited for qfratio
//   return GSL_CONTINUE ; // edited for qfratio
// } // edited for qfratio

int
gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel)
{
  const double tolerance = epsabs + epsrel * fabs(x1)  ;

  if (epsrel < 0.0)
    GSL_ERROR ("relative tolerance is negative", GSL_EBADTOL);
  
  if (epsabs < 0.0)
    GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
  
  if (fabs(x1 - x0) < tolerance || x1 == x0)
    return GSL_SUCCESS;
  
  return GSL_CONTINUE ;
}

// int // edited for qfratio
// gsl_root_test_residual (double f, double epsabs) // edited for qfratio
// { // edited for qfratio
//   if (epsabs < 0.0) // edited for qfratio
//     GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL); // edited for qfratio
//   // edited for qfratio
//   if (fabs(f) < epsabs) // edited for qfratio
//     return GSL_SUCCESS; // edited for qfratio
//    // edited for qfratio
//   return GSL_CONTINUE ; // edited for qfratio
// } // edited for qfratio

