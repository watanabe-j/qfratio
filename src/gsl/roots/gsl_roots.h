// This file is taken from GSL version 2.7.1 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.

/* roots/gsl_roots.h
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

#ifndef __GSL_ROOTS_H__
#define __GSL_ROOTS_H__

#include <stdlib.h>
#include "../gsl_types.h" // edited for qfratio
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

typedef struct
  {
    const char *name;
    size_t size;
    int (*set) (void *state, gsl_function * f, double * root, double x_lower, double x_upper);
    int (*iterate) (void *state, gsl_function * f, double * root, double * x_lower, double * x_upper);
  }
gsl_root_fsolver_type;

typedef struct
  {
    const gsl_root_fsolver_type * type;
    gsl_function * function ;
    double root ;
    double x_lower;
    double x_upper;
    void *state;
  }
gsl_root_fsolver;

// typedef struct // edited for qfratio
//   { // edited for qfratio
//     const char *name; // edited for qfratio
//     size_t size; // edited for qfratio
//     int (*set) (void *state, gsl_function_fdf * f, double * root); // edited for qfratio
//     int (*iterate) (void *state, gsl_function_fdf * f, double * root); // edited for qfratio
//   } // edited for qfratio
// gsl_root_fdfsolver_type; // edited for qfratio

// typedef struct // edited for qfratio
//   { // edited for qfratio
//     const gsl_root_fdfsolver_type * type; // edited for qfratio
//     gsl_function_fdf * fdf ; // edited for qfratio
//     double root ; // edited for qfratio
//     void *state; // edited for qfratio
//   } // edited for qfratio
// gsl_root_fdfsolver; // edited for qfratio

gsl_root_fsolver *
gsl_root_fsolver_alloc (const gsl_root_fsolver_type * T);
void gsl_root_fsolver_free (gsl_root_fsolver * s);

int gsl_root_fsolver_set (gsl_root_fsolver * s,
                          gsl_function * f, 
                          double x_lower, double x_upper);

int gsl_root_fsolver_iterate (gsl_root_fsolver * s);

const char * gsl_root_fsolver_name (const gsl_root_fsolver * s);
double gsl_root_fsolver_root (const gsl_root_fsolver * s);
double gsl_root_fsolver_x_lower (const gsl_root_fsolver * s);
double gsl_root_fsolver_x_upper (const gsl_root_fsolver * s);


// gsl_root_fdfsolver * // edited for qfratio
// gsl_root_fdfsolver_alloc (const gsl_root_fdfsolver_type * T); // edited for qfratio

// int // edited for qfratio
// gsl_root_fdfsolver_set (gsl_root_fdfsolver * s,  // edited for qfratio
//                          gsl_function_fdf * fdf, double root); // edited for qfratio

// int // edited for qfratio
// gsl_root_fdfsolver_iterate (gsl_root_fdfsolver * s); // edited for qfratio

// void // edited for qfratio
// gsl_root_fdfsolver_free (gsl_root_fdfsolver * s); // edited for qfratio

// const char * gsl_root_fdfsolver_name (const gsl_root_fdfsolver * s); // edited for qfratio
// double gsl_root_fdfsolver_root (const gsl_root_fdfsolver * s); // edited for qfratio

int
gsl_root_test_interval (double x_lower, double x_upper, double epsabs, double epsrel);

int
gsl_root_test_residual (double f, double epsabs);

int
gsl_root_test_delta (double x1, double x0, double epsabs, double epsrel);

GSL_VAR const gsl_root_fsolver_type  * gsl_root_fsolver_bisection;
GSL_VAR const gsl_root_fsolver_type  * gsl_root_fsolver_brent;
GSL_VAR const gsl_root_fsolver_type  * gsl_root_fsolver_falsepos;
// GSL_VAR const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_newton; // edited for qfratio
// GSL_VAR const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_secant; // edited for qfratio
// GSL_VAR const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson; // edited for qfratio

__END_DECLS

#endif /* __GSL_ROOTS_H__ */
