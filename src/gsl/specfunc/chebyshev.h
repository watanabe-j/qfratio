// This file is taken from GSL version 2.7.1 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.

/* specfunc/chebyshev.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* data for a Chebyshev series over a given interval */

#ifndef __GSL__CHEBYSHEV_H__ // added for qfratio
#define __GSL__CHEBYSHEV_H__ // added for qfratio

struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;


#endif // added for qfratio


