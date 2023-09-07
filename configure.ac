## This file is part of qfratio
##
## This is for a configure script to generate src/config.h
## which defines necessary macros for Eigen and GSL.
##
## Run autoreconf to (re-)generate configure and src/config.h.in
## (or autoconf to generate configure and autoheader to generate config.h.in)
##
## GSL-related parts in this file were taken from configure.ac in
## GSL version 2.7.1.
## The original line numbers are specified at the top of each chunk.
##
## Copyright (C) 2023 Junya Watanabe
##           (C) ND   Free Software Foundation, Inc. (for GSL-originated part)
## Licensed under the GNU GPL (>=3). See LISENSE.md

AC_INIT([qfratio], [1.0.1.9000], [https://github.com/watanabe-j/qfratio/issues])

AC_CONFIG_SRCDIR([R/qfratio-package.R])

## config.guess and config.sub necessary for AC_CANONICAL_HOST below
## These are in tools/ as recommended in  'Writing R Extensions'
## Run src/scripts/get_config_aux.sh to download the latest versions
AC_CONFIG_AUX_DIR(tools)


#####
## Use the same C compiler used in R, from 'Writing R Extensions'
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`


dnl #####
dnl ## Use the same C++ compiler used in R, from 'Writing R Extensions'
dnl ## This is unnecessary because configure only concerns C
dnl CXX=`"${R_HOME}/bin/R" CMD config CXX`
dnl if test -z "$CXX"; then
dnl   AC_MSG_ERROR([No C++ compiler is available])
dnl fi
dnl CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
dnl CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
dnl dnl AC_LANG(C++)


#####
## Get host system type for HAVE_EXTENDED_PRECISION_REGISTERS
AC_CANONICAL_HOST

AC_LANG(C)

AC_PROG_CC
dnl AC_PROG_CPP

AC_SUBST(CFLAGS)
AC_SUBST(LIBM)


#####
## Check availability of inline;
## HAVE_INLINE required to define gsl_poly_eval()
## From lines 135-165 in configure.ac from GSL
if test "$ac_cv_c_inline" != no ; then 
dnl Check for "extern inline", using a modified version of the test
dnl for AC_C_INLINE from acspecific.mt
dnl
   AC_CACHE_CHECK([for GNU-style extern inline], ac_cv_c_extern_inline,
   [ac_cv_c_extern_inline=no
   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern $ac_cv_c_inline double foo(double x);
   extern $ac_cv_c_inline double foo(double x) { return x + 1.0 ; } ;
   double foo (double x) { return x + 1.0 ; };]], [[  foo(1.0)  ]])],[ac_cv_c_extern_inline="yes"],[])
   ])

   if test "$ac_cv_c_extern_inline" != no ; then
      AC_DEFINE(HAVE_INLINE,[1],[Define if you have inline])
   else
      AC_CACHE_CHECK([for C99-style inline], ac_cv_c_c99inline,
      [ac_cv_c_c99inline=no
      dnl next line is a necessary condition
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[extern inline void* foo() { foo(); return &foo ; };]], 
      [[  return foo() != 0 ]])],[ac_cv_c_c99inline="yes"],[])
      dnl but not sufficient, extern must work but inline on its own should not
      if test "$ac_cv_c_c99inline" != no ; then
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[inline void* foo() { foo(); return &foo ; };]], 
            [[  return foo() != 0 ]])],[],ac_cv_c_c99inline="no")
      fi
      ])
      if test "$ac_cv_c_c99inline" != no ; then
         AC_DEFINE(HAVE_INLINE,[1],[Define if you have inline])
         AC_DEFINE(HAVE_C99_INLINE,[1],[Define if you have inline with C99 behavior])
      fi
   fi
fi


#####
## Use libm; necessary to compile the check code for
## IEEE comparisons (with exp()) below
## From lines 216-220 in configure.ac from GSL
dnl Use alternate libm if specified by user

if test "x$LIBS" = "x" ; then
  AC_CHECK_LIB(m, cos)
fi


#####
## Check for HAVE_DECL_ISFINITE, HAVE_DECL_FINITE, and HAVE_DECL_ISNAN;
## HAVE_DECL_ISFINITE, HAVE_DECL_FINITE or HAVE_IEEE_COMPARISONS required
## to define gsl_finite();
## HAVE_DECL_ISNAN or HAVE_IEEE_COMPARISONS required to define gsl_isnan()
## From lines 236-241 in configure.ac from GSL
AC_CHECK_DECLS(isfinite,,,[#include <math.h>])
AC_CHECK_DECLS(finite,,,[#include <math.h>
#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif])
AC_CHECK_DECLS(isnan,,,[#include <math.h>])


#####
## HAVE_EXTENDED_PRECISION_REGISTERS influences the definition of
## GSL_COERCE_DBL(), which is used in gsl_sf_multiply_e()
## From lines 269-300 in configure.ac from GSL
AC_CACHE_CHECK([for extended floating point registers],ac_cv_c_extended_fp,
[case "$host" in
    *sparc*-*-*)
        ac_cv_c_extended_fp=no
        ;;     
    *powerpc*-*-*)
        ac_cv_c_extended_fp=no
        ;;      
    *hppa*-*-*)
        ac_cv_c_extended_fp=no
        ;;      
    *alpha*-*-*)
        ac_cv_c_extended_fp=no
        ;;      
    *68k*-*-*)
        ac_cv_c_extended_fp=yes
        ;;      
    *86-*-*)
        ac_cv_c_extended_fp=yes
        ;;      
    x86_64-*-*)
        ac_cv_c_extended_fp=yes
        ;;      
    *) 
        ac_cv_c_extended_fp=unknown
        ;;
esac
])

if test $ac_cv_c_extended_fp != "no" ; then
    AC_DEFINE(HAVE_EXTENDED_PRECISION_REGISTERS,1,[Defined on architectures with excess floating-point precision])
fi


#####
## Check for HAVE_IEEE_COMPARISONS; see above for where this is used
## From lines 441-458 in configure.ac from GSL
dnl Check IEEE comparisons, whether "x != x" is true for NaNs
dnl
AC_CACHE_CHECK([for IEEE comparisons], ac_cv_c_ieee_comparisons,
[AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <math.h>
int main (void) 
{ 
   int status; double inf, nan;
   inf = exp(1.0e10);
   nan = inf / inf ;
   status = (nan == nan);
   exit (status);
}]])],[ac_cv_c_ieee_comparisons="yes"],[ac_cv_c_ieee_comparisons="no"],[ac_cv_c_ieee_comparisons="yes"])
])

if test "$ac_cv_c_ieee_comparisons" != no ; then
  AC_DEFINE(HAVE_IEEE_COMPARISONS,1,[Define this if IEEE comparisons work correctly (e.g. NaN != NaN)])
fi


dnl #####
dnl ## To define GSL_COERCE_DBL() in config.h
dnl ## From lines 507-525 in configure.ac from GSL
AH_BOTTOM([/* Define a rounding function which moves extended precision values
   out of registers and rounds them to double-precision. This should
   be used *sparingly*, in places where it is necessary to keep
   double-precision rounding for critical expressions while running in
   extended precision. For example, the following code should ensure
   exact equality, even when extended precision registers are in use,

      double q = GSL_COERCE_DBL(3.0/7.0) ;
      if (q == GSL_COERCE_DBL(3.0/7.0)) { ... } ;

   It carries a penalty even when the program is running in double
   precision mode unless you compile a separate version of the
   library with HAVE_EXTENDED_PRECISION_REGISTERS turned off. */

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif])


dnl #####
dnl ## To define RETURN_IF_NULL() in config.h
dnl ## From lines 587-588 in configure.ac from GSL
AH_BOTTOM([#define RETURN_IF_NULL(x) if (!x) { return ; }
])


dnl #####
dnl ## To add user-defined methods in Eigen classes in config.h
AH_TOP([// To add user-defined methods in Eigen classes
#ifndef EIGEN_ARRAYBASE_PLUGIN
#define EIGEN_ARRAYBASE_PLUGIN "ArrayBaseAddons.h"
#endif
])


AC_CONFIG_HEADERS([src/config.h])

AC_OUTPUT
