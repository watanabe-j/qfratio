// To add user-defined methods in Eigen classes
#ifndef EIGEN_ARRAYBASE_PLUGIN
#define EIGEN_ARRAYBASE_PLUGIN "ArrayBaseAddons.h"
#endif

#define HAVE_DECL_ISFINITE 1
#define HAVE_DECL_FINITE 1
#define HAVE_DECL_ISNAN 1
#define HAVE_EXTENDED_PRECISION_REGISTERS 1
#define HAVE_IEEE_COMPARISONS 1
#define HAVE_INLINE 1
#define HAVE_C99_INLINE 1
#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif

#define RETURN_IF_NULL(x) if (!x) { return ; }
