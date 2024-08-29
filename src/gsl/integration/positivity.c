// This file is taken from GSL version 2.8 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.
// - 2023 Junya Watanabe

/* Compare the integral of f(x) with the integral of |f(x)|
   to determine if f(x) covers both positive and negative values */

#include <math.h> // added for qfratio
#include "../gsl_machine.h" // added for qfratio

static inline int
test_positivity (double result, double resabs);

static inline int
test_positivity (double result, double resabs)
{
  int status = (fabs (result) >= (1 - 50 * GSL_DBL_EPSILON) * resabs);

  return status;
}
