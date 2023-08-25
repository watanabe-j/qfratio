// This file is taken from GSL version 2.7.1 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.

#include "gsl_integration.h" // added for qfratio

static inline void
reset_nrmax (gsl_integration_workspace * workspace);

static inline void
reset_nrmax (gsl_integration_workspace * workspace)
{
  workspace->nrmax = 0;
  workspace->i = workspace->order[0] ;
}
