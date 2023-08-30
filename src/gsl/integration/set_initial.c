// This file is taken from GSL version 2.7.1 and distributed as part of qfratio
// with modification, in accordance with the GNU General Public License
// version 3.  All modified lines are marked with comments.
// - 2023 Junya Watanabe

#include "gsl_integration.h" // added for qfratio

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
                         double result, double error);

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
                         double result, double error)
{
  workspace->size = 1;
  workspace->rlist[0] = result;
  workspace->elist[0] = error;
}
