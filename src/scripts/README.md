# To update GSL

Files from GSL are copied and modified with shell scripts:
1. copy files using `src/scripts/copy_gsl.sh`; list of files is `src/scripts/gsl_files.mkf`
2. modify copied files using `src/scripts/process_gsl.sh`

```
# usage (from package root):
src/scripts/copy_gsl.sh [source_dir] src/gsl src/scripts/gsl_files.mkf
src/scripts/process_gsl.sh src/gsl
```

The scripts were originally written with GSL 2.7.1
and later updated for GSL 2.8.
Modification may be necessary for newer versions.

## Notes for GSL 2.8 (2024.07.19)
Most files used remain unchanged from 2.7.1; changes were limited to:
* Addition of "#include <stdlib.h>" in exp.c, hyperg_U.c, poch.c
  (all in specfun/)
* Use of abs() or fabs((double)()) instead of fabs() for int parameters;
  cast_fabs_arg_double() in process_gsl.sh is no longer necessary
* Addition of Lebedev quadrature-related funs in integration/gsl_integration.h;
  to be commented out
