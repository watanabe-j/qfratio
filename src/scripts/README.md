# To update GSL

Files from GSL are copied and modified with shell scripts:
1. copy files using `src/scripts/copy_gsl.sh`; list of files is `src/scripts/gsl_files.mkf`
2. modify copied files using `src/scripts/process_gsl.sh`

```
# usage (from package root):
src/scripts/copy_gsl.sh [source_dir] src/gsl src/scripts/gsl_files.mkf
src/scripts/process_gsl.sh src/gsl
```

The scripts were written with GSL 2.7.1.
Modification may be necessary for newer versions.
