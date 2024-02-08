## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... [6s/12s] Note_to_CRAN_maintainers
Maintainer: ‘Junya Watanabe <jw2098@cam.ac.uk>’

* checking installed package size ... NOTE
  installed size is 132.6Mb
  sub-directories of 1Mb or more:
    libs  131.3Mb

- This package involves many C++ ('Eigen') functions, which are essential to
  the package functionality. This typically causes a large installation size
  in many environments (apparently except for R >= 4.2 on Windows).
  If possible, please ignore the check for installation size.


## Reverse dependencies

Currently, there are no reverse dependencies on this package.
