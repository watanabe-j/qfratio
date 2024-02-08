## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... [10s/32s] NOTE
Maintainer: ‘Junya Watanabe <Junya.Watanabe@uab.cat>’

New maintainer:
  Junya Watanabe <Junya.Watanabe@uab.cat>
Old maintainer(s):
  Junya Watanabe <jw2098@cam.ac.uk>

- The maintainer email address has changed because he changed affiliations.
  A confirmation will be sent from the old maintainer email.

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
