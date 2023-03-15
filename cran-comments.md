This is a resubmission

## Reply to previous submission comments

* "Using foo:::f instead of foo::f allows access to unexported objects.
  This is generally not recommended, as the semantics of unexported objects may be changed by the package author in routine maintenance."
  Please omit one colon.
  Used ::: in documentation:
       man/KiK.Rd:
          qfratio:::KiK(S1)
       man/KiK.Rd:
          qfratio:::KiK(S2)
       man/KiK.Rd:
          qfratio:::KiK(S2, tol = 1e-20)

* You have examples for unexported functions.
  Please either omit these examples or export these functions.
  Examples for unexported function
    htil2_pj() in:
       KiK.Rd

- Reply: Both were caused by examples in the unexported function KiK().
  The example is now omitted.

## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Junya Watanabe <jw2098@cam.ac.uk>’

New submission


* checking installed package size ... NOTE
  installed size is 128.2Mb
  sub-directories of 1Mb or more:
    libs  127.3Mb

- This package involves many C++ ('Eigen') functions, which are essential to
  the package functionality. This typically causes a large installation size
  on Linux environments (but apparently not on Windows or macOS).
