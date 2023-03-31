This is a resubmission of a recently (2023-03-31) archived package.


## Comments to previous CRAN check results

* Version: 1.0.0
  Check: installed package size
  Result: NOTE
       installed size is 108.5Mb
       sub-directories of 1Mb or more:
       libs 107.5Mb
  Flavors: r-release-macos-arm64, r-release-macos-x86_64, r-oldrel-macos-arm64,
  r-oldrel-macos-x86_64, r-oldrel-windows-ix86+x86_64

- Reply: This large installation size is not due to error, but from
  a large number of C++ functions involved.  Please see further comments below.


### clang-ASAN
* checking whether package ‘qfratio’ can be installed ... [1094s/590s] WARNING
  Found the following significant warnings:
    :25: warning: 'result_of<Eigen::internal::scalar_product_op<double> (const double &, const double &)>' is deprecated [-Wdeprecated-declarations]
  See ‘/data/gannet/ripley/R/packages/tests-clang-SAN/qfratio.Rcheck/00install.out’ for details.

- Reply: I did some experiments over this deprecation warning, and found
  that this is from the combination of clang, libc++ (used in clang-ASAN), 
  C++17 (default in R-devel), and 'Eigen' 3.3 (included in 'RcppEigen'
  0.3.3.9.3).
  This could be cleared by specifying C++14 as the standard, but this is not
  done for portability.  Hopefully this will be solved by a future update of
  'RcppEigen' coming with the current version of 'Eigen' (3.4).
  In any case, this warning does not seem to happen in the latest development
  version of R (2023-03-30 r84127).


* checking examples ... ERROR
  Running examples in ‘qfratio-Ex.R’ failed
  The error most likely occurred in:
  ... #15 0x7f59b428b705 in ... qfratio/src/dk_funs.cpp:2062:44 ...

- Reply: The problem has been spotted and fixed in src/dk_funs.cpp.


### gcc-ASAN, valgrind

- These detected memory leak at the same point as above, which has been fixed.


## R CMD check results

0 errors | 0 warnings | 2 notes

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Junya Watanabe <jw2098@cam.ac.uk>’

New submission

Package was archived on CRAN

CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2023-03-31 as issues were not corrected
    in time.


* checking installed package size ... NOTE
  installed size is 105.3b
  sub-directories of 1Mb or more:
    libs  104.4Mb

- This package involves many C++ ('Eigen') functions, which are essential to
  the package functionality. This typically causes a large installation size
  in many environments (apparently except for R >= 4.2 on Windows).
  If possible, please ignore the check for installation size.
