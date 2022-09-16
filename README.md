
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qfratio: R Package for Moment of Ratios of Quadratic Forms

<!-- badges: start -->
<!-- badges: end -->

This package provides functions to evaluate moments of ratios of
quadratic forms in normal variables, specifically using recursive
algorithms developed by Chikuse, Hillier, Bao, Kan and colleagues.

There exists a couple of `Matlab` programs developed by Raymond Kan
(available from <https://www-2.rotman.utoronto.ca/~kan/>), but this `R`
package has been developed independently with added functionalities.
This has originally been developed for a biological application,
specifically for evaluating average evolvability measures in
quantitative genetics.

This project is under active development. Most planned functionalities
have been implemented, but have not been extensively documented. At
present, substantial restructuring might be possible.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("watanabe-j/qfratio")
```

This package has the following dependencies:

    Imports: Rcpp (>= 1.0.8.3)
    LinkingTo: Rcpp, RcppEigen
    Suggests: gsl, mvtnorm

## Examples

Here are some simple examples:

``` r
## Simple matrices
nv <- 4
A <- diag(1:nv)
B <- diag(sqrt(nv:1))

## Expectation of (x^T A x)^2 / (x^T x)^2 where x ~ N(0, I)
qfrm(A, p = 2)
#> 
#>  Moment of ratio of quadratic form
#> 
#> Moment = 6.6667
#> This value is exact

## Compare with Monte Carlo mean
mean(rqfr(1000, A = A, p = 2))
#> [1] 6.663473

## Expectation of (x^T A x)^1/2 / (x^T x)^1/2
(mom_A0.5 <- qfrm(A, p = 1/2))
#> 
#>  Moment of ratio of quadratic form
#> 
#> Moment = 1.5672, Error = -6.3358e-19
#> Possible range:
#>  1.567224 1.567224

## Monte Carlo mean
mean(rqfr(1000, A = A, p = 1/2))
#> [1] 1.55763

plot(mom_A0.5)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r

## Expectation of (x^T x) / (x^T A^-1 x)
##   = "average conditional evolvability"
(avr_cevoA <- qfrm(diag(nv), solve(A)))
#> 
#>  Moment of ratio of quadratic form
#> 
#> Moment = 2.1168, Error = 2.7756e-15
#> Possible range:
#>  2.11678 2.11678

mean(rqfr(1000, A = diag(nv), B = solve(A), p = 1))
#> [1] 2.124956
plot(avr_cevoA)
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" />

``` r

## Expectation of (x^T x)^2 / (x^T A x) (x^T A^-1 x) 
##   = "average autonomy"
(avr_autoA <- qfmrm(diag(nv), A, solve(A), p = 2, q = 1, r = 1))
#> 
#>  Moment of ratio of quadratic form
#> 
#> Moment = 0.84166
#> Error bound unavailable; recommended to inspect plot() of this object

mean(rqfmr(1000, A = diag(nv), B = A, D = solve(A), p = 2, q = 1, r = 1))
#> [1] 0.8387947
plot(avr_autoA)
```

<img src="man/figures/README-unnamed-chunk-3-3.png" width="100%" />

``` r

## Expectation of (x^T A B x) / ((x^T A^2 x) (x^T B^2 x))^1/2 
##   = "average vector correlation"
## whose Monte Carlo evaluation is called the "random skewers" analysis,
## while this is essentially an analytic solution (with slight truncation error)
(avr_rcorA <- qfmrm(crossprod(A, B), crossprod(A), crossprod(B), 
                    p = 1, q = 1/2, r = 1/2))
#> 
#>  Moment of ratio of quadratic form
#> 
#> Moment = 0.84622
#> Error bound unavailable; recommended to inspect plot() of this object

mean(rqfmr(1000, A = crossprod(A, B), B = crossprod(A), D = crossprod(B), 
           p = 1, q = 1/2, r = 1/2))
#> [1] 0.8459347
plot(avr_rcorA)
```

<img src="man/figures/README-unnamed-chunk-3-4.png" width="100%" />

``` r


## More complex (but arbitrary) example
## Expectation of (x^T A x)^2 / (x^T B x)^3 where x ~ N(mu, Sigma)
mu <- 1:nv / nv
Sigma <- diag(runif(nv) * 3)
(mom_A2B3 <- qfrm(A, B, p = 2, q = 3, mu = mu, Sigma = Sigma, 
                  m = 500, use_cpp = TRUE))
#> 
#>  Moment of ratio of quadratic form
#> 
#> Moment = 0.51095, Error = 0
#> Possible range:
#>  0.510947 0.510947
plot(mom_A2B3)
```

<img src="man/figures/README-unnamed-chunk-3-5.png" width="100%" />

## References

Bao, Y. & Kan, R. (2013). On the moments of ratios of quadratic forms in
normal random variables. *Journal of Multivariate Analysis*, **117**,
229–245. doi:
[10.1016/j.jmva.2013.03.002](https://doi.org/10.1016/j.jmva.2013.03.002).

Hillier, G., Kan, R, & Wang, X. (2009). Computationally efficient
recursions for top-order invariant polynomials with applications.
*Econometric Theory*, **25**, 211–242. doi:
[10.1017/S0266466608090075](https://doi.org/10.1017/S0266466608090075).

Hillier, G., Kan, R, & Wang, X. (2014). Generating functions and short
recursions, with applications to the moments of quadratic forms in doi:
[10.1017/S0266466613000364](https://doi.org/10.1017/S0266466613000364).

Smith, M. D. (1989). On the expectation of a ratio of quadratic forms in
normal variables. *Journal of Multivariate Analysis*, **31**, 244–257.
doi:
[10.1016/0047-259X(89)90065-1](https://doi.org/10.1016/0047-259X(89)90065-1)

Smith, M. D. (1993). Expectations of ratios of quadratic forms in normal
variables: evaluating some top-order invariant polynomials. *Australian
Journal of Statistics*, **35**, 271–282. doi:
[10.1111/j.1467-842X.1993.tb01335.x](https://doi.org/10.1111/j.1467-842X.1993.tb01335.x).
