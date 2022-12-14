% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dk_funs.R
\name{arl}
\alias{arl}
\alias{arl.default}
\alias{arl.matrix}
\title{Recursion for a_{r,l}}
\usage{
arl(L, ...)

\method{arl}{default}(L, D, m = 10L)

\method{arl}{matrix}(A, mu, m = 10L)
}
\arguments{
\item{L}{Eigenvalues of the argument matrix; vector of \eqn{\lambda_i}}

\item{...}{Additional arguments passed to internal methods}

\item{D}{Squared norm of the mean vector projected on the eigenvalues of
the argument matrix: vectors of \eqn{\delta_i}}

\item{m}{Scalar to specify the desired order}

\item{A}{Argument matrix.  Assumed to be symmetric.}

\item{mu}{Mean vector \eqn{\bm{\mu}} for \eqn{\mathbf{x}}}
}
\description{
\code{arl()} is an internal function to calculate \eqn{a_{r,l}} as defined
in Hillier et al. (2014; eq. 24), which is used in the calculation of
the moment of such a ratio of quadratic forms in normal variables
where the denominator matrix is identity.
}
\details{
This function implements the super-short recursion described in
Hillier et al. (2014  eqs. 31--32).
Note that \eqn{w_{r,i}} there should be understood as \eqn{w_{r,l,i}} with
the index \eqn{l} fixed for each \eqn{a_{r,l}}.

The \code{matrix} method just calculates \code{L} and \code{D} from
\code{A} and \code{mu} and passes them to the \code{default} method.
}
\seealso{
\code{\link{qfrm}}; this function is used in \code{qfrm_ApIq_int()}
(for noncentral cases only)
}
