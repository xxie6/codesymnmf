% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codesymnmf.R
\name{codesymnmf}
\alias{codesymnmf}
\title{This is an implementation of the coordinate descent symNMF method from
A. Vandaele, N. Gillis, Q. Lei, K. Zhong, I. Dhillon in R. Given an n-by-n
symmetric nonnegative matrix A and a factorization rank r, it outputs an
n-by-r nonnegative matrix H which solves min_{H >0} \|A - HH'\|_F^2}
\usage{
codesymnmf(A, r, init_type = "zeros", maxiter = 100, rand_permutation = FALSE)
}
\arguments{
\item{A}{An n-by-n symmetric matrix}

\item{r}{The number of columns in H}

\item{init_type}{The type of initialization to use for H.}

\item{maxiter}{The maximum number of iterations}

\item{rand_permutation}{A True/False value indicating whether to randomly
permute the columns of H when running coordinate descent}
}
\value{
A list containing the estimate for H and the objective function value
}
\description{
This is an implementation of the coordinate descent symNMF method from
A. Vandaele, N. Gillis, Q. Lei, K. Zhong, I. Dhillon in R. Given an n-by-n
symmetric nonnegative matrix A and a factorization rank r, it outputs an
n-by-r nonnegative matrix H which solves min_{H >0} \|A - HH'\|_F^2
}
