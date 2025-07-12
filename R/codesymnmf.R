#' This is an implementation of the coordinate descent symNMF method from
#' A. Vandaele, N. Gillis, Q. Lei, K. Zhong, I. Dhillon in R. Given an n-by-n
#' symmetric nonnegative matrix A and a factorization rank r, it outputs an
#' n-by-r nonnegative matrix H which solves min_{H >0} \|A - HH'\|_F^2
#'
#' @param A An n-by-n symmetric matrix
#' @param r The number of columns in H
#' @param init_type The type of initialization to use for H.
#' @param maxiter The maximum number of iterations
#' @param rand_permutation A True/False value indicating whether to randomly
#' permute the columns of H when running coordinate descent
#'
#' @return A list containing the estimate for H and the objective function value
#' @export
#'
codesymNMF <- function(A, r, init_type = 'zeros', maxiter = 100, rand_permutation = FALSE){
  H_init <- codesymnmf_init(A, r, init_type)
  codesymnmf_results <- codesymnmf_fit(A, H_init, maxiter, rand_permutation)
  return(codesymnmf_results)
}
