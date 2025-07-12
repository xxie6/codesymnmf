#' The initialization function for CoDesymNMF
#'
#' @param A An n-by-n symmetric matrix
#' @param r The number of columns in H
#' @param init_type The type of initialization to use for H.
#' If \code{init_type = 'zeros'}, then H is set to an n-by-r matrix of zeros.
#' Otherwise, H is set to an n-by-r matrix of entries drawn from standard uniform.
#'
#' @return Return an initial value for H
#'
#' @importFrom stats runif
#'
#' @export
#'

codesymnmf_init <- function(A, r, init_type){
  n <- ncol(A)
  if(init_type == 'zeros'){
    H <- matrix(rep(0,n*r), ncol = r)
  }else{
    H <- matrix(runif(n*r), ncol = r)
  }
}
