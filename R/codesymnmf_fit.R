#' Solves min_{H >0} \|A - HH'\|_F^2
#'
#' @param A An n-by-n symmetric matrix
#' @param H An initial value for the n-by-r matrix H
#' @param maxiter The maximum number of iterations
#' @param rand_permutation A True/False value indicating whether to randomly
#' permute the columns of H when running coordinate descent
#'
#' @return A list containing the estimate for H and the objective function value
#' @export
#'
codesymnmf_fit <- function(A, H, maxiter = 100, rand_permutation = FALSE){
  n <- nrow(H)
  r <- ncol(H)

  # scaling (if H is not the zero matrix)
  if (max(H) > 0){
    HHt <- tcrossprod(H)
    beta <- sum(A*HHt)/sum(HHt^2) #scaling factor
    H <- sqrt(beta)*H
  }

  iter <- 1
  x <- rep(0,4)
  # main loop
  while(iter <= maxiter){
    # randomly permute columns if rand_permutation is TRUE
    if (rand_permutation == TRUE){
      col_permutation <- sample(r, size = r, replace = FALSE)
      H <- H[,col_permutation]
    }

    # calculate residual
    R <- A - tcrossprod(H)
    for (k in 1:r){
      R <- R + tcrossprod(H[,k])
      diag_R <- diag(R)
      HtH <- crossprod(H[,k]) #scalar, equals sum(H[,k]^2)
      for (i in 1:n){
        # coefficients a and b of x^3 + ax + b
        HtH <- HtH - H[i,k]^2 #scalar
        a <- HtH - diag_R[i] #scalar
        b <- -(t(H[,k]) %*% R[,i] - H[i,k]*diag_R[i]) #scalar

        Delta <- 4*(a^3) + 27*(b^2) # note: delta can be negative!
        d <- 0.5 * (-b + sqrt(as.complex(Delta/27))) #might be complex

        if (Delta <= 0) {
          r3 <- 2*(abs(d)^(1/3))
          th3 <- atan2(Im(d),Re(d))/3
          x[2] <- r3*cos(th3)
          x[3] <- r3*cos(th3 + (2*pi/3))
          x[4] <- r3*cos(th3 + (4*pi/3))
          x <- x[x >= 0]
          ind <- which.min((x^4)/4 + c(a)*(x^2)/2 + c(b)*x)
          H[i,k] <- x[ind]
          HtH <- HtH + H[i,k]^2
        }else{
          d <- Re(d) # d is just a real number in this case
          z <- sign(d)*(abs(d))^(1/3)
          val <- z - (a/(3*z))
          if(((val^4)/4 + a*(val^2/2) + b*val < 0) & (val >= 0)){
            HtH <- HtH + val^2
            H[i,k] <- val
          }else{
            H[i,k] <- 0
          }
        }
      } # end of 1:n for loop
      R <- R - tcrossprod(H[,k])
    } # end of 1:r for loop
    iter <- iter + 1
  } # end of maxiter while loop
  # return output
  return(list(H=H, obj_func = sum((A - tcrossprod(H))^2)))
}
