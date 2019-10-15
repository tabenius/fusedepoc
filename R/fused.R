dyn.load('src/fusedepoc.so')
#dyn.load('src/splitbregman.so')

#' Penalty matrix
#'
#' This function allows you to create a penalty matrix
#' @param p The number of predictors
#' @param K The number of classes
#' @importFrom utils combn
#' @keywords penalty
#' @export 
penmatrix <- function(p,K) {
  cc <- combn(1:K,2)
  m <- dim(cc)[2]*p
  L <- array(0,dim=c(m,K*p))
  mm <- 0
  for (i in 1:dim(cc)[2]) {
    p1 <- cc[1,i]
    p2 <- cc[2,i]
    for (i in 1:p) {
      L[mm+i,(p1-1)*p+i] <- 1
      L[mm+i,(p2-1)*p+i] <- -1
    }
    mm <- mm + p
  }
  list(L=L,m=m)
}

