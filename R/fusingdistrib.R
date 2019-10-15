#' Show distribution of fusing over genes
#'
#' Plot a heatmap showing how fusing is distributed over genes
#' and lambda2. This can be used for a broad overview to identify genes 
#' behaving similar or different in two classes.
#' @param result Fused EPoC result set
#' @param ks the two classes to compare
#' @param l1 index of lambda1 vector
#' @import stats
#' @importFrom graphics image
#' @importFrom grDevices hcl.colors
#' @examples
#' D <- syntheticPerturbZap(1000,10)
#' lambda1 <- 0.01
#' lambda2 <- 500*exp(-(1:150)/5)
#' result <- fusedepoc(Y=list(D$r1$Y, D$r2$Y), U=list(D$r2$U,D$r2$U), 
#'   lambda1=lambda1,lambda2=lambda2, method='A')
#' H <- fuseplot(result, c(1,2), l1=1)
#' @export
fuseplot <- function(result, ks, l1=1) {
  nl2 <- length(result$lambda2)
  p <- dim(result$coefficients)[1]
  if (length(ks) != 2) stop("ks should be a vector of two classes (integers)")
  xs <- NULL
  ys <- NULL
  H <- array(dim=c(nl2,p))
  for (l2 in 1:nl2) {
    A <- result$coefficients[,,ks[1],l1,l2]
    B <- result$coefficients[,,ks[2],l1,l2]
    meandistance <- mean(abs(A-B))
    xs <- c(xs, result$lambda2[l2])
    ys <- c(ys, meandistance)
    dd <- sapply(1:p, function(i) mean(abs(A[i,] - B[i,])))
    H[l2,] <- dd
  }
  image(1:nl2, 1:p, H, xlab='fusing parameter index', ylab='node number', main=paste('Distribution of mean distance of link magnitude \nper gene vs lambda2 with lambda1 = ',result$lambda1[l1]), col=grDevices::hcl.colors(12, "YlOrRd", rev = FALSE))
  invisible(H)
}

