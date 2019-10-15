genDataFromMasks <- function(mask, mask2, N) {
  p <- dim(mask)[1]
  Y <- array(stats::rnorm(N*p), dim=c(N,p))
  U <- array(stats::rnorm(N*p), dim=c(N,p))
  for(i in 1:p) {
    for(j in 1:p) {
      # mRNA expression levels affected by other mRNA expression
      # levels
      if (mask[i,j] > 0) {
        Y[,i] <- Y[,i] + mask[i,j]*Y[,j] + stats::rnorm(N, sd=0.01)
      }
      # mRNA expression levels also affected by DNA CNA
      if (mask2[i,j] > 0) {
        Y[,i] <- Y[,i] + mask2[i,j]*U[,j] + stats::rnorm(N, sd=0.01)
      }
    }
  }
  list(Y=Y, U=U)
}
genData1 <- function(N,p) {
  sparsity = .2
  mask = rand(p)
  mask = mask * (mask < sparsity)

  sparsity = .2
  mask2 = rand(p)
  mask2 = mask2 * (mask2 < sparsity)
  diag(mask2) <- stats::runif(p) # self-CNA to mRNA

  r <- genDataFromMasks(mask, mask2, N)
  list(mask=mask, mask2=mask2, Y=r$Y,U=r$U)
}
genData2 <- function(N,p) {
  r <- genData1(N,p)
  maskOther <- r$mask
  maskOther2 <- r$mask2
  nodiffering <- p*p * 0.1
  # perturb 10% of the coefficients randomly
  ii <- round(stats::runif(nodiffering, min=1, max=p))
  jj <- round(stats::runif(nodiffering, min=1, max=p))
  for (k in 1:nodiffering) {
    old <- r$mask[ii[k], jj[k]]
    new <- old * stats::rnorm(1)
    maskOther[ii[k], jj[k]] <- new

    old <- r$mask2[ii[k], jj[k]]
    new <- old * stats::rnorm(1)
    maskOther2[ii[k], jj[k]] <- new
  }
  r2 <- genDataFromMasks(maskOther, maskOther2, N)
  list(r1=r, r2=list(mask=maskOther, mask2=maskOther2, Y=r2$Y, U=r2$U), perturb.ii=ii, perturb.jj=jj)
}

#' Generate synthetic data for Fused EPoC
#'
#' This function allows you to create randomized synthetic
#' data which can be used by Fused EPoC. Two networks are
#' generated where the second is a modified version of the
#' first with both perturbed and zapped links.
#'
#' @param N number of samples
#' @param p number of nodes = genes in the Gene Regulatory
#' Network
#' @import pracma
#' @import stats
#' @import graph
#' @examples
#' D <- syntheticPerturbZap(1000,10)
#' lambda1 <- seq(from=0.01, to=1, length.out=10)
#' lambda2 <- 500*exp(-(1:150)/5)
#' result <- fusedepoc(Y=list(D$r1$Y, D$r2$Y), U=list(D$r2$U,
#' D$r2$U), lambda1=lambda1, lambda2=lambda2, method='A')
#' gl <- getGraph(2,40,result)
#' \dontrun{Rgraphviz::plot(gl$ng)}
#' \dontrun{graph::plot(gl$graph)}
#' @export
syntheticPerturbZap <- function(N,p) {
#  require('pracma')
  r <- genData1(N,p)
  maskOther <- r$mask
  maskOther2 <- r$mask2
  nodiffering <- p*p * 0.1
  # perturb 10% of the coefficients randomly
  ii <- round(stats::runif(nodiffering, min=1, max=p))
  jj <- round(stats::runif(nodiffering, min=1, max=p))
  for (k in 1:nodiffering) {
    old <- r$mask[ii[k], jj[k]]
    new <- old * stats::rnorm(1)
    maskOther[ii[k], jj[k]] <- new

    old <- r$mask2[ii[k], jj[k]]
    new <- old * stats::rnorm(1)
    maskOther2[ii[k], jj[k]] <- new
  }
  # zap 10% of the coefficients randomly
  zap.ii <- round(stats::runif(nodiffering, min=1, max=p))
  zap.jj <- round(stats::runif(nodiffering, min=1, max=p))
  for (k in 1:nodiffering) {
    maskOther[ii[k], jj[k]] <- 0
    maskOther2[ii[k], jj[k]] <- 0
  }
  r2 <- genDataFromMasks(maskOther, maskOther2, N)
  list(r1=r, r2=list(mask=maskOther, mask2=maskOther2, Y=r2$Y, U=r2$U), zap.ii=zap.ii, zap.jj=zap.jj, perturb.ii=ii, perturb.jj=jj)
}


