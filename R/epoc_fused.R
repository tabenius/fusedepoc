#source('R/fused.R')

#' Fused EPoC 
#'
#' From several (K) data sets with data-levels mRNA and CNA
#' estimate a multi-graph representing gene regulatory networks with fusing showing similarities and differences between sets of cancer tumour samples
#' @param Y list of K mRNA expression level matrices
#' classes of size N x p
#' @param U list of K Copy Number Aberration (CNA) matrices
#' of size N x p
#' @param lambda1 LASSO penalty parameter lambda1 penalizes
#' the L1 norm of the resulting coefficients. This regulates
#' sparsity
#' @param lambda2 Fused LASSO penalty parameter lambda2
#' penalizes the L1 norm of the difference between
#' coefficients in the K different classes
#' @param method 'G' or 'A' as defined for the EPoC method.
#' In short, 'G' uses CNA as predictors while 'A' uses mRNA as
#' predictors.
#' @param trace The verbosity of information while running
#' @param optim Whether to use the optimized version or not
#' @import epoc
#' @import lassoshooting
#' @examples
#' \dontrun{lambda1 <- 0.01
#'lambda2 <- 500*exp(-(1:150)/5)
#'result <- fusedepoc(Y=list(Y1, Y2), U=list(U1, U2), 
#'  lambda1=lambda1,lambda2=lambda2, method='A')
#' }
#'
#' @export
#' @useDynLib fusedepoc
 
fusedepoc <- function(Y,U,lambda1,lambda2,method='G',trace=1,optim=TRUE) {
#  if (!require(epoc)) stop("require package epoc")
  cl <- match.call()
  epsilon <- 1e-6
  if (typeof(Y)!="list") stop("not a list")
  if (typeof(U)!="list") stop("not a list")
  K <- length(Y)
  cat(paste("K",K,sep=": "),"\n")
  psame <- length(unique(sapply(Y,function(x) dim(x)[2]))) == 1
  if (psame != TRUE) stop("not same p")
  p <- dim(Y[[1]])[2]
  cat(paste("p",p,sep=": "),"\n")
  Ns <- sapply(Y,function(x) dim(x)[1])
  print(Ns)
  hasU <- !is.null(U)
  YonU <- list()
  Yres <- list()
  d <- array(0,dim=c(p,K))
  muYres <- array(0,dim=c(p,K))
  muY <- array(0,dim=c(p,K))
  muU <- array(0,dim=c(p,K))
  for (k in 1:K) {
    # alloc
    YonU[[k]] <- array(0,dim=c(Ns[k],p))
    Yres[[k]] <- array(0,dim=c(Ns[k],p))
    # /alloc
    if (trace > 0) cat("Centering...")
    muY[,k] <- colMeans(Y[[k]])
    Y[[k]] <- Y[[k]] - muY[,k]
    if (hasU) {
      muU[,k] <- colMeans(U[[k]])
      U[[k]] <- U[[k]] - muU[,k]
    } else {
      U <- list(array(0,dim=c(N,p)))
    }
    if (trace > 0) cat("DONE\n")
    regf <- function(i) {
      return(coef(lsfit(U[[k]][,i],Y[[k]][,i],intercept=T))[2])
    }
    d[,k] <- sapply(1:p,regf)
    d[,k] <- pmax(d[,k],0)
    for (i in 1:p) YonU[[k]][,i] <- d[i,k]*U[[k]][,i]
    Yres[[k]] <- Y[[k]] - YonU[[k]]
    muYres[,k] <- colMeans(Yres[[k]])
    Yres[[k]] <- Yres[[k]] - muYres[,k]
  }
  if(method=='G') {
    pred <- U
    resp <- Yres
  } else {
    pred <- Y
    resp <- Yres
  }
  cat('Finding lambdamax\n')
  lambdamax <- array(0,dim=K)
  for (k in 1:K) {
    inorms <- epoc.lambdamax(pred[[k]],resp[[k]],getall=T)
    lambdamax[k] <- max(inorms)
    scal <- sqrt(lambdamax[k]) # *(Ns[k]/max(Ns))**(-0.25)
    pred[[k]] <- pred[[k]] / scal
    resp[[k]] <- resp[[k]] / scal
  }
#  lambda1 <- max(lambdamax) * lambda1

  cat('Creating penalty matrix\n')
  pen <- penmatrix(p,K)
  L <- pen$L
  m <- pen$m
  rm(pen)

  # stack       / X1  0   0  \ 
  #        X = |  0   X2  0   |
  #             \ 0   0   X3 /
  X <- NULL
  for (k in 1:K) {
    if (k==1) {
      OL <- NULL
      OR <- array(0,dim=c(Ns[k],p*(K-1)))
    } else {
      OL <- array(0,dim=c(Ns[k],p*(k-1)))
      OR <- array(0,dim=c(Ns[k],p*(K-k)))
    }
    XR <- cbind(OL,pred[[k]],OR)
    X <- rbind(X,XR)
  }
#  print(dim(L))
#  print(dim(X))

  B <- array(0,dim=c(p,p,K,length(lambda1),length(lambda2)))
  if (optim == TRUE) {
    xtx <- crossprod(X)
  }
  cat('Running solver for row \n')
  for (i in 1:p) {
    cat(paste(i,'\n',sep=''))
    X1 <- X
    X1[,i] <- 0
    y <- unlist(sapply(1:K, function(k) resp[[k]][,i]))
#    print(length(y))
    fail <- FALSE
    err <- ""
    N <- sum(Ns)

    lambdaElastic <- as.double(0)
    Blasso <- array(rnorm(p),dim=p*K)
    if (optim == TRUE) {
      for (k in 1:K) {
        X1L <- pred[[k]]
        X1L[,i] <- 0
        lfit <- lassoshooting(X1L,resp[[k]][,i],max(lambda1))
        Blasso[(p*(k-1)+1):(p*(k-1)+p)] <- lfit$coefficients
      }
#      filename <- paste(strftime(format="fusedump_%Y-%m-%d_%H_%M_%S",x=Sys.time()),"_i",i,".rda",sep="")
#      save(err, m, K, N, Ns, p, i, L, xtx, X1, X, y, lambda1, lambda2,file=filename,compress=TRUE)
      xtx1 <- xtx
      xtx1[i,] <- 0; xtx1[,i] <- 0
      tryCatch(
        res <- .External('fusedopt',dims=as.integer(c(m,K,sum(Ns),p*K,length(lambda1),length(lambda2),length(lambdaElastic))),trace=as.integer(0),maxit=as.integer(1000),L=as.double(L),XtX=as.double(xtx1),X=as.double(X1),y=as.double(y),lasso=as.double(lambda1),fused=as.double(lambda2), elastic=as.double(lambdaElastic),beta=as.double(array(Blasso,dim=p*K*length(lambda1)*length(lambda2)))),
        error=function (error) {
          fail <<- TRUE
          err <<- error
        })
    } else {
      tryCatch(
        res <- .C('fused',as.integer(c(m,K,sum(Ns),p*K,length(lambda1),length(lambda2))),as.integer(c(0,1000)),as.double(L),as.double(X1),as.double(y),as.double(lambda1),as.double(lambda2),as.double(array(Blasso,dim=p*K*length(lambda1)*length(lambda2)))),
        error=function (error) {
          fail <<- TRUE
          err <<- error
      })
    }
    if (fail) {
      filename <- paste(strftime(format="fusedfaildump_%Y-%m-%d_%H_%M_%S",x=Sys.time()),"_i",i,".rda",sep="")
      save(err, m, K, N, p, i, L, X1, y, lambda1, lambda2,file=filename,compress=TRUE)
      B[,i,,,] <- NA
      stop(err)
    } else {
      if (optim)
        B[,i,,,] <- res$coefficients
      else
        B[,i,,,] <- res[[8]]
    }
    B[i,i,,,] <- d[i,]
    B[,i,,,][abs(B[,i,,,]) < epsilon] <- 0
  }
  cat("Returning\n")
  ret <- list(call=cl,coefficients=B,lambdamax=lambdamax,d=d,Y.mean=muY,U.mean=muU,Yres.mean=muYres, lambda1=lambda1, lambda2=lambda2)
  class(ret) <- "FusedEPOC"
  return(ret)
}

#' Net size
#'
#' This function show the number of links for an individual network matrix
#' @param A a network matrix
#' @param thresh coefficients with a magnitude smaller than
#' threshold will not be counted
#' @keywords network links
#' @export
netsize <- function(A, thresh=1e-2) {
  p <- dim(A)[1]
  di <- seq(1,p*p,by=p+1)
  sum(abs(A[-di]) > thresh)
}
