#' Get colored multi-graph from Fused EPoC results
#'
#' Given parameter indices for lambda1 and lambda2, get a
#' multi-graph from a result set of Fused EPoC
#' @param i lambda1 index (LASSO penalty)
#' @param j lambda2 index (Fused penalty)
#' @param A Fused EPoC result set or a coefficient array
#' @param labs The names of genes/species/nodes
#' @param threshIsLink threshold for abs(coefficient) to be
#' included as a link
#' @param threshIsSame threshold for difference between
#' coefficients to be fused
#' @import Rgraphviz
#' @import graph
#' @importFrom methods new
#' @export
getGraph <- function(i, j, A, labs=NULL, threshIsLink=1e-2, threshIsSame=1e-3) {
  if (class(A) == 'FusedEPOC' || class(A) == 'FusedEPOCstabilized')
    A <- A$coefficients

  if (i > dim(A)[4]) stop("i outside range")
  if (j == 0) {
    A <- A[,,,i,]
  } else {
    if (j > dim(A)[5]) stop("j outside range")
    A <- A[,,,i,j]
  }
  K <- dim(A)[3]
  p <- dim(A)[1]
  if (!is.null(labs))
    labsc <- labs
  else
    labsc <- paste(1:p)

  

  # diagonal self-edges are not showed here
  for (i in (1:K)) {
    diag(A[,,i]) <- 0
  }

  fno <- sapply(1:p, function(ii) any(abs(A[ii,,]) >= threshIsLink))
  tno <- sapply(1:p, function(ii) any(abs(A[,ii,]) >= threshIsLink))
  mask <- fno | tno

  nodes <- labsc[1:p][mask]
  gw <- new("graphNEL", nodes=nodes,edgemode="directed")
  splines <- NULL
  edgetype <- NULL
  edgecolor <- NULL
  thick <- NULL
  pwfact <- 1
  en <- NULL
  shape <- rep((foo="rectangle"),length(nodes))
  width <- rep((foo=1.6),length(nodes))
  totedges <- 0
  for (ii in (1:p)[fno]) {
    for (jj in (1:p)[tno]) {
      if (ii != jj && any(abs(A[ii,jj,]) > threshIsLink)) {
	#               from       to
	edgename <- paste(labsc[ii],"~",labsc[jj],sep="")
	en <- c(en,edgename)
	splines <- c(splines, splines="ortho")
	sameall <- all(apply(combn(1:K,2),2,function (netpair) {
            return (abs(A[ii,jj,netpair[1]] - A[ii,jj,netpair[2]]) < threshIsSame)
          }
	))
  ew <- 0
	if (sameall) {
	  edgecolor <- c(edgecolor,color="black",recursive=T)
	  ew <- abs(mean(A[ii,jj,]))
	} else if (sum(abs(A[ii,jj,]) > threshIsLink) == K) {
	  edgecolor <- c(edgecolor,color="orange",recursive=T)
	  ew <- abs(mean(A[ii,jj,]))
	} else {
	  if (abs(A[ii,jj,1]) > threshIsLink && abs(A[ii,jj,2]) > threshIsLink ) {
	    edgecolor <- c(edgecolor,color="yellow",recursive=T)
	    ew <- mean(abs(A[ii,jj,c(1,2)]))
	  } else if (abs(A[ii,jj,1]) > threshIsLink) {
	    edgecolor <- c(edgecolor,color="red",recursive=T)
	    ew <- abs(A[ii,jj,1])
	  } else if (abs(A[ii,jj,2]) > threshIsLink) {
	    edgecolor <- c(edgecolor,color="green",recursive=T)
	    ew <- abs(A[ii,jj,2])
	  } else if (abs(A[ii,jj,3]) > threshIsLink && abs(A[ii,jj,1]) > threshIsLink ) {
	    edgecolor <- c(edgecolor,color="purple",recursive=T)
	    ew <- mean(abs(A[ii,jj,c(1,3)]))
	  } else if (abs(A[ii,jj,3]) > threshIsLink) {
	    edgecolor <- c(edgecolor,color="blue",recursive=T)
	    ew <- mean(abs(A[ii,jj,c(1,3)]))
	  } else if (abs(A[ii,jj,3]) > threshIsLink && abs(A[ii,jj,2]) > threshIsLink ) {
	    edgecolor <- c(edgecolor,color="cyan",recursive=T)
      ew <- mean(abs(A[ii,jj,c(2,3)]))
	  } else {
      cat("problem with this edge (could be that more than 3 cancer types are unsupported:")
      print(edgename)
    }
	}
  print(c(edgename,ew))
  if (ew > 0) {
    ew <- ew * pwfact
    thick <- c(thick, ew)
    gw <- addEdge(labsc[ii], labsc[jj], gw, ew) #A[ii,jj,1])

    et <- ifelse(mean((A[ii,jj,])) > 0, "normal", "tee")
    edgetype <- c(edgetype, arrowType=et, recursive=T)
  } else {
      #print(edgename)
  }
      }
    }
    totedges <- totedges + 1
  }
  if (totedges == 0) {
    stop("no edges in this matrix")
  }
  names(shape) <- nodes
  names(width) <- nodes
  lw2 <- thick
  names(lw2) <- en
  names(splines) <- en
  names(edgecolor) <- en
  names(edgetype) <- en
  lwd <- 1 + 2*thick
  names(lwd) <- en
  l <- list(graph=gw,edgeAttrs=list(color=edgecolor, arrowhead=edgetype, lwd=lwd),nodeAttrs=list(width=width,shape=shape),mask=mask, en=en)
	l$ng <- layoutGraph(l$graph, edgeAttrs=l$edgeAttrs,nodeAttrs=l$nodeAttrs,name="gene regulatory network")
  return(l$ng)
}
getEdgeAttr <- function(graph, attribute, ID1, ID2) {
  idfunc <- function(x) paste0( sort(x), collapse="~" )
  all.ids <- sapply( AgEdge(graph), 
    function(e) idfunc( c( attr(e, "head"), attr( e, "tail" ))))
  sel.ids <- apply(cbind( ID1, ID2 ), 1, idfunc )
  if(!all(sel.ids %in% all.ids)) stop( "only existing edges, please" )
  sel <- match( sel.ids, all.ids )
  for(i in 1:length(sel)) {
    return (attr( attr( graph, "AgEdge" )[[ sel[i] ]], attribute ))
  }
}
# from https://stackoverflow.com/questions/29282522/setting-edge-width-in-rgraphviz
setEdgeAttr <- function( graph, attribute, value, ID1, ID2 ) {
  idfunc <- function(x) paste0( sort(x), collapse="~" )
  all.ids <- sapply( AgEdge(graph), 
    function(e) idfunc( c( attr(e, "head"), attr( e, "tail" ))))
  sel.ids <- apply(cbind( ID1, ID2 ), 1, idfunc )
  if(!all(sel.ids %in% all.ids)) stop( "only existing edges, please" )
  sel <- match( sel.ids, all.ids )
  for(i in 1:length(sel)) {
    attr( attr( graph, "AgEdge" )[[ sel[i] ]], attribute ) <- value[i]
  }
  return(graph)
}

