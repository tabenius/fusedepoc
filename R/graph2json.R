# { description: "ov"
#   parameters: {lambda1: 0.1, lambda2: 0.04},
#   adjacencyMatrix: [[...]]
# }
#
examplejson <- function(N) {
  I <- 2
  J <- 3
  K <- 4
  dat <- runif(n=(I*J*N*N*K), min=-2, max=2)
  FF <- array(data=dat, dim=c(I,J,N,N,K))
  lambda1s <- seq(from=0.1, to=1, length.out=I)
  lambda2s <- seq(from=0.1, to=1, length.out=J)
  invisible(result2json3(letters[1:K], lambda1s, lambda2s, LETTERS[1:N], FF))
}
result2json2 <- function(netnames,lambda1s,lambda2s,species,FF) {
  #require('jsonlite')
  l <- list(descriptions=netnames, parameters=list(lambda1=lambda1s, lambda2=lambda2s), species=species, A=FF)
  invisible(toJSON(l))
}

#' Creates a JSON string from fused EPoC results usable in the
#' GUI
#'
#' @param result Fused EPoC result
#' @param species vector of names of species / genes
#' @param netnames vector of names of cancers
#' @return a JSON string comprising of all that is needed for
#'the GUI
#' @export
#' @import jsonlite
#' @examples
#'
#' D <- syntheticPerturbZap(1000,10)
#' lambda1 <- seq(from=0.01, to=1, length.out=10)
#' lambda2 <- 500*exp(-(1:150)/5)
#' result <- fusedepoc(Y=list(D$r1$Y, D$r2$Y), U=list(D$r2$U,
#'   D$r2$U), lambda1=lambda1, lambda2=lambda2, method='A')
#' 
#' result.json <- fusedepocToJSON(result, paste(1:10),
#'   c('type A','type B'))
#'
fusedepocToJSON <- function(result, species, netnames) {
  l <- result2json3(netnames, result$lambda1, result$lambda2, species, result$coefficients)
  invisible(toJSON(l))
}

result2json3 <- function(netnames,lambda1s,lambda2s,species,FF, thresh=1e-1) {
  #require('jsonlite')
  I <- length(lambda1s)
  J <- length(lambda2s)
	K <- length(netnames)
	N <- dim(FF)[1]
	sparse <- list()
	for(i in 1:I) {
		for(j in 1:J) {
			el <- list()
			for(n1 in 1:N) {
				for(n2 in 1:N) {
					if(any(abs(FF[n1,n2,,i,j]) > thresh)) {
						el[[length(el)+1]] <- list(s=n1,t=n2,w=FF[n1,n2,,i,j])
					}
				}
			}
			sparse[[length(sparse)+1]] <- list(lambda1ix=i, lambda2ix=j, edgelist=el)
		}
	}
  l <- list(descriptions=netnames, parameters=list(lambda1=lambda1s, lambda2=lambda2s), species=species, sparse=sparse)
 # invisible(toJSON(l))
}
result2json1 <- function(netnames,lambda1s,lambda2s,FF) {
  #require('jsonlite')
  I <- dim(FF)[1]
  J <- dim(FF)[2]
  K <- dim(FF)[5]
  if(length(netnames) != K) stop("n.o. netnames != K")
  if(length(lambda1s) != I) stop("n.o. lambda1s != I")
  if(length(lambda2s) != J) stop("n.o. lambda2s != J")

  iii <- 1
  l <- list()
  for (k in 1:K) {
    for (i in 1:I) {
      lambda1 <- lambda1s[i]
      for (j in 1:J) {
        lambda2 <- lambda2s[j]
        l[[iii]] <- list(description=netnames[k], paramaters=list(lambda1=lambda1, lambda2=lambda2), adjacencyMatrix=FF[i,j,,,k])
        iii <- iii + 1
      }
    }
  }
  invisible(toJSON(l, auto_unbox=TRUE))
}

adjacencyMatrix2json <- function(A) {
  #require('jsonlite')
  txt <- toJSON(A, auto_unbox=TRUE)
  invisible(txt)
}
# TODO: export several graphs, one for each cancer
#       export colors
graph2json <- function (g,filename=NULL,format="jsmodule") {
  #require('jsonlite')
  es <- edgeNames(g)
  nodes <- lapply(unique(sort(unlist(lapply(es, function(e) strsplit(e,"~"))))), function(n) {
    return (list(data=list(id=as.character(n))))
  })
  edges <- lapply(es,function(ee) { 
    e<-ee[[1]];
    ids <- strsplit(e,"~"); 
    src <- ids[[1]][1]
    tgt <- ids[[1]][2]
    ed <- edgeData(g,from=src,to=tgt)[[1]]
    print(ed)
    return (list(data=list(id=as.character(e[1]),source=src,target=tgt,color="#F22", w=ed$weight))) 
  })
  l <- c(nodes,edges)
  txt <- toJSON(l, auto_unbox=TRUE)
  if (length(filename) > 0) {
    fc <- file(filename)
    if (format == "js") {
      writeLines(paste("var cytoscapedata = ",txt,";", sep=""), con=fc)
    } else if (format == "jsmodule") {
      writeLines(paste("module.exports = ",txt,";", sep=""), con=fc)
    } else if (format == "json") {
      writeLines(txt, con=fc)
    } else {
      stop("unknown format \"",format,"\". Supported formats are \"json\", \"jsmodule\" and \"js\".")
    }
  } 
  return(invisible(txt))
}
#' Save a JSON string as a file
#'
#' @param json JSON string
#' @param filename filename
#' @param format jsmodule, json or js
#' @export
savejson <- function(json, filename=NULL, format="jsmodule") {
  txt <- json
  if (length(filename) > 0) {
    fc <- file(filename)
    if (format == "js") {
      writeLines(paste("var cytoscapedata = ",txt,";", sep=""), con=fc)
    } else if (format == "jsmodule") {
      writeLines(paste("module.exports = ",txt,";", sep=""), con=fc)
    } else if (format == "json") {
      writeLines(txt, con=fc)
    } else {
      stop("unknown format \"",format,"\". Supported formats are \"json\", \"jsmodule\" and \"js\".")
    }
    close(fc)
  }  else {
    stop("filename required")
  }
}
