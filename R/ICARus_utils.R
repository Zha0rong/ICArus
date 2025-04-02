#' fasterICA_def
#' @description A faster implementation of fastICA_def from fastICA package
#' @importFrom matrixStats rowMeans2
#' @importFrom Rfast mat.mult
fasterICA_def <- function (X, n.comp, tol, fun, alpha, maxit, w.init){
  p <- ncol(X)
  W <- matrix(0, n.comp, n.comp)
  for (i in 1:n.comp) {
    w <- matrix(w.init[i,], n.comp, 1)
    if (i > 1) {
      t <- w
      t[1:length(t)] <- 0
      for (u in 1:(i - 1)) {
        k <- sum2(w * W[u, ])
        t <- t + k * W[u, ]
      }
      w <- w - t
    }
    w <- w/sqrt(sum2(w^2))
    lim <- rep(1000, maxit)
    it <- 1
    if (fun == "logcosh") {
      while (lim[it] > tol && it < maxit) {
        wx <- Rfast::mat.mult(Rfast::transpose(w), X)
        gwx <- tanh(alpha * wx)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X * gwx
        v1 <- matrixStats::rowMeans2(xgwx)
        g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
        v2 <- mean(g.wx) * w
        w1 <- v1 - v2
        w1 <- matrix(w1, n.comp, 1)
        it <- it + 1
        if (i > 1) {
          t <- w1
          t[1:length(t)] <- 0
          for (u in 1:(i - 1)) {
            k <- sum(w1 * W[u, ])
            t <- t + k * W[u, ]
          }
          w1 <- w1 - t
        }
        w1 <- w1/sqrt(sum(w1^2))
        lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
        
        w <- matrix(w1, n.comp, 1)
      }
    }
    if (fun == "exp") {
      while (lim[it] > tol && it < maxit) {
        wx <- Rfast::mat.mult(Rfast::transpose(w), X)
        gwx <- wx * exp(-(wx^2)/2)
        gwx <- matrix(gwx, n.comp, p, byrow = TRUE)
        xgwx <- X * gwx
        v1 <- matrixStats::rowMeans2(xgwx)
        g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
        v2 <- mean(g.wx) * w
        w1 <- v1 - v2
        w1 <- matrix(w1, n.comp, 1)
        it <- it + 1
        if (i > 1) {
          t <- w1
          t[1:length(t)] <- 0
          for (u in 1:(i - 1)) {
            k <- sum(w1 * W[u, ])
            t <- t + k * W[u, ]
          }
          w1 <- w1 - t
        }
        w1 <- w1/sqrt(sum(w1^2))
        lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
        w <- matrix(w1, n.comp, 1)
      }
    }
    W[i, ] <- w
  }
  W
}
#' fasterICA_par
#' @description A faster implementation of fastICA_par from fastICA package
#' @importFrom matrixStats rowMeans2
#' @importFrom Rfast mat.mult
fasterICA_par <- function (X, n.comp, tol, fun, alpha, maxit, verbose, w.init) {
  Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
  p <- ncol(X)
  W <- w.init
  sW <- rsvd::rsvd(W)
  W <- Rfast::mat.mult(sW$u, Rfast::mat.mult(Diag(1/sW$d), Rfast::mat.mult(Rfast::transpose(sW$u), W)))
  W1 <- W
  lim <- rep(1000, maxit)
  it <- 1
  if (fun == "logcosh") {
    while (lim[it] > tol && it < maxit) {
      wx <- Rfast::mat.mult(W, X)
      gwx <- tanh(alpha * wx)
      v1 <- Rfast::mat.mult(gwx, Rfast::transpose(X)/p)
      g.wx <- alpha * (1 - (gwx)^2)
      v2 <- Rfast::mat.mult(Diag(matrixStats::rowMeans2(g.wx)), W)
      W1 <- v1 - v2
      sW1 <- rsvd::rsvd(W1)
      W1 <- mat.mult(sW1$u, mat.mult(Diag(1/sW1$d), Rfast::mat.mult(Rfast::transpose(sW1$u),W1)))
      lim[it + 1] <- max(Mod(Mod(diag(Rfast::mat.mult(W1, Rfast::transpose(W)))) - 1))
      W <- W1
      it <- it + 1
    }
  }
  if (fun == "exp") {
    
    while (lim[it] > tol && it < maxit) {
      wx <- Rfast::mat.mult(W, X)
      gwx <- wx * exp(-(wx^2)/2)
      v1 <- Rfast::mat.mult(gwx, Rfast::transpose(X)/p)
      g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
      v2 <- Rfast::mat.mult(Diag(matrixStats::rowMeans2(g.wx)), W)
      W1 <- v1 - v2
      sW1 <- rsvd::rsvd(W1)
      W1 <- Rfast::mat.mult(sW1$u, mat.mult(Diag(1/sW1$d), Rfast::mat.mult(Rfast::transpose(sW1$u),W1)))
      lim[it + 1] <- max(Mod(Mod(diag(Rfast::mat.mult(W1, Rfast::transpose(W)))) - 1))
      W <- W1
      it <- it + 1
    }
  }
  W
}

#' faster_ICA_whitening
#' @description This function generates the whitened matrix for fasterICA
#' @importFrom Rfast mat.mult
#' @import rsvd
faster_ICA_whitening <- function (X) {
  require(Rfast)
  require(rsvd)
  n <- nrow(X)
  p <- ncol(X)
  X = Rfast::transpose(X)
  
  
  V <- Rfast::mat.mult(X, Rfast::transpose(X)/n)
  
  s <- rsvd::rsvd(V)
  D <- diag(c(1/sqrt(s$d)))
  
  K <- as.matrix(Rfast::mat.mult(D , Rfast::transpose(s$u)))
  return(list(X=X,K=K,n=n,p=p))
}

#' faster_ICA
#' @description A faster implementation of fastICA from fastICA package
#' @importFrom matrixStats rowMeans2
#' @importFrom Rfast mat.mult
#' @import Rcpp
faster_ICA <- function (whitening_list,n.comp, alg.typ = c("parallel","deflation"),
                        fun = c("logcosh", "exp"),row.norm = FALSE, maxit = 200, tol = 1e-04,alpha=1,
                        w.init=NULL,method='R') {
  alg.typ <- match.arg(alg.typ)
  fun <- match.arg(fun)
  
  if(is.null(w.init)){
    w.init <- matrix(rnorm(n.comp^2),n.comp,n.comp)
  }
  X=whitening_list$X
  K=whitening_list$K
  n=whitening_list$n
  p=whitening_list$p
  X1=whitening_list$X1
  if (method=='R'){
    if (alg.typ == "deflation"){
      a <- fasterICA_def(X1, n.comp, tol = tol, fun = fun,
                         alpha = alpha, maxit = maxit, w.init = w.init)
    }
    else if (alg.typ == "parallel"){
      a <- fasterICA_par(X1, n.comp, tol = tol, fun = fun,
                         alpha = alpha, maxit = maxit, verbose = verbose, w.init = w.init)
    }
    w <- Rfast::mat.mult(a , K)
    S <- Rfast::mat.mult(w , X)
    A <- Rfast::mat.mult(Rfast::transpose(w), solve(mat.mult(w, Rfast::transpose(w))))
  }
  
  return(list(X = Rfast::transpose(X), K = t(K), W = t(a), A = t(A), S = Rfast::transpose(S)))
}










#' ParaICA
#' @description This function runs ICA using multiple threads.
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import truncnorm
#' @import Rfast
#' @import snow

ParaICA <- function(CountMatrix,numberofcomponents,iteration,numberofcores=2,scale=T,...) {
  #faster_whiten=faster_ICA_whitening(CountMatrix)
  #faster_whiten$K <- matrix(faster_whiten$K[1:numberofcomponents, ], numberofcomponents, faster_whiten$p)
  
  #faster_whiten$X1 <- Rfast::mat.mult(faster_whiten$K, faster_whiten$X)
  
  cl <- snow::makeCluster(numberofcores)
  doSNOW::registerDoSNOW(cl)
  #pb <- txtProgressBar(max = iteration, style = 3)
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress = progress)
  x=foreach::foreach(i=seq(1,iteration),.packages=c('fastICA','Rfast')
                    # ,.options.snow = opts
                     ) %dopar% {
    resICA=fastICA(CountMatrix,n.comp = numberofcomponents,row.norm = scale,method='C',...)
    Affiliation.Matrix=(as.matrix(resICA$A))
    Signature.Matrix=(as.matrix(resICA$S))
    rownames(Affiliation.Matrix)=paste0('n.',seq(1,nrow(Affiliation.Matrix)))
    colnames(Affiliation.Matrix)=colnames(CountMatrix)
    colnames(Signature.Matrix)=paste0('n','.',seq(1,ncol((Signature.Matrix))))
    rownames(Signature.Matrix)=rownames(CountMatrix)
    Affiliation.Matrix=t(Affiliation.Matrix)
    Results=list()
    Results[['Affiliation.Matrix']]=Affiliation.Matrix
    Results[['Signature.Matrix']]=Signature.Matrix
    return(Results)
  }
  #close(pb)
  snow::stopCluster(cl)
  Affiliation.Matrix=list()
  Signature.Matrix=list()
  for (i in 1:length(x)) {
    tmp=x[[i]][['Affiliation.Matrix']]
    colnames(tmp)=paste0('iteration.',i,'_',colnames(tmp))
    Affiliation.Matrix[[i]]=tmp
    tmp=x[[i]][['Signature.Matrix']]
    colnames(tmp)=paste0('iteration.',i,'_',colnames(tmp))
    Signature.Matrix[[i]]=tmp
  }
  Affiliation.Matrix=do.call(cbind,Affiliation.Matrix)
  Signature.Matrix=do.call(cbind,Signature.Matrix)
  return(list('Signature.Matrix'=Signature.Matrix,
              'Affiliation.Matrix'=Affiliation.Matrix))
}

#' Cluster_Stability_Calculation
#' @description This function generates the stability index for ICA
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import snow
Cluster_Stability_Calculation <- function(Correlation_Matrix,Clustering_identity,numberofcores=2) {
  Correlation_Matrix=as.matrix(Correlation_Matrix)
  Correlation_Matrix=apply(Correlation_Matrix,c(1,2),as.numeric)
  names=names(Clustering_identity)
  Clustering_identity=as.character(Clustering_identity)
  names(Clustering_identity)=names
  Clusters=unique(Clustering_identity)
  cl <- snow::makeCluster(numberofcores)
  registerDoSNOW(cl)
  
  #pb <- txtProgressBar(max = length(Clusters), style = 3)
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress = progress)
  x=foreach(i=seq(1,length(Clusters)),.packages=c('fastICA','JADE')
            #, .options.snow = opts
            ,.inorder = T) %dopar% {
    cluster=Clusters[i]
    Ck=length((Clustering_identity[Clustering_identity==cluster]))
    krij=as.numeric(Correlation_Matrix[names(Clustering_identity[Clustering_identity==cluster]),
                                       names(Clustering_identity[Clustering_identity==cluster])])
    Cl=length((Clustering_identity[Clustering_identity!=cluster]))
    lrij=as.numeric(Correlation_Matrix[names(Clustering_identity[Clustering_identity==cluster]),
                                       names(Clustering_identity[Clustering_identity!=cluster])])
    stability_index=((1/(Ck^2))*(sum(krij)))-((1/(Ck*Cl))*(sum(lrij)))
    return(stability_index)
    
  }
  #close(pb)
  snow::stopCluster(cl)
  x=unlist(x)
  stability_indices=x
  names(stability_indices)=Clusters
  return(stability_indices)
}


#' Direction_correction
#' @description This function is used to correct the direction of the signature from different runs.
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import doSNOW
Direction_correction = function (Signature.Matrix,Affiliation.Matrix,Group) {
  Reference=unique(Group)[1]
  Reference.Matrix.A=Affiliation.Matrix[,names(Group)[Group==Reference]]
  Reference.Matrix.S=Signature.Matrix[,names(Group)[Group==Reference]]
  
  Results.A=list()
  Results.S=list()
  for (i in unique(Group)[unique(Group)!=Reference]) {
    tested.matrix.A=Affiliation.Matrix[,c(names(Group)[Group==i])]
    tested.matrix.S=Signature.Matrix[,c(names(Group)[Group==i])]
    
    correlation=cor(as.matrix(cbind(tested.matrix.S,Reference.Matrix.S)))
    correlation=correlation[names(Group)[Group==Reference],c(names(Group)[Group==i])]
    for (j in names(Group)[Group==i]) {
      tested.matrix.A[,j]=tested.matrix.A[,j]*sign(correlation[,j][abs(correlation[,j])==max(abs(correlation[,j]))])
      tested.matrix.S[,j]=tested.matrix.S[,j]*sign(correlation[,j][abs(correlation[,j])==max(abs(correlation[,j]))])
      
    }
    Results.A[[i]]=tested.matrix.A
    Results.S[[i]]=tested.matrix.S
    
  }
  Results.A=do.call(cbind,Results.A)
  Results.S=do.call(cbind,Results.S)
  
  Results.A=cbind(Reference.Matrix.A,Results.A)
  Results.S=cbind(Reference.Matrix.S,Results.S)
  
  return(list(Results.A=Results.A,
              Results.S=Results.S))
}






#' kneedle
#' @description The kneedle algorithm implemented by Emmet Tam. Pulled directly from https://github.com/etam4260/kneedle.
#' @param x A vector of x coordinates.
#' @param y A vector of y coordinates.
#' @param sensitivity How sensitive should knee detection be? Defaults to 1. Can be a float/double.
#' @import quantmod
#' @returns This function returns a x,y coordinate pair that corresponds to the knee point

kneedle <- function(x, y, decreasing, concave, sensitivity = 1) {
  # Make sure inputs are correct
  if(length(x) == 0 || length(y) == 0) {
    stop("Make sure size of both inputs x and y are greater than 0")
  }
  if(typeof(x) == "list"|| typeof(y) == "list" || is.data.frame(x) ||
     is.data.frame(y) || is.array(x) || is.array(y) || is.matrix(x) || is.matrix(y)) {
    stop("Make sure both inputs x and y are vectors")
  }
  if(length(x) != length(y)) {
    stop("Make sure size of both inputs x and y are equal")
  }
  
  data <- matrix(unlist(list(x, y)), ncol = 2)
  # This decreasing value has nothing to do with the user inputted value.
  data <- data[order(data[,1], decreasing = FALSE), ]
  
  
  # If both decreasing and concave are not specified, then this algorithm will
  # take a guess at those parameters instead of defaulting to certain values.
  # One method is to take the derivative of the starting point to the ending point
  # from min x value to max x value. This should be similar to taking the average of all
  # derivatives from xi to xi+1 where i = 1 to i = length(xvalues)
  if(missing(decreasing)) {
    # Increasing discrete data
    if( (data[(nrow(data)), 2] - data[1, 2]) >= 0 ) {
      decreasing = FALSE
      # Decreasing discrete data
    } else {
      decreasing = TRUE
    }
  }
  
  # To determine concavity we need to look at the second derivative of the
  # entire set of discrete data from xi to xi+1 where i = 1 to i = length(xvalues)
  # Taking the average of all the second derivatives, if greater or equal to 0
  # then concave up. If less than 0 then concave down.
  if(missing(concave)) {
    secondderiv <- diff(diff(data[, 2]) / diff(data[, 1]))
    if(mean(secondderiv) > 0) {
      concave = TRUE
    } else {
      concave = FALSE
    }
  }
  
  maxy <- max(y)
  miny <- min(y)
  maxx <- max(x)
  minx <- min(x)
  data[ ,1] <- (data[, 1]- min(data[, 1]))/(max(data[ ,1])- min(data[, 1]))
  data[ ,2] <- (data[, 2]- min(data[, 2]))/(max(data[ ,2])- min(data[, 2]))
  
  if(concave && !decreasing) {
    differ <- abs(c(data[ ,2] - data[ ,1]))
  } else if(concave && decreasing) {
    differ <- abs(c(data[ ,2] - (1 - data[ ,1])))
  } else if(!concave && !decreasing) {
    differ <- abs(c(data[ ,2] - data[ ,1]))
  } else if(!concave && decreasing) {
    differ <- abs(c(data[ ,2] - (1 - data[ ,1])))
  }
  
  
  peak.indices <- findPeaks(differ) - 1
  
  data <- cbind(data, differ)
  
  diffx = diff(data[, 1])
  T.lm.x.s <- sensitivity * mean(diffx)
  knee = NULL
  
  for(i in 1:length(peak.indices)) {
    T <- data[peak.indices[i] ,3] - (T.lm.x.s)
    
    y.value <- data[peak.indices[i] ,3]
    
    for(j in (peak.indices[i]):if(i+1 < length(peak.indices)) peak.indices[i+1] else length(differ)) {
      if(differ[j] < T) {
        knee = peak.indices[i];
        break;
      }
    }
    if(!is.null(knee)) {
      break;
    }
  }
  
  # Returns the x,y coordinate values
  x <- ((maxx - minx) * (data[knee, 1])) + minx
  y <- ((maxy - miny) * (data[knee, 2])) + miny
  
  return(c(as.numeric(x),as.numeric(y)))
}

