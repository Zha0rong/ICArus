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
        wx <- Rfast::mat.mult(t(w), X)
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
        wx <- Rfast::mat.mult(t(w), X)
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
  W <- Rfast::mat.mult(sW$u, Rfast::mat.mult(Diag(1/sW$d), Rfast::mat.mult(t(sW$u), W)))
  W1 <- W
  lim <- rep(1000, maxit)
  it <- 1
  if (fun == "logcosh") {
    while (lim[it] > tol && it < maxit) {
      wx <- Rfast::mat.mult(W, X)
      gwx <- tanh(alpha * wx)
      v1 <- Rfast::mat.mult(gwx, t(X)/p)
      g.wx <- alpha * (1 - (gwx)^2)
      v2 <- Rfast::mat.mult(Diag(matrixStats::rowMeans2(g.wx)), W)
      W1 <- v1 - v2
      sW1 <- rsvd::rsvd(W1)
      W1 <- mat.mult(sW1$u, mat.mult(Diag(1/sW1$d), Rfast::mat.mult(t(sW1$u),W1)))
      lim[it + 1] <- max(Mod(Mod(diag(Rfast::mat.mult(W1, t(W)))) - 1))
      W <- W1
      it <- it + 1
    }
  }
  if (fun == "exp") {
    
    while (lim[it] > tol && it < maxit) {
      wx <- Rfast::mat.mult(W, X)
      gwx <- wx * exp(-(wx^2)/2)
      v1 <- Rfast::mat.mult(gwx, t(X)/p)
      g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
      v2 <- Rfast::mat.mult(Diag(matrixStats::rowMeans2(g.wx)), W)
      W1 <- v1 - v2
      sW1 <- rsvd::rsvd(W1)
      W1 <- Rfast::mat.mult(sW1$u, mat.mult(Diag(1/sW1$d), Rfast::mat.mult(t(sW1$u),W1)))
      lim[it + 1] <- max(Mod(Mod(diag(Rfast::mat.mult(W1, t(W)))) - 1))
      W <- W1
      it <- it + 1
    }
  }
  W
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
    A <- Rfast::mat.mult(t(w), solve(mat.mult(w, t(w))))
  }
  
  return(list(X = t(X), K = t(K), W = t(a), A = t(A), S = t(S)))
  
}

#' ParaICA
#' @description This function runs ICA using multiple threads.
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import truncnorm
#' @import Rfast
#' @import snow

ParaICA <- function(CountMatrix,faster_whiten=NULL,numberofcomponents,iteration,numberofcores=2,...) {

  cl <- snow::makeCluster(numberofcores)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = iteration, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach::foreach(i=seq(1,iteration),.packages=c('fastICA','Rfast'),.export = c('faster_ICA','fasterICA_par','fasterICA_def'), .options.snow = opts) %dopar% {
    if (is.null(faster_whiten)){resICA=fastICA(CountMatrix,n.comp = numberofcomponents,...)}
    if (!is.null(faster_whiten)){resICA=faster_ICA(whitening_list = faster_whiten,n.comp = numberofcomponents,...)
    }
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
  close(pb)
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

#' faster_ICA_whitening
#' @description This function generates the whitened matrix for fasterICA
#' @importFrom Rfast mat.mult
#' @import rsvd
faster_ICA_whitening <- function (X) {
  require(Rfast)
  require(rsvd)
  n <- nrow(X)
  p <- ncol(X)
  X = t(X)
  
  
  V <- Rfast::mat.mult(X, t(X)/n)
  
  s <- rsvd::rsvd(V)
  D <- diag(c(1/sqrt(s$d)))
  
  K <- Rfast::mat.mult(D , t(s$u))
  return(list(X=X,K=K,n=n,p=p))
}

#' IR_Calculation
#' @description This function generates the IR index for ICA
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import snow

IR_Calculation <- function(Correlation_Matrix,Clustering_identity,numberofcores=2) {
  Correlation_Matrix=as.matrix(Correlation_Matrix)
  Correlation_Matrix=apply(Correlation_Matrix,c(1,2),as.numeric)
  names=names(Clustering_identity)
  Clustering_identity=as.character(Clustering_identity)
  names(Clustering_identity)=names
  Clusters=unique(Clustering_identity)
  stability_indices=c()
  cl <- snow::makeCluster(numberofcores)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = length(Clusters), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(1,length(Clusters)),.packages=c('fastICA'), .options.snow = opts) %dopar% {
    cluster=Clusters[i]
    Ck=length((Clustering_identity[Clustering_identity==cluster]))
    krij=as.numeric(Correlation_Matrix[names(Clustering_identity[Clustering_identity==cluster]),
                                       names(Clustering_identity[Clustering_identity==cluster])])
    list=c()
    for (j in Clusters[Clusters!=cluster]) {
      Cl=length(Clustering_identity[Clustering_identity==j])
      lrij=as.numeric(Correlation_Matrix[names(Clustering_identity[Clustering_identity==cluster]),
                                         names(Clustering_identity[Clustering_identity!=cluster])])
      list=c(list,(((1/(Ck))*(1/(Cl)))*sum(lrij)))
      list[list==0]=.Machine$double.xmin
    }
    stability_index=((1/(Ck^2))*(sum(krij)))/min(list)
    return(stability_index)
  }
  
  close(pb)
  snow::stopCluster(cl)
  x=unlist(x)
  stability_indices=sum(x)/length(Clusters)
  return(stability_indices)
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
  
  pb <- txtProgressBar(max = length(Clusters), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(1,length(Clusters)),.packages=c('fastICA','JADE'), .options.snow = opts,.inorder = T) %dopar% {
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
  close(pb)
  snow::stopCluster(cl)
  x=unlist(x)
  stability_indices=x
  names(stability_indices)=Clusters
  return(stability_indices)
}

#' Individual_Clustering
#' @description This function cluster the ICA results from iterations
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import doSNOW
Individual_Clustering <- function(Matrix,Group,ncluster,distance_measure=c('pearson','euclidean'),method='complete') {
  selected=sample(unique(Group),1)
  selected.matrix=Matrix[,names(Group)[Group==selected]]
  ncluster=ncluster
  reference=seq(1,ncluster)
  names(reference)=colnames(selected.matrix)
  results=reference
  for (i in unique(Group)[unique(Group)!=selected]) {
    temp=cbind(selected.matrix,Matrix[,names(Group)[Group==i]])
    if (distance_measure=='pearson') {
      temp=WGCNA::adjacency(as.matrix(temp),power = 1)
      temp=1-abs(temp)
    } else if (distance_measure=='euclidean') {
      temp=as.matrix(dist(t(temp)))
    }
    clustering=hclust(as.dist(temp),method = method)
    clustering=cutree(clustering,k=ncluster)
    print(hclust(as.dist(temp),method = method))
    testing.object=colnames(Matrix[,names(Group)[Group==i]])
    testing.results=c()
    for (j in testing.object) {
      ref=(clustering)[names(clustering)==j]
      ref=clustering[clustering==ref]
      ref=names(ref)[names(ref)!=j]
      testing.results=c(testing.results,reference[names(reference)==ref])
    }
    print(length(testing.results))
    names(testing.results)=testing.object
    results=c(results,testing.results)
  }
  return(results)
}




#' Individual_Clustering
#' @description This function cluster the ICA results from iterations
#' @import reshape2
#' @import foreach
#' @import doParallel
#' @import doSNOW
Direction_correction = function (Matrix,Group) {
  Reference=unique(Group)[1]
  Reference.Matrix=Matrix[,names(Group)[Group==Reference]]
  Results=list()
  for (i in unique(Group)[unique(Group)!=Reference]) {
    tested.matrix=Matrix[,c(names(Group)[Group==i])]
    correlation=cor(as.matrix(cbind(tested.matrix,Reference.Matrix)))#WGCNA::adjacency(,power = 1,type = 'signed')
    correlation=correlation[names(Group)[Group==Reference],c(names(Group)[Group==i])]
    for (j in names(Group)[Group==i]) {
      tested.matrix[,j]=tested.matrix[,j]*sign(correlation[,j][abs(correlation[,j])==max(abs(correlation[,j]))])
    }
    Results[[i]]=tested.matrix
  }
  Results=do.call(cbind,Results)
  Results=cbind(Reference.Matrix,Results)
  return(Results)
}

