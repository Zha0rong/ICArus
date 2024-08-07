#' ParaICA
#' @description This function runs ICA using multiple threads.
#' @import foreach
#' @import doParallel
#' @import doSNOW
#' @import truncnorm
#' @import Rfast
#' @import snow

ParaICA <- function(CountMatrix,numberofcomponents,iteration,numberofcores=2,...) {

  cl <- snow::makeCluster(numberofcores)
  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(max = iteration, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach::foreach(i=seq(1,iteration),.packages=c('fastICA','Rfast'), .options.snow = opts) %dopar% {
    resICA=fastICA(CountMatrix,row.norm = F,n.comp = numberofcomponents,method = 'C',...)
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

#' Individual_Matching
#' @description This function matches signatures from different iterations using galeShapley.marriageMarket
#' @importFrom matchingR galeShapley.marriageMarket
Individual_Matching <- function(Matrix,Group,ncluster) {
  
  selected=sample(unique(Group),1)
  selected.matrix=Matrix[,names(Group)[Group==selected]]
  ncluster=ncluster
  reference=seq(1,ncluster)
  names(reference)=colnames(selected.matrix)
  results=reference
  for (i in unique(Group)[unique(Group)!=selected]) {
    match.reference=Matrix[names(Group)[Group==selected],names(Group)[Group==i]]
    match.test=t(match.reference)
    match=matchingR::galeShapley.marriageMarket(match.reference,match.test)
    match=as.integer(match$proposals)
    names(match)=colnames(match.reference)
    results=c(results,match)
  }
  return(results)
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

