#' ICARus
#' @description 
#' This function allows parallel running of multiple Independent component analysis (ICA) at the same time,
#' which can be used to evaluate the robustness and stability of the independent components being extracted.
#' 
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @param numberofcomponents Number of independent component to compute, must be an integer larger than 1.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param clustering_algorithm The hierarchical clustering algorithm for clustering the independent component analysis. Default is 'complete'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @import GDAtools
#' @import coop
#' @import matrixStats
#' @import WGCNA
#' @import Rfast
#' @export
ICARus <- function(Matrix,numberofcomponents,iteration=100,numberofcores,clustering_algorithm="complete",...) {

  WGCNA::enableWGCNAThreads(nThreads = numberofcores)

  faster_whiten=faster_ICA_whitening(Matrix)

  
  ICAResults=ParaICA(Matrix,numberofcomponents = numberofcomponents,iteration = iteration,numberofcores = numberofcores,...)
  print('Finished ParaICA')

  Signature.Matrix=as.matrix(ICAResults$Signature.Matrix)
  Affiliation.Matrix=as.matrix(ICAResults$Affiliation.Matrix)
  correlation=WGCNA::adjacency(as.matrix(Affiliation.Matrix),power = 1)
  Disimilarity.fixed=1-abs(correlation)
  Disimmilarity.Results=list()

  Group=stringr::str_split_fixed(colnames(Signature.Matrix),pattern = '_',n=2)[,1]
  names(Group)=colnames(Signature.Matrix)
  Matrix=Signature.Matrix
  Disimmilarity.Results$Clustering.results.item$clustering=Individual_Clustering(Matrix=Matrix,Group=Group,ncluster=numberofcomponents,method=clustering_algorithm)
  
  Medoids=GDAtools::medoids(as.dist(Disimilarity.fixed), Disimmilarity.Results$Clustering.results.item$clustering)
  Medoids=names(Disimmilarity.Results$Clustering.results.item$clustering)[Medoids]
  Clustered.Signature.matrix=Signature.Matrix[,Medoids]
  Clustered.Affiliation.matrix=Affiliation.Matrix[,Medoids]
  colnames(Clustered.Affiliation.matrix)=seq(1,ncol(Clustered.Affiliation.matrix))
  colnames(Clustered.Signature.matrix)=seq(1,ncol(Clustered.Signature.matrix))
  colnames(Clustered.Signature.matrix)=paste('signature.',colnames(Clustered.Signature.matrix),sep = '')
  colnames(Clustered.Affiliation.matrix)=paste('signature.',colnames(Clustered.Affiliation.matrix),sep = '')
  a=Cluster_Stability_Calculation(abs(correlation),Clustering_identity = Disimmilarity.Results$Clustering.results.item$clustering,numberofcores = 6)
  a=data.frame(ICs=rep(numberofcomponents,length(a)),ClusterNumber=names(a),QualityIndex=a)
  b=data.frame(table(Disimmilarity.Results$Clustering.results.item$clustering))
  rownames(b)=b$Var1
  colnames(b)=c('ClusterNumber','SignatureNumber')
  a=merge(a,b,by='ClusterNumber')
  Cluster.Quality=a
  rm(a,b)
  
  WGCNA::disableWGCNAThreads()
  
  for (i in 1:nrow(Cluster.Quality)) {
    if (Cluster.Quality$SignatureNumber[i]>iteration) {
      Cluster.Quality$QualityIndex[i]=Cluster.Quality$QualityIndex[i]#/(Cluster.Quality$SignatureNumber[i])
    }
    if (Cluster.Quality$SignatureNumber[i]<iteration) {
      Cluster.Quality$QualityIndex[i]=Cluster.Quality$QualityIndex[i]#(Cluster.Quality$SignatureNumber[i])
    } 
  }
  
  
  
  return(list(Clustered.Signature.matrix=Clustered.Signature.matrix,
              Clustered.Affiliation.matrix=Clustered.Affiliation.matrix,
              Cluster.Quality=Cluster.Quality,
              clustering=Disimmilarity.Results$Clustering.results.item$clustering
              ))
}


