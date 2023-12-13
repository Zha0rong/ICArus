#' ICARus
#' @description 
#' This function allows parallel running of multiple Independent component analysis (ICA) at the same time,
#' which can be used to evaluate the robustness and stability of the independent components being extracted.
#' 
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @param numberofcomponents Number of independent component to compute, must be an integer larger than 1.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param numberofcores Number of threads to use. The default is 2.
#' @param clustering_algorithm The hierarchical clustering algorithm for clustering the independent component analysis. Default is 'complete'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @import GDAtools
#' @import coop
#' @importFrom matrixStats rowMeans2
#' @import WGCNA
#' @import Rfast
#' @export
ICARus <- function(Matrix,numberofcomponents,iteration=100,numberofcores=2,clustering_algorithm="complete",...) {

  WGCNA::enableWGCNAThreads(nThreads = numberofcores)

  faster_whiten=faster_ICA_whitening(Matrix)

  
  ICAResults=ParaICA(Matrix,numberofcomponents = numberofcomponents,iteration = iteration,numberofcores = numberofcores,...)
  print('Finished ParaICA')

  Signature.Matrix=as.matrix(ICAResults$Signature.Matrix)
  Affiliation.Matrix=as.matrix(ICAResults$Affiliation.Matrix)
  correlation=WGCNA::adjacency(as.matrix(Signature.Matrix),power = 1)
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
  return(list(Clustered.Signature.matrix=Clustered.Signature.matrix,
              Clustered.Affiliation.matrix=Clustered.Affiliation.matrix,
              Cluster.Quality=Cluster.Quality,
              clustering=Disimmilarity.Results$Clustering.results.item$clustering
              ))
}

#' ICARus_est
#' @description 
#' This function can be used to evaluate different parameters for ICA. In the parameter "parameter_set" user can input multiple integers and compare which parameter is a better one.
#' 
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @param parameter_set A vector of integers. Every one of the integer in the vector needs to be larger than 2.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param numberofcores Number of threads to use. The default is 2.
#' @param clustering_algorithm The hierarchical clustering algorithm for clustering the independent component analysis. Default is 'complete'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @import GDAtools
#' @import coop
#' @importFrom matrixStats rowMeans2
#' @import WGCNA
#' @import Rfast
#' @export
ICARus_est <- function(Matrix,parameter_set,iteration=100,numberofcores=2,clustering_algorithm='complete',...) {
  ICA.Tests.list=list()
  enableWGCNAThreads(nThreads = numberofcores)
  
  
  QC.Metrics=data.frame(ICs=parameter_set)
  QC.Metrics$IR=0
  QC.Metrics$MSE=0
  faster_whiten=faster_ICA_whitening(Matrix)
  
  for (i in parameter_set){
    print(i)
    temp=faster_whiten
    temp$K <- matrix(temp$K[1:i, ], i, temp$p)
    temp$X1 <- mat.mult(temp$K, temp$X)
    ICAResults=ParaICA(CountMatrix = Matrix,faster_whiten =  temp,numberofcomponents = i,iteration=iteration,numberofcores = numberofcores,...)
    Signature.Matrix=ICAResults$Signature.Matrix
    Affiliation.Matrix=ICAResults$Affiliation.Matrix
    correlation=WGCNA::adjacency(as.matrix(Signature.Matrix),power = 1)
    Disimilarity.fixed=1-abs(correlation)
    Disimmilarity.Results=list()
    Group=stringr::str_split_fixed(colnames(Signature.Matrix),pattern = '_',n=2)[,1]
    names(Group)=colnames(Signature.Matrix)
    Disimmilarity.Results$Clustering.results.item$clustering=Individual_Clustering(Matrix=Signature.Matrix,Group=Group,ncluster=i,method=clustering_algorithm)
    a=Cluster_Stability_Calculation(abs(correlation),Clustering_identity = Disimmilarity.Results$Clustering.results$clustering,numberofcores = numberofcores)
    a=data.frame(ICs=rep(i,length(a)),ClusterNumber=names(a),QualityIndex=a)
    b=data.frame(table(Disimmilarity.Results$Clustering.results$clustering))
    rownames(b)=b$Var1
    colnames(b)=c('ClusterNumber','SignatureNumber')
    a=merge(a,b,by='ClusterNumber')
    ICA.Tests.list[[i]]=a
  }
  disableWGCNAThreads()
  ICA.Tests=ICA.Tests.list
  ICA.Tests=do.call(rbind,ICA.Tests)
  ICA.Tests.median=aggregate(QualityIndex~ICs,data=ICA.Tests,median)
  ICA.Tests.mean=aggregate(QualityIndex~ICs,data=ICA.Tests,mean)
  return(list(ICA.Tests.results=ICA.Tests,
              ICA.Tests.median=ICA.Tests.median,
              ICA.Tests.mean=ICA.Tests.mean,
              QC.Metrics=QC.Metrics))
  
  
}
