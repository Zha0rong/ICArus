#' ICARus
#' @description 
#' This function allows parallel running of multiple Independent component analysis (ICA) at the same time,
#' which can be used to evaluate the robustness and stability of the independent components being extracted.
#' 
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @param numberofcomponents Number of independent component to compute, must be an integer larger than 1.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param numberofcores Number of cores to use. The default is 2.
#' @param clustering_algorithm The hierarchical clustering algorithm for clustering the independent component analysis. Default is 'complete'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @importFrom GDAtools medoids
#' @import coop
#' @importFrom matrixStats rowMeans2
#' @import WGCNA
#' @import Rfast
#' @import kneedle
#' @import parallelDist
#' @export
ICARus <- function(Matrix,numberofcomponents,iteration=100,numberofcores=2,distance_measure=c('pearson','euclidean'),clustering_algorithm="complete",...) {

  WGCNA::enableWGCNAThreads(nThreads = numberofcores)

  faster_whiten=faster_ICA_whitening(Matrix)

  
  ICAResults=ParaICA(Matrix,numberofcomponents = numberofcomponents,iteration = iteration,numberofcores = numberofcores,...)

  Signature.Matrix=as.matrix(ICAResults$Signature.Matrix)
  Affiliation.Matrix=as.matrix(ICAResults$Affiliation.Matrix)
  if (distance_measure=='pearson') {
    gene_variance=matrixStats::rowVars(Signature.Matrix,useNames = T)
    gene_variance=gene_variance[order(gene_variance,decreasing = T)]
    selected_genes=names(gene_variance)#[seq(1,kneedle::kneedle(seq(1,length(gene_variance)),gene_variance)[1])]
    correlation=WGCNA::adjacency(as.matrix(Signature.Matrix[selected_genes,]),power = 1)
    Disimilarity.fixed=1-abs(correlation)
    
    Disimmilarity.Results=list()
    
    Group=stringr::str_split_fixed(colnames(Signature.Matrix),pattern = '_',n=2)[,1]
    names(Group)=colnames(Signature.Matrix)
    Disimmilarity.Results$Clustering.results.item$clustering=Individual_Clustering(Matrix=Signature.Matrix,Group=Group,ncluster=numberofcomponents,method=clustering_algorithm)
    
    
    
  } else if (distance_measure=='euclidean') {
    Group=stringr::str_split_fixed(colnames(Signature.Matrix),pattern = '_',n=2)[,1]
    names(Group)=colnames(Signature.Matrix)
    Corrected.Signature.Matrix=Direction_correction(Signature.Matrix,Group)
    rownames(Corrected.Signature.Matrix)=rownames(Signature.Matrix)
    colnames(Corrected.Signature.Matrix)=colnames(Signature.Matrix)
    print(colnames(Corrected.Signature.Matrix))
    gene_variance=matrixStats::rowVars(Corrected.Signature.Matrix,useNames = T)
    gene_variance=gene_variance[order(gene_variance,decreasing = T)]
    selected_genes=names(gene_variance)[seq(1,kneedle::kneedle(seq(1,length(gene_variance)),gene_variance)[1])]
    distance=parallelDist::parallelDist(t(Corrected.Signature.Matrix[selected_genes,]))
    Disimilarity.fixed=distance
    
    Disimmilarity.Results=list()

    Disimmilarity.Results$Clustering.results.item$clustering=Individual_Clustering(Matrix=Corrected.Signature.Matrix,Group=Group,ncluster=numberofcomponents,method=clustering_algorithm)
    
  }


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
              clustering=Disimmilarity.Results$Clustering.results.item$clustering,
              Disimilarity.fixed=Disimilarity.fixed
              ))
}

#' ICARus_est
#' @description 
#' This function can be used to evaluate different parameters for ICA. In the parameter "parameter_set" user can input multiple integers and compare which parameter is a better one.
#' 
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @param parameter_set A vector of integers. Every one of the integer in the vector needs to be larger than 2.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param numberofcores Number of cores to use. The default is 2.
#' @param clustering_algorithm The hierarchical clustering algorithm for clustering the independent component analysis. Default is 'complete'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @import coop
#' @importFrom matrixStats rowMeans2
#' @import WGCNA
#' @import Rfast
#' @import kneedle
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
    
    gene_variance=matrixStats::rowVars(Signature.Matrix,useNames = T)
    gene_variance=gene_variance[order(gene_variance,decreasing = T)]
    selected_genes=names(gene_variance)[seq(1,kneedle::kneedle(seq(1,length(gene_variance)),gene_variance)[1])]
    
    
    correlation=WGCNA::adjacency(as.matrix(Signature.Matrix[selected_genes,]),power = 1)
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


#' PCA.Estimation
#' @description 
#' This function performs PCA analysis on a given matrix, and estimate the best number of independent components by finding the elbow point in the variance explained elbow plot.
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @return A list of one vector and one ggplot,
#' @import ggplot2
#' @import ggrepel
#' @import kneedle
#' @export
PCA.Estimation <- function(Matrix=NULL) {
  Results=list()
  PCA=prcomp(t(Matrix),center=F,scale.=F)
  PCA.summary=summary(PCA)
  PCA.candidates=PCA.summary$importance[1,][order(PCA.summary$importance[1,],decreasing = T)]
  Results$ElbowPoint=kneedle(seq(1,length(PCA.candidates)),y = PCA.candidates)[1]
  
  plot_data=data.frame(index=seq(1,length(PCA.candidates)),stdev=PCA.candidates)
  plot_data$label=ifelse(plot_data$index==Results$ElbowPoint,yes=paste0('Elbow Point: ',Results$ElbowPoint),no='')
  Results$plot=ggplot(plot_data,aes(x=index,y=stdev,label=label))+geom_point(data = plot_data[plot_data$label == "",])+geom_text_repel(point.size =5)+geom_point(data = plot_data[plot_data$label != "",])+geom_vline(xintercept = Results$ElbowPoint,linewidth=1.0,colour='red',linetype="dotted")
  return(Results)
}




#' Signature_Hierarchical_Clustering
#' @description 
#' This function estimate the optimal number of clusters in a group of independent component analysis results.
#' @param Disimmilarity A dissimilarity matrix.
#' @param Affiliation.Matrix The Affiliation Matrix output by ICARus.
#' @param Signature.Matrix The Signature Matrix output by ICARus.
#' @param min_cluster Minimum number of clusters to be tested, default is 2.
#' @param max_cluster Maximum number of clusters to be tested, default is half of number of columns in  Signature Matrix.
#' @param numberofcores Number of cores to use. The default is 2.
#' @return A list of one vector and one ggplot,
#' @importFrom cluster silhouette
#' @import doParallel
#' @import ggplot2
#' @import doParallel
#' @import doSNOW
#' @import foreach
#' @import kneedle
#' @importFrom GDAtools medoids
#' @export
Signature_Hierarchical_Clustering <- function(Disimmilarity,Affiliation.Matrix,Signature.Matrix,min_cluster=2,max_cluster=ncol(Disimmilarity)/2,numberofcores=2,...) {
  clustering_results=hclust(as.dist(Disimmilarity),...)
  max_cluster=ncol(Disimmilarity)/2
  
  cl <- makeCluster(numberofcores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(seq(min_cluster,max_cluster)), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(min_cluster,max_cluster),.packages=c('cluster'), .options.snow = opts) %dopar% {
    Clustering.results=cutree(clustering_results,k = i)
    sil=silhouette(Clustering.results,dist = Disimmilarity)
    sil=mean(sil[,3])
    return(sil)
  }
  close(pb)
  stopCluster(cl)
  
  
  sils=unlist(x)
  
  silhouettes=data.frame(silhouettes=(sils),resolution=seq(min_cluster,max_cluster))
  elbowpoint=kneedle::kneedle(x=silhouettes$resolution,y=silhouettes$silhouettes)[1]
  print(kneedle::kneedle(x=silhouettes$resolution,y=silhouettes$silhouettes))
  
  #figure=ggplot(silhouettes,aes(x=resolution,y=silhouettes))+geom_point()+geom_vline(xintercept = silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouettes)])+ggtitle('Averaged Silhouettes Graph')
  figure=ggplot(silhouettes,aes(x=resolution,y=silhouettes))+geom_point()+geom_vline(xintercept = elbowpoint)+ggtitle('Averaged Silhouettes Graph')
  
  #Clustering.results=cutree(clustering_results,k = min(silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouette)]))
  Clustering.results=cutree(clustering_results,k = elbowpoint)
  
  clustering.results=Clustering.results
  clustering.results=data.frame(clustering.results)
  
  
  
  Medoids=GDAtools::medoids(as.dist(Disimmilarity),
                            Clustering.results)
  
  Medoids=names(Clustering.results)[Medoids]
  print(Medoids)
  Clustered.Signature.matrix=Signature.Matrix[,Medoids]
  Clustered.Affiliation.matrix=Affiliation.Matrix[,Medoids]
  
  colnames(Clustered.Affiliation.matrix)=seq(1,ncol(Clustered.Affiliation.matrix))
  colnames(Clustered.Signature.matrix)=seq(1,ncol(Clustered.Signature.matrix))
  colnames(Clustered.Signature.matrix)=paste('signature.',colnames(Clustered.Signature.matrix),sep = '')
  colnames(Clustered.Affiliation.matrix)=paste('signature.',colnames(Clustered.Affiliation.matrix),sep = '')
  
  
  
  
  
  
  
  
  
  
  return(list(
    'silhouettes'=silhouettes,
    'silhouettes_figure'=figure,
    'Clustered.Signature.matrix'=Clustered.Signature.matrix,
    'Clustered.Affiliation.matrix'=Clustered.Affiliation.matrix,
    'Clustering.results'=Clustering.results))
}

