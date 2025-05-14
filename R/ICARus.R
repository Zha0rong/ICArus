#' ICARus
#' @description 
#' This function allows parallel running of multiple Independent component analysis (ICA) at the same time,
#' which can be used to evaluate the robustness and stability of the independent components being extracted.
#' 
#' @param Matrix  A Matrix where rows are features and columns are observations.
#' @param numberofcomponents Number of independent component to compute, must be an integer larger than 1.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param numberofcores Number of cores to use. The default is 2.
#' @param clustering_algorithm Choose which clustering algorithm to use. Currently the pipeline supports Hierarchical clustering ('Hierarchical') and a clustering method based on Gale Shapely Algorithm ('MatchMaking'). The default is 'Hierarchical'.
#' @param distance_measure Choose which distance measurement to use for signatures Currently the pipeline supports pearson correlation coefficient ('pearson') and euclidean distance ('euclidean'). The default is 'pearson'.
#' @param Hierarchical.clustering.method Choose which hierarchical clustering algorithm to use. The default is 'ward.D2'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @importFrom GDAtools medoids
#' @importFrom matrixStats rowMeans2
#' @import Rfast
#' @import fastICA
#' @export
ICARus <- function(Matrix,numberofcomponents,iteration=100,numberofcores=2,distance_measure=c('pearson','euclidean'),
  clustering_algorithm=c('Hierarchical'),Hierarchical.clustering.method=c('ward.D2','ward.D','single',"single", "complete", "average","mcquitty","median","centroid"),
  scale=T,
  ...) {
  distance_measure=match.arg(distance_measure)
  clustering_algorithm=match.arg(clustering_algorithm)
  Hierarchical.clustering.method=match.arg(Hierarchical.clustering.method)
  ICAResults=ParaICA(Matrix,numberofcomponents = numberofcomponents,iteration = iteration,numberofcores = numberofcores,scale=scale,...)

  Signature.Matrix=as.matrix(ICAResults$Signature.Matrix)
  Affiliation.Matrix=as.matrix(ICAResults$Affiliation.Matrix)
  Group=stringr::str_split_fixed(colnames(Signature.Matrix),pattern = '_',n=2)[,1]
  names(Group)=colnames(Signature.Matrix)
  Corrected=Direction_correction(Signature.Matrix,Affiliation.Matrix,Group)
  Signature.Matrix=Corrected$Results.S
  Affiliation.Matrix=Corrected$Results.A
  Disimmilarity.Results=list()
  
  if (distance_measure=='pearson') {
    correlation=Rfast::cora(Signature.Matrix)
    Disimilarity.fixed=1-abs(correlation)
    cluster=hclust(d = as.dist(Disimilarity.fixed),method = Hierarchical.clustering.method)
    cluster=cutree(cluster,numberofcomponents)
    Disimmilarity.Results$Clustering.results.item$clustering=cluster
    
    } else if (distance_measure=='euclidean') {
      correlation=as.matrix(Rfast::Dist(t(Signature.Matrix),method = 'euclidean'))
      colnames(correlation)=colnames(Signature.Matrix)
      rownames(correlation)=colnames(Signature.Matrix)
      
      Disimilarity.fixed=correlation
      cluster=hclust(as.dist(correlation),method = Hierarchical.clustering.method)
      cluster=cutree(cluster,numberofcomponents)
      Disimmilarity.Results$Clustering.results.item$clustering=cluster
      correlation=Rfast::cora(Signature.Matrix)
      
  }


  Medoids=GDAtools::medoids(as.dist(Disimilarity.fixed), Disimmilarity.Results$Clustering.results.item$clustering)
  Medoids=names(Disimmilarity.Results$Clustering.results.item$clustering)[Medoids]
  Clustered.Signature.matrix=Signature.Matrix[,Medoids]
  Clustered.Affiliation.matrix=Affiliation.Matrix[,Medoids]
  colnames(Clustered.Affiliation.matrix)=seq(1,ncol(Clustered.Affiliation.matrix))
  colnames(Clustered.Signature.matrix)=seq(1,ncol(Clustered.Signature.matrix))
  colnames(Clustered.Signature.matrix)=paste('signature.',colnames(Clustered.Signature.matrix),sep = '')
  colnames(Clustered.Affiliation.matrix)=paste('signature.',colnames(Clustered.Affiliation.matrix),sep = '')
  
  Cluster_Quality=Cluster_Stability_Calculation(correlation,Clustering_identity = Disimmilarity.Results$Clustering.results.item$clustering,numberofcores = numberofcores)
  Cluster_Quality=data.frame(ICs=rep(numberofcomponents,length(Cluster_Quality)),
                             ClusterNumber=names(Cluster_Quality),
                             QualityIndex=Cluster_Quality)
  Cluster_Size=data.frame(table(Disimmilarity.Results$Clustering.results.item$clustering))
  rownames(Cluster_Size)=Cluster_Size$Var1
  colnames(Cluster_Size)=c('ClusterNumber','SignatureNumber')
  Cluster_Quality=merge(Cluster_Quality,Cluster_Size,by='ClusterNumber')
  Cluster.Quality=Cluster_Quality
  rm(Cluster_Quality,Cluster_Size)
  
  return(list(Clustered.Signature.matrix=Clustered.Signature.matrix,
              Clustered.Affiliation.matrix=Clustered.Affiliation.matrix,
              Cluster.Quality=Cluster.Quality,
              clustering=Disimmilarity.Results$Clustering.results.item$clustering,
              Disimilarity.fixed=correlation
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
#' @param clustering_algorithm Choose which clustering algorithm to use. Currently the pipeline supports Hierarchical clustering ('Hierarchical') and a clustering method based on Gale Shapely Algorithm ('MatchMaking'). The default is 'Hierarchical'.
#' @param distance_measure Choose which distance measurement to use for signatures Currently the pipeline supports pearson correlation coefficient ('pearson') and euclidean distance ('euclidean'). The default is 'pearson'.
#' @param Hierarchical.clustering.method Choose which hierarchical clustering algorithm to use. The default is 'ward.D2'.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @importFrom matrixStats rowMeans2
#' @import Rfast
#' @import fastICA
#' @export
ICARus_est <- function(Matrix,parameter_set,iteration=100,numberofcores=2,distance_measure=c('pearson','euclidean'),clustering_algorithm=c('Hierarchical'),Hierarchical.clustering.method=c('ward.D2','ward.D','single',"single", "complete", "average","mcquitty","median","centroid"),...) {
  distance_measure=match.arg(distance_measure)
  clustering_algorithm=match.arg(clustering_algorithm)
  Hierarchical.clustering.method=match.arg(Hierarchical.clustering.method)
  ICA.Tests.list=list()
  QC.Metrics=data.frame(ICs=parameter_set)
  QC.Metrics$IR=0
  QC.Metrics$MSE=0
  for (i in parameter_set){
    print(i)
    ICAResults=ParaICA(CountMatrix = Matrix,
                       numberofcomponents = i,iteration=iteration,numberofcores = numberofcores,...)
    Signature.Matrix=ICAResults$Signature.Matrix
    Affiliation.Matrix=ICAResults$Affiliation.Matrix
    Disimmilarity.Results=list()
    Group=stringr::str_split_fixed(colnames(Signature.Matrix),pattern = '_',n=2)[,1]
    names(Group)=colnames(Signature.Matrix)
    Corrected=Direction_correction(Signature.Matrix,Affiliation.Matrix,Group)
    Signature.Matrix=Corrected$Results.S
    Affiliation.Matrix=Corrected$Results.A
    rm(Corrected)
    if (distance_measure=='pearson') {
      correlation=Rfast::cora(as.matrix(Signature.Matrix))
      Disimilarity.fixed=1-(correlation)
        cluster=hclust(d = as.dist(Disimilarity.fixed),method = Hierarchical.clustering.method)
        cluster=cutree(cluster,i)
        Disimmilarity.Results$Clustering.results.item$clustering=cluster
    }
    Cluster_Quality=Cluster_Stability_Calculation((correlation),Clustering_identity = Disimmilarity.Results$Clustering.results$clustering,numberofcores = numberofcores)
    Cluster_Quality=data.frame(ICs=rep(i,length(Cluster_Quality)),ClusterNumber=names(Cluster_Quality),QualityIndex=Cluster_Quality)
    Cluster_Size=data.frame(table(Disimmilarity.Results$Clustering.results$clustering))
    rownames(Cluster_Size)=Cluster_Size$Var1
    colnames(Cluster_Size)=c('ClusterNumber','SignatureNumber')
    Cluster_Quality=merge(Cluster_Quality,Cluster_Size,by='ClusterNumber')
    ICA.Tests.list[[i]]=Cluster_Quality
  }
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
#' @param measure How to estimate the number of independent components. Default is cumulative_proportion.
#' @param scale Whether scale by column. Default is T
#' @return A list of one vector and one ggplot,
#' @import ggplot2
#' @export
PCA.Estimation <- function(Matrix=NULL,measure=c('cumulative_proportion','standard_deviation'),scale=T) {
  measure=match.arg(measure)
  
  Normalized=Matrix
  Normalized=scale(Normalized,center = T,scale = scale)
  rownames(Matrix)=rownames(Matrix)
  colnames(Matrix)=colnames(Matrix)

  Results=list()
  PCA=prcomp(t(Normalized),center=F,scale.=F)
  PCA.summary=summary(PCA)
  PCA.candidates=smooth(PCA.summary$importance[ifelse(measure=='standard_deviation',yes = 1,no = 3),])
  PCA.candidates=PCA.candidates
  Results$ElbowPoint=ifelse(measure=='standard_deviation',
    kneedle(seq(1,length(PCA.candidates)),PCA.candidates)[1]-1,
    kneedle(seq(1,length(PCA.candidates)),1/PCA.candidates)[1])
  
  plot_data=data.frame(index=seq(1,length(PCA.candidates)),stdev=PCA.candidates)
  plot_data$label=ifelse(plot_data$index==Results$ElbowPoint,yes=paste0('Elbow Point: ',Results$ElbowPoint),no='')
  Results$plot=ggplot(plot_data,aes(x=index,y=stdev,label=label))+geom_point(data = plot_data[plot_data$label == "",])+geom_point(data = plot_data[plot_data$label != "",])+geom_vline(xintercept = Results$ElbowPoint,linewidth=1.0,colour='red',linetype="dotted")
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
#' @importFrom GDAtools medoids
#' @import fastICA
#' @import pheatmap
#' @export
Signature_Hierarchical_Clustering <- function(Disimmilarity,Affiliation.Matrix,Signature.Matrix,min_cluster=2,max_cluster=ncol(Disimmilarity)/2,numberofcores=2,...) {
  clustering_results=hclust(as.dist(Disimmilarity),...)
  max_cluster=ncol(Disimmilarity)/2
  
  cl <- makeCluster(numberofcores)
  registerDoSNOW(cl)

  #opts <- list(progress = progress)
  x=foreach(i=seq(min_cluster,max_cluster),.packages=c('cluster')
            #, .options.snow = opts
            ) %dopar% {
    Clustering.results=cutree(clustering_results,k = i)
    sil=silhouette(Clustering.results,dist = Disimmilarity)
    sil=mean(sil[,3])
    return(sil)
  }
  stopCluster(cl)
  
  
  sils=unlist(x)
  
  silhouettes=data.frame(silhouettes=(sils),resolution=seq(min_cluster,max_cluster))
  elbowpoint=silhouettes$resolution[which(silhouettes$silhouettes==max(silhouettes$silhouettes))]

  #figure=ggplot(silhouettes,aes(x=resolution,y=silhouettes))+geom_point()+geom_vline(xintercept = silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouettes)])+ggtitle('Averaged Silhouettes Graph')
  figure=ggplot(silhouettes,aes(x=resolution,y=silhouettes))+geom_point()+geom_vline(xintercept = elbowpoint)+ggtitle('Averaged Silhouettes Graph')
  
  #Clustering.results=cutree(clustering_results,k = min(silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouette)]))
  Clustering.results=cutree(clustering_results,k = elbowpoint)
  
  clustering.results=Clustering.results
  clustering.results=data.frame(clustering.results)
  
  
  
  Medoids=GDAtools::medoids(as.dist(Disimmilarity),
                            Clustering.results)
  
  Medoids=names(Clustering.results)[Medoids]
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

#' ICARus_complete
#' @description 
#' This function is the main function of the package ICARus. The function automatically estimates the number of independent components to extract and output the most stable ones.
#' @param Matrix  A Matrix where rows are features and columns are observations. The count matrix needs to be normalized and log transformed before being input.
#' @param measure How to estimate the number of independent components. Default is cumulative_proportion.
#' @param iteration The number of iterations of ICA to be run, default is 100.
#' @param numberofcores Number of cores to use. The default is 2.
#' @param clustering_algorithm Choose which clustering algorithm to use. Currently the pipeline supports Hierarchical clustering ('Hierarchical') and a clustering method based on Gale Shapely Algorithm ('MatchMaking'). The default is 'Hierarchical'.
#' @param distance_measure Choose which distance measurement to use for signatures Currently the pipeline supports pearson correlation coefficient ('pearson') and euclidean distance ('euclidean'). The default is 'pearson'.
#' @param Hierarchical.clustering.method Choose which hierarchical clustering algorithm to use. The default is 'ward.D2'.
#' @param tolerance The tolerance input for fastICA to converge, the default is 1e-10.
#' @param max.iteration The maximum number of iteration of fastICA. Default is 10000.
#' @param upperbound For each parameter x the fastICA is performed y times specified by user. Then the clustering algorithm will cluster the resultsinto x clusters. Ideally each cluster will have the size of y, but that is not always the case. This parameter sets the upper bound for the size of the cluster for the signature to be considered 'stable'.
#' @param lowerbound For each parameter x the fastICA is performed y times specified by user. Then the clustering algorithm will cluster the resultsinto x clusters. Ideally each cluster will have the size of y, but that is not always the case. This parameter sets the lower bound for the size of the cluster for the signature to be considered 'stable'.
#' @param quality.index.threshold For each parameter, the signatures obtained by fastICA are evaluated using the index proposed by icasso. The index ranges from 0 to 1, where 0 means unstable at all and 1 means perfectly stable. The threshold specify the lower bound of index for each signature.
#' @return Three Matrix: 1. Stability of independent components. 2. The "A" matrix from ICA. 3. The "S" matrix from ICA.
#' @importFrom GDAtools medoids
#' @importFrom matrixStats rowMeans2
#' @import Rfast
#' @import pheatmap
#' @import fastICA
#' @export
ICARus_complete <- function(Matrix,measure=c('cumulative_proportion','standard_deviation'),iteration=100,numberofcores=4,
                            numbers_of_parameter_for_reproducibility_test=10,
                            distance_measure=c('pearson','euclidean'),
                            clustering_algorithm=c('Hierarchical'),
                            Hierarchical.clustering.method=c('ward.D2','ward.D','single',"single", "complete", "average","mcquitty","median","centroid"),
                            scale=T,
                            tolerance=1e-6,
                            max.iteration=1000,
                            upperbound=100,
                            quality.index.threshold=0.75,
                            lowerbound=50) {
  distance_measure=match.arg(distance_measure)
  clustering_algorithm=match.arg(clustering_algorithm)
  Hierarchical.clustering.method=match.arg(Hierarchical.clustering.method)
  measure=match.arg(measure)
  
  Overall.Results=list()

  Normalized=Matrix
  Estimation=PCA.Estimation(Matrix = Normalized,measure=measure,scale=scale)
  Overall.Results[["PCA_Elbow_Plot"]]=Estimation$plot
  optimal=Estimation$ElbowPoint
  
  
  Results=list()
  pb = txtProgressBar(min = 1, max = length(seq(optimal,(optimal+numbers_of_parameter_for_reproducibility_test-1))),style = 3) 
  int=1
  for (i in seq(optimal,(optimal+numbers_of_parameter_for_reproducibility_test-1))) {
    ICAResults=ICARus(Matrix = Matrix,numberofcomponents = i,
      iteration = iteration,scale=scale,numberofcores = numberofcores,clustering_algorithm = clustering_algorithm,
      distance_measure = distance_measure,Hierarchical.clustering.method = Hierarchical.clustering.method,
      tol=tolerance,maxit=max.iteration)
    Results[[paste0('IC.',i)]]=ICAResults
    rm(ICAResults)
    int=int+1
    setTxtProgressBar(pb,int)
  }
  close(pb)
  
  Overall.Results[['Raw.Results']]=Results
  Filtered.Results=Results
  
  
  for (run in names(Filtered.Results)) {
    Run=Filtered.Results[[run]]
    quality=Run$Cluster.Quality
    num.passed.checks=sum(quality$SignatureNumber==iteration&quality$QualityIndex>=quality.index.threshold)
    if (num.passed.checks==0) {
      Filtered.Results[[run]]=NULL
    } else {
      Run$Clustered.Signature.matrix=as.matrix(Run$Clustered.Signature.matrix[,
                                                                                paste0('signature.',quality$ClusterNumber[quality$SignatureNumber==iteration&quality$QualityIndex>quality.index.threshold])])
        colnames(Run$Clustered.Signature.matrix)=paste0('signature.',quality$ClusterNumber[quality$SignatureNumber==iteration&quality$QualityIndex>quality.index.threshold])
        Run$Clustered.Affiliation.matrix=as.matrix(Run$Clustered.Affiliation.matrix[,
                                                                                paste0('signature.',quality$ClusterNumber[quality$SignatureNumber==iteration&quality$QualityIndex>quality.index.threshold])])
        colnames(Run$Clustered.Affiliation.matrix)=paste0('signature.',quality$ClusterNumber[quality$SignatureNumber==iteration&quality$QualityIndex>quality.index.threshold])
        Run$Disimilarity.fixed=NULL
        Filtered.Results[[run]]=Run
    }
  }
  
  if (length(names(Filtered.Results))<(numbers_of_parameter_for_reproducibility_test/2)) {
    stop('There is not enough signatures pass quality index threshold. Please increase the max.iteration (10000 is recommended) and decrease the tolerance (1e-6 at least is recommended) for fastICA.')
  }
    
  Clustered.Signature.matrix=list()
  Clustered.Affiliation.matrix=list()
  if (length(names(Filtered.Results)))
  for (run in names(Filtered.Results)) {
    Run=Filtered.Results[[run]]
    #CS is the Clustered.Signature.matrix in this run
    CS=Run$Clustered.Signature.matrix
    colnames(CS)=paste0(run,'.',colnames(CS))
    Clustered.Signature.matrix[[run]]=CS
    #CA is the Clustered.Affiliation.matrix in this run
    
    CA=Run$Clustered.Signature.matrix
    colnames(CA)=paste0(run,'.',colnames(CA))
    Clustered.Affiliation.matrix[[run]]=CA
  }
  Clustered.Signature.matrix=do.call(cbind,Clustered.Signature.matrix)
  Clustered.Affiliation.matrix=do.call(cbind,Clustered.Affiliation.matrix)
  
  
  correlation=Rfast::cora(as.matrix(Clustered.Signature.matrix))
  
  Disimilarity.fixed=1-abs(correlation)
  
  Reproducibility.clustering=Signature_Hierarchical_Clustering(Disimilarity.fixed,Affiliation.Matrix = Clustered.Affiliation.matrix,Signature.Matrix = Clustered.Signature.matrix,
                                                               min_cluster = 2,max_cluster = ifelse(nrow(Disimilarity.fixed)/2>100,yes=100,no=nrow(Disimilarity.fixed)-1),
                                                               numberofcores = numberofcores,method=Hierarchical.clustering.method)
  
  
  
  Cluster_Quality=Cluster_Stability_Calculation(correlation,Clustering_identity = Reproducibility.clustering$Clustering.results,numberofcores = numberofcores)
  Cluster_Quality=data.frame(ICs=rep(optimal,length(Cluster_Quality)),ClusterNumber=names(Cluster_Quality),QualityIndex=Cluster_Quality)
  Cluster_Size=data.frame(table(Reproducibility.clustering$Clustering.results))
  rownames(Cluster_Size)=Cluster_Size$Var1
  colnames(Cluster_Size)=c('ClusterNumber','SignatureNumber')
  Cluster_Quality=merge(Cluster_Quality,Cluster_Size,by='ClusterNumber')
  
  Overall.Results[["Quality_Index_of_Reproducible_Signatures"]]=Cluster_Quality
  
  
  
  anno=data.frame(cluster=(Reproducibility.clustering$Clustering.results))
  
  ordering_of_cluster=table(Reproducibility.clustering$Clustering.results)
  ordering_of_cluster=ordering_of_cluster[order(ordering_of_cluster,decreasing = T)]
  
  orderofcolumn=c()
  for (i in names(ordering_of_cluster)) {
    orderofcolumn=c(orderofcolumn,
                    Reproducibility.clustering$Clustering.results[Reproducibility.clustering$Clustering.results==i])
  }
  
  cols = colorRampPalette(c( 'darkred',"white"))(50)
  
  ph=pheatmap::pheatmap(Disimilarity.fixed[names(orderofcolumn),
                                 names(orderofcolumn)],
              cluster_rows = F,cluster_cols = F,show_rownames = F,
              show_colnames = F,annotation_names_col = F,annotation_names_row = F,color = cols,silent=T,annotation_row = data.frame(cluster=as.factor(Reproducibility.clustering$Clustering.results),
                                                                                                                                    row.names = names(Reproducibility.clustering$Clustering.results)),
              annotation_col = data.frame(cluster=as.factor(Reproducibility.clustering$Clustering.results),
                                          row.names = names(Reproducibility.clustering$Clustering.results)))
  
  Overall.Results[["Reproducibility_Heatmap"]]=ph
  
  usable.cluster=Reproducibility.clustering$Clustering.results[Reproducibility.clustering$Clustering.results%in%names(table(Reproducibility.clustering$Clustering.results)
                                                                                                                      [table(Reproducibility.clustering$Clustering.results)>=(numbers_of_parameter_for_reproducibility_test/2)])]
  usable.cluster=names(usable.cluster)
  usable.cluster=usable.cluster[grepl(paste0('IC.',optimal),usable.cluster)]
  usable.cluster=gsub(paste0('IC.',optimal,'.'),'',usable.cluster)
  Ideal_Results=Results[[paste0('IC.',optimal)]]
  
  Consensus.Affiliation.matrix=Ideal_Results$Clustered.Affiliation.matrix[,usable.cluster]
  Consensus.Signature.matrix=Ideal_Results$Clustered.Signature.matrix[,usable.cluster]
  Overall.Results[["Reproducible_Signature_Matrix"]]=Consensus.Signature.matrix
  Overall.Results[["Reproducible_Affiliation_Matrix"]]=Consensus.Affiliation.matrix
  
  return(Overall.Results)
}


