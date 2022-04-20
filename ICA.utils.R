ParaICA <- function(CountMatrix,numberofcomponents,iteration,numberofcores=4,...) {
  require(foreach)
  require(doParallel)
  require(doSNOW)
  cl <- makeCluster(numberofcores)
  registerDoSNOW(cl)

  pb <- txtProgressBar(max = iteration, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(1,iteration),.packages=c('fastICA','JADE','MASS'), .options.snow = opts) %dopar% {
    resICA=fastICA(CountMatrix,n.comp = numberofcomponents,...)
    Affiliation.Matrix=(as.matrix(resICA$A))
    Signature.Matrix=as.matrix(resICA$S)
    Demixing.Matrix=t(ginv(Affiliation.Matrix))
    colnames(Affiliation.Matrix)=colnames(CountMatrix)
    rownames(Affiliation.Matrix)=paste0('n.',seq(1,nrow(Affiliation.Matrix)))
    colnames(Signature.Matrix)=paste0('n','.',seq(1,ncol((Signature.Matrix))))
    
    rownames(Demixing.Matrix)=paste0('n','.',seq(1,nrow((Demixing.Matrix))))
    colnames(Demixing.Matrix)=colnames(CountMatrix)
    
    Results=list()
    #Make sure the matrix:
    #Affiliation: samples * Signature
    #Signature: genes * Signature
    #Unmixing: samples * Signature
    Results[['Affiliation.Matrix']]=t(Affiliation.Matrix)
    Results[['Signature.Matrix']]=Signature.Matrix
    Results[['Demixing.Matrix']]=Demixing.Matrix
    
    return(Results)
  }
  close(pb)
  stopCluster(cl)
  Affiliation.Matrix=list()
  Signature.Matrix=list()
  Demixing.Matrix=list()
  for (i in 1:length(x)) {
    tmp=x[[i]][['Affiliation.Matrix']]
    colnames(tmp)=paste0('iteration.',i,'_',colnames(tmp))
    Affiliation.Matrix[[i]]=tmp
    
    tmp=x[[i]][['Signature.Matrix']]
    colnames(tmp)=paste0('iteration.',i,'_',colnames(tmp))
    Signature.Matrix[[i]]=tmp
    
    tmp=x[[i]][['Demixing.Matrix']]
    rownames(tmp)=paste0('iteration.',i,'_',rownames(tmp))
    Demixing.Matrix[[i]]=tmp
  }
  Affiliation.Matrix=do.call(cbind,Affiliation.Matrix)
  Signature.Matrix=do.call(cbind,Signature.Matrix)
  Demixing.Matrix=t(do.call(rbind,Demixing.Matrix))
  return(list('Signature.Matrix'=Signature.Matrix,
              'Affiliation.Matrix'=Affiliation.Matrix,
              'Demixing.Matrix'=Demixing.Matrix))
}

Stability_Score_Calculation <- function(Correlation_Matrix,Clustering_identity,numberofcores=6) {
  library(reshape2)
  Clusters=unique(Clustering_identity)
  melted=melt(Correlation_Matrix)
  melted=melted[melted$Var1!=melted$Var2,]
  stability_indices=c()
  
  require(foreach)
  require(doParallel)
  require(doSNOW)
  cl <- makeCluster(numberofcores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(Clusters), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(1,length(Clusters)),.packages=c('fastICA','JADE','MASS'), .options.snow = opts) %dopar% {
    cluster=Clusters[i]
    Ck=length(Clustering_identity[Clustering_identity==cluster])
    krij=melted$value[melted$Var1%in%names(Clustering_identity[Clustering_identity==cluster])&
                        melted$Var2%in%names(Clustering_identity[Clustering_identity==cluster])]
    list=c()
    for (j in Clusters[Clusters!=cluster]) {
      Cl=length(Clustering_identity[Clustering_identity==j])
      lrij=melted$value[melted$Var1%in%names(Clustering_identity[Clustering_identity==cluster])&
                          melted$Var2%in%names(Clustering_identity[Clustering_identity==j])]
      list=c(list,(((1/(Ck))*(1/(Cl)))*sum(lrij)))
      list[list==0]=.Machine$double.xmin
    }
    stability_index=((1/(Ck^2))*(sum(krij)))/min(list)
    return(stability_index)
  }
  
  close(pb)
  stopCluster(cl)
  x=unlist(x)
  stability_indices=sum(x)/length(Clusters)
  names(stability_indices)='1'
  return(stability_indices)
}

Signature_Leiden_Clustering <- function(Disimmilarity,Affiliation.Matrix,Signature.Matrix,min_res=0.1,max_res=10,numberofcores=6) {
  require('Seurat')
  require('igraph')
  require('ggplot2')
  require('cluster')
  Graph=FindNeighbors(Disimmilarity,k.param = 20,compute.SNN = T,verbose = F)
  Graph.snn=(Graph$snn)
  tmp=dim(Graph.snn)
  Graph.snn=graph_from_adjacency_matrix(Graph.snn[1:tmp[1],1:tmp[2]], mode='undirected', weighted=T)
  rm(Graph)
  
  cl <- makeCluster(numberofcores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(seq(min_res,max_res,0.1)), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(min_res,max_res,0.1),.packages=c('Seurat','igraph','cluster'), .options.snow = opts) %dopar% {
    clusterlouvain <- cluster_leiden(Graph.snn,objective_function='Modularity',resolution_parameter=i)
    Clustering.results=clusterlouvain$membership
    names(Clustering.results)=clusterlouvain$names
    sil=silhouette(Clustering.results,dist = Disimmilarity)
    sil=mean(sil[,3])
    return(sil)
  }
  close(pb)
  stopCluster(cl)
  sils=unlist(x)
  silhouettes=data.frame(silhouettes=sils,resolution=seq(min_res,max_res,0.1))
  
  
  
  
  #sils=c()
  #for (i in seq(min_res,max_res,0.1)){
  #  clusterlouvain <- cluster_leiden(Graph.snn,objective_function='Modularity',resolution_parameter=i)
  #  Clustering.results=clusterlouvain$membership
  #  names(Clustering.results)=clusterlouvain$names
  #  sil=silhouette(Clustering.results,dist = Disimmilarity)
  #  sil=mean(sil[,3])
  #  sils=c(sils,sil)
  #  rm(Clustering.results,clusterlouvain,sil)
  #}
  #silhouettes=data.frame(silhouettes=sils,resolution=seq(min_res,max_res,0.1))
  figure=ggplot(silhouettes,aes(x=resolution,y=silhouettes))+geom_point()+geom_vline(xintercept = silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouettes)])+ggtitle('Averaged Silhouettes Graph')
  clusterlouvain <- cluster_leiden(Graph.snn,objective_function='Modularity',
                                   resolution_parameter=max(silhouettes$resolution[silhouettes$silhouettes==min(silhouettes$silhouettes)]))

  Clustering.results=clusterlouvain$membership
  names(Clustering.results)=clusterlouvain$names
  clustering.results=Clustering.results
  clustering.results=data.frame(clustering.results)
  Clustered.Signature.matrix=matrix(rep(0,length(unique(clustering.results$clustering.results))*nrow(Signature.Matrix)),
                                    nrow = nrow(Signature.Matrix)
  )
  rownames(Clustered.Signature.matrix)=rownames(Signature.Matrix)
  colnames(Clustered.Signature.matrix)=unique(clustering.results$clustering.results)
  
  Clustered.Affiliation.matrix=matrix(rep(0,length(unique(clustering.results$clustering.results))*nrow(Affiliation.Matrix)),
                                      nrow = nrow(Affiliation.Matrix)
  )
  rownames(Clustered.Affiliation.matrix)=rownames(Affiliation.Matrix)
  colnames(Clustered.Affiliation.matrix)=unique(clustering.results$clustering.results)
  
  
  
  for (i in 1:ncol(Clustered.Signature.matrix)) {
    if (length(rownames(clustering.results)[clustering.results$clustering.results==i])>1){
    tmp_sig=Signature.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
    tmp_sig=apply(tmp_sig,1,median)
    Clustered.Signature.matrix[,i]=tmp_sig
    tmp_sig=Affiliation.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
    tmp_sig=apply(tmp_sig,1,median)
    Clustered.Affiliation.matrix[,i]=tmp_sig}
    else {
      tmp_sig=Signature.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
      Clustered.Signature.matrix[,i]=tmp_sig
      tmp_sig=Affiliation.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
      Clustered.Affiliation.matrix[,i]=tmp_sig
      
    }
  }
  
  
  
  
  
  return(list('Graph.snn'=Graph.snn,
              'silhouettes'=silhouettes,
              'silhouettes_figure'=figure,
              'Clustered.Signature.matrix'=Clustered.Signature.matrix,
              'Clustered.Affiliation.matrix'=Clustered.Affiliation.matrix,
              'Clustering.results'=Clustering.results))
}

Signature_Hierarchical_Clustering <- function(Disimmilarity,Affiliation.Matrix,Signature.Matrix,min_cluster=2,max_cluster=200,numberofcores=8,...) {
  clustering_results=hclust(as.dist(Disimmilarity),...)
  
  
  cl <- makeCluster(numberofcores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = length(seq(min_cluster,max_cluster)), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  x=foreach(i=seq(min_cluster,max_cluster),.packages=c('Seurat','igraph','cluster'), .options.snow = opts) %dopar% {
    Clustering.results=cutree(clustering_results,k = i)
    sil=silhouette(Clustering.results,dist = Disimmilarity)
    sil=mean(sil[,3])
    return(sil)
  }
  close(pb)
  stopCluster(cl)
  
  
  sils=unlist(x)
  silhouettes=data.frame(silhouettes=sils,resolution=seq(min_cluster,max_cluster))
  figure=ggplot(silhouettes,aes(x=resolution,y=silhouettes))+geom_point()+geom_vline(xintercept = silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouettes)])+ggtitle('Averaged Silhouettes Graph')
  
  Clustering.results=cutree(clustering_results,k = silhouettes$resolution[silhouettes$silhouettes==max(silhouettes$silhouette)])
  
  clustering.results=Clustering.results
  clustering.results=data.frame(clustering.results)
  Clustered.Signature.matrix=matrix(rep(0,length(unique(clustering.results$clustering.results))*nrow(Signature.Matrix)),
                                    nrow = nrow(Signature.Matrix)
  )
  rownames(Clustered.Signature.matrix)=rownames(Signature.Matrix)
  colnames(Clustered.Signature.matrix)=unique(clustering.results$clustering.results)
  
  Clustered.Affiliation.matrix=matrix(rep(0,length(unique(clustering.results$clustering.results))*nrow(Affiliation.Matrix)),
                                      nrow = nrow(Affiliation.Matrix)
  )
  rownames(Clustered.Affiliation.matrix)=rownames(Affiliation.Matrix)
  colnames(Clustered.Affiliation.matrix)=unique(clustering.results$clustering.results)
  
  
  
  for (i in 1:ncol(Clustered.Signature.matrix)) {
    if (length(rownames(clustering.results)[clustering.results$clustering.results==i])>1){
      tmp_sig=Signature.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
      tmp_sig=apply(tmp_sig,1,median)
      Clustered.Signature.matrix[,i]=tmp_sig
      tmp_sig=Affiliation.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
      tmp_sig=apply(tmp_sig,1,median)
      Clustered.Affiliation.matrix[,i]=tmp_sig}
    else {
      tmp_sig=Signature.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
      Clustered.Signature.matrix[,i]=tmp_sig
      tmp_sig=Affiliation.Matrix[,rownames(clustering.results)[clustering.results$clustering.results==i]]
      Clustered.Affiliation.matrix[,i]=tmp_sig
      
    }
  }
  
  
  
  
  
  return(list(
              'silhouettes'=silhouettes,
              'silhouettes_figure'=figure,
              'Clustered.Signature.matrix'=Clustered.Signature.matrix,
              'Clustered.Affiliation.matrix'=Clustered.Affiliation.matrix,
              'Clustering.results'=Clustering.results))
}

