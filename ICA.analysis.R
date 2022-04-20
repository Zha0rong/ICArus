library(BUSpaRse)
library(tximport)
library(EBSeq)
library(fastICA)
library(JADE)
library(cluster)
library(factoextra)
library(parallelDist)
library(coop)
HPV16.Study=metadata[metadata$infection%in%c('HPV16 Infected','HPV16 Partial Existence','Uninfected'),]
files=paste0('Counts//',HPV16.Study$File.ID,'.isoforms.results')
names(files)=(HPV16.Study$Case.ID)

Counts=tximport(files,type = 'rsem',countsFromAbundance = 'lengthScaledTPM',tx2gene = ProteinCoding)


CountsNormData = Counts$abundance
CountsNormData=CountsNormData[rowSums(CountsNormData)>50,]
#size=MedianNorm(CountsNormData)
#CountsNormData=GetNormalizedMat(CountsNormData,size)
CountsNormData = log(CountsNormData+1)





Overall_stability_score=c()
for (i in 3:50){
  ICAResults=ParaICA(CountsNormData,i,100,8)
  Signature.Matrix=ICAResults$Signature.Matrix
  Affiliation.Matrix=ICAResults$Affiliation.Matrix
  Disimmilarity=pcor(Signature.Matrix)
  Correlation=Disimmilarity
  Disimmilarity=(1-(Disimmilarity))/2
  
  clusteringresults=hclust(as.dist(Disimmilarity))
  silhouettes=c()
  for (cluster in 2:200) {
    sil_cl <- silhouette(x=cutree(clusteringresults, k=cluster) ,dist=as.dist(Disimmilarity),full = T)
    silhouettes=c(silhouettes,mean(sil_cl[,3]))
    rm(cluster,sil_cl)
  }
  names(silhouettes)=as.character(seq(2,200))
  silhouettes=data.frame(silhouettes)
  silhouettes$ICs=as.numeric(rownames(silhouettes))
  clustering.results=cutree(clusteringresults, k=silhouettes$ICs[silhouettes$silhouettes==max(silhouettes$silhouettes)])

  stability_score=Stability_Score_Calculation(Correlation,Clustering_identity = clustering.results)
  Overall_stability_score=c(Overall_stability_score,stability_score)
  
  
}










