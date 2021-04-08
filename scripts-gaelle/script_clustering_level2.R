###
# Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
# This content is licensed under CC BY 4.0
###

# Import functions
source("graphFunctions.R")

# data reading
expMatrix = as.matrix(read.table("Mito_Genes.txt", header = T, row.names = 1))

# number of clusters
N = 4
# distance (EUCLIDEAN or CORRELATION)
distance = "EUCLIDEAN"
# algorithm for clustering (KEMANS or HCL)
algo = "KMEANS"


# Kmeans, Eclidean distance or correlation distance
if(algo == "KMEANS"){
  
  if(distance == "EUCLIDEAN"){
    
    res = kmeans(expMatrix, N)
    vecCluster = res$cluster
    
    for(i in 1:N){
      cluster = expMatrix[vecCluster == i,]
      plotGenes(cluster, yMax = ceiling(max(cluster)), title = paste("Cluster #", i, sep = ""))
    }
    
  }else if(distance == "CORRELATION"){

    matDist = as.dist(1 - cor(t(expMatrix)))
    res = kmeans(matDist, N)
    vecCluster = res$cluster
    
    for(i in 1:N){
      cluster = expMatrix[vecCluster == i,]
      plotGenes(cluster, yMax = ceiling(max(cluster)), title = paste("Cluster #", i, sep = ""))
    }
  }
}


# HCL, Eclidean distance or correlation distance
if(algo == "HCL"){
  
  if(distance == "EUCLIDEAN"){
    
    matDist = dist(expMatrix)
    res = hclust(matDist)
    vecCluster = cutree(res, N)
    
    for(i in 1:N){
      cluster = expMatrix[vecCluster == i,]
      plotGenes(cluster, yMax = ceiling(max(cluster)), title = paste("Cluster #", i, sep = ""))
    }
    
  }else if(distance == "CORRELATION"){
    
    matDist = as.dist(1 - cor(t(expMatrix)))
    res = hclust(matDist)
    vecCluster = cutree(res, N)
    
    for(i in 1:N){
      cluster = expMatrix[vecCluster == i,]
      plotGenes(cluster, yMax = ceiling(max(cluster)), title = paste("Cluster #", i, sep = ""))
    }
  }
}
