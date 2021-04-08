###
# Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
# This content is licensed under CC BY 4.0
###

# Function 1
dataReading <- function(FileName, folder = "./"){
  
  data = read.table(paste(folder, FileName, sep = ""), header = T, row.names = 1)
  expMatrix = as.matrix(data)
  
  return(expMatrix)
  
# end of function dataReading()  
}

# Function 2
distCalc <- function(expMatrix, distName){
  
  if(distName == "EUCLIDEAN"){
    
    distMat = dist(expMatrix)
    return(distMat)
    
  }else if(distName == "CORRELATION"){

    distMat = as.dist(1 - cor(t(expMatrix)))
    return(distMat)
    
  }else{
    
    print(paste("The distance", distName, "is unknown"))
  }

# End of function distCalc()    
}

# Function 3
createClusters <- function(distMat, algo, N){
  
  if(algo == "KMEANS"){
    
    res = kmeans(distMat, N)
    vecClusters = res$cluster
    return(vecClusters)
    
  }else if(algo == "HCL"){
    
    res = hclust(distMat)
    vecClusters = cutree(res, N)
    return(vecClusters)
    
  }else{

    print(paste("The algorithm", algo, "is unknown"))
  }
  
# End of function createClusters()  
}

# Function 4
extractCluster <- function(expMat, vecCluster, number){
  
  clusterMat = expMat[which(vecCluster == number),]
  return(clusterMat)
  
}

# Function 5
drawGraph <- function(expMat, type, title){
  
  source("graphFunctions.R")
  
  if(type == "PROFILES"){
    
    plotGenes(expMat, yMax = ceiling(max(expMat)), title = title)
        
  }else if(type == "HEATMAP"){
    
    heatmapGenes(expMat, title = title)
    
  }else{

    print(paste("The representation", type, "is unknown"))
  }
  
# End of function drawClusters()  
}


# Main function
myClustering <- function(fileName, distance, algo, N, outputGraph){

  mat      = dataReading(fileName)
  matDist  = distCalc(mat, distance)
  clustRes = createClusters(matDist, algo, N)

  for(i in 1:N){
    clust = extractCluster(mat, clustRes, i)
    drawGraph(clust, outputGraph, title = paste("Cluster #", i, sep = ""))
    
  }
  
# End of function myClustering()
}


# Tests of the function
myClustering(fileName = "Mito_Genes.txt", 
             distance = "EUCLIDEAN",
             algo = "KMEANS",
             N    = 4,
             outputGraph = "PROFILES")

myClustering(fileName = "Mito_Genes.txt", 
             distance = "EUCLIDEAN",
             algo = "HCL",
             N    = 4,
             outputGraph = "PROFILES")

myClustering(fileName = "Mito_Genes.txt", 
             distance = "CORRELATION",
             algo = "KMEANS",
             N    = 4,
             outputGraph = "PROFILES")

myClustering(fileName = "Mito_Genes.txt", 
             distance = "CORRELATION",
             algo = "KMEANS",
             N    = 4,
             outputGraph = "PROFILES")


