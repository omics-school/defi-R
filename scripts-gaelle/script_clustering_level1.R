###
# Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
# This content is licensed under CC BY 4.0
###

# Import functions
source("graphFunctions.R")
par(mfrow = c(2,2))

# data reading
expMatrix = as.matrix(read.table("Mito_Genes.txt", header = T, row.names = 1))

# 4 clusters, Kmeans, Euclidean distance
res = kmeans(expMatrix, 4)
vecCluster = res$cluster
    
cluster1 = expMatrix[vecCluster == 1,]
cluster2 = expMatrix[vecCluster == 2,]
cluster3 = expMatrix[vecCluster == 3,]
cluster4 = expMatrix[vecCluster == 4,]

plotGenes(cluster1, yMax = ceiling(max(cluster1)), title = "Cluster1")
plotGenes(cluster2, yMax = ceiling(max(cluster2)), title = "Cluster2")
plotGenes(cluster3, yMax = ceiling(max(cluster3)), title = "Cluster3")
plotGenes(cluster4, yMax = ceiling(max(cluster4)), title = "Cluster4")

# 4 clusters, HCL, Euclidean distance
matDist = dist(expMatrix)
res = hclust(matDist)
vecCluster = cutree(res, 4)

cluster1 = expMatrix[vecCluster == 1,]
cluster2 = expMatrix[vecCluster == 2,]
cluster3 = expMatrix[vecCluster == 3,]
cluster4 = expMatrix[vecCluster == 4,]

plotGenes(cluster1, yMax = ceiling(max(cluster1)), title = "Cluster1")
plotGenes(cluster2, yMax = ceiling(max(cluster2)), title = "Cluster2")
plotGenes(cluster3, yMax = ceiling(max(cluster3)), title = "Cluster3")
plotGenes(cluster4, yMax = ceiling(max(cluster4)), title = "Cluster4")

# 4 clusters, Kmeans, Correlation distance
matDist = as.dist(1 - cor(t(expMatrix)))
res = kmeans(matDist, 4)
vecCluster = res$cluster

cluster1 = expMatrix[vecCluster == 1,]
cluster2 = expMatrix[vecCluster == 2,]
cluster3 = expMatrix[vecCluster == 3,]
cluster4 = expMatrix[vecCluster == 4,]

plotGenes(cluster1, yMax = ceiling(max(cluster1)), title = "Cluster1")
plotGenes(cluster2, yMax = ceiling(max(cluster2)), title = "Cluster2")
plotGenes(cluster3, yMax = ceiling(max(cluster3)), title = "Cluster3")
plotGenes(cluster4, yMax = ceiling(max(cluster4)), title = "Cluster4")

# 4 clusters, HCL, Correlation distance
matDist = as.dist(1 - cor(t(expMatrix)))
res = hclust(matDist)
vecCluster = cutree(res, 4)

cluster1 = expMatrix[vecCluster == 1,]
cluster2 = expMatrix[vecCluster == 2,]
cluster3 = expMatrix[vecCluster == 3,]
cluster4 = expMatrix[vecCluster == 4,]

plotGenes(cluster1, yMax = ceiling(max(cluster1)), title = "Cluster1")
plotGenes(cluster2, yMax = ceiling(max(cluster2)), title = "Cluster2")
plotGenes(cluster3, yMax = ceiling(max(cluster3)), title = "Cluster3")
plotGenes(cluster4, yMax = ceiling(max(cluster4)), title = "Cluster4")
