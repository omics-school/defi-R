###
# Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
# This content is licensed under CC BY 4.0
###

# This function is useful to draw gene expression profiles
plotGenes <- function(expData, title = "", yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(0, ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}

#######################################
# Clustering project - Part 1
#######################################

# 1) Read the data file
expMatrix = read.table("Mito_Genes.txt", header = T, row.names = 1)

# 2) Create 4 clusters, using Kmeans algorithm and default distance
# method (Euclidean)
# -- use the function kmeans()
res = kmeans(expMatrix, 4)
# --> "res" is a complex object, it is a list in which several other objects
# are associated. The $ is used to access only one object from the list
# (here the "cluster" information).
vecCluster = res$cluster
# --> "vecCluster" is a vector, in which the values are the cluster name
# (here from 1 to 4). The values are labeled with the names of genes.

# 3) Find all genes that belong to each cluster, extract their associated
# gene expression profiles and draw them in a graph (using the function
# "PlotGenes()", see below).

# Cluster 1
# -- search for the genes which belong to cluster 1
geneCluster1 = names(which(vecCluster == 1))
# -- extract genes expression profiles for these genes
cluster1 = expMatrix[geneCluster1,]
# -- draw gene expression profiles using the function "PlotGenes()"
plotGenes(cluster1, title = "Cluster 1", yMax = max(expMatrix))

########### ----------------------------------------------------------
# TO DO : Repeat the process for cluster 2, 3 and 4.
# Create a PDF file with the obtained four graphics.
########### ----------------------------------------------------------

# Write the code here

# 4) Create 4 clusters, using HCL algorithm and Euclidean distance
# method
# -- Distance matrix is calculated with the dist() function
matDist = dist(expMatrix)
# -- HCL algorithm is applied, starting from the distance matrix
res = hclust(matDist)
# -- Clusters are created, in a second step
vecCluster = cutree(res, 4)

########### ----------------------------------------------------------
# TO DO : As you did with the results obtained with "kmeans()"
# draw the expression profiles for genes in clusters 1, 2, 3 and 4.
# Create a PDF file with the obtained four graphics.
########### ----------------------------------------------------------

# Write the code here

# 5) Create 4 clusters, Kmeans, Correlation distance
# -- Distance matrix is calculated with the dist() and cor() functions
matDist = as.dist(1 - cor(t(expMatrix)))
# -- Kmeans function can take a distance matrix as first argument
res = kmeans(matDist, 4)
vecCluster = res$cluster

########### ----------------------------------------------------------
# TO DO : As you did with the results obtained with euclidean distance,
# draw the expression profiles for genes in clusters 1, 2, 3 and 4.
# Create a PDF file with the obtained four graphics.
########### ----------------------------------------------------------

# Write the code here

# 6) Create 4 clusters, HCL, Correlation distance

########### ----------------------------------------------------------
# TO DO : As you did with the results obtained with kmeans(),
# draw the expression profiles for genes in clusters 1, 2, 3 and 4.
# Create a PDF file with the obtained four graphics.
########### ----------------------------------------------------------

# Write the code here

# This is the end of the first part.