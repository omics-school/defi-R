###
# Gaëlle Lelandais <gaelle.lelandais@universite-paris-saclay.fr>
# This content is licensed under CC BY 4.0
###

# Gene profiles
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

# heatmap
heatmapGenes <- function(expData, title = "my graph"){
  
  heatmap(as.matrix(expData[nrow(expData):1,]), Rowv = NA, Colv = NA, 
          margin = c(5,5), cexRow = 0.8, cexCol = 0.8, main = title)
  
  # End of function heatmapGenes()  
}
