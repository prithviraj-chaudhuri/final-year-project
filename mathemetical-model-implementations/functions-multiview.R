

#This Program is for feature selection
#This File Contans all the functions used
#Author - Prithviraj Chaudhuri


#Helper Function to print a value



#Function to read Data
readData <- function(NMIfilename, isHeader){
  data <- read.table(NMIfilename, header=isHeader)
  print("Data Read")
  flush.console()
  return(data)
}




#Function to get MI from CSV Data
#getNMI <- funtion(CSVDataMatrix){
#  NMIDataMatrix <- makemim(CSVDataMatrix)
#  return(NMIDataMatrix)
#}




#function to scale the CSV Data
scaleing <- function(CSVDataMatrix){
  
  numOfNodes <- ncol(CSVDataMatrix)
  numOfSamples <- nrow(CSVDataMatrix)
  
  minScaledVal <- 0
  maxScaledVal <- 1
  minVal <- 0
  maxVal <- 0
  scaledCsvDataMatrix <- matrix( c(0), nrow = numOfSamples, ncol = numOfNodes) 
  
  for(j in 1:numOfNodes){
    
    minVal <- sapply(CSVDataMatrix[j], min, na.rm = TRUE)
    maxVal <- sapply(CSVDataMatrix[j], max, na.rm = TRUE)
    
    for(i in 1:numOfSamples){
      scaledCsvDataMatrix[i,j] <- ((CSVDataMatrix[i,j]-minVal)/(maxVal-minVal))*(maxScaledVal-minScaledVal)+minScaledVal
    }
  }
  return(scaledCsvDataMatrix)
}






#Funtion to get mean of a vector
getMean <- function(scaledCsvDataMatrix, col){
  
  numOfSamples <- nrow(scaledCsvDataMatrix)
  
  sum <- 0
  for(i in 1:numOfSamples){
    sum <- sum + scaledCsvDataMatrix[i,col]
  }
  mean <- sum/numOfSamples
  return(mean)
}






#Function to calculate the variance of the scaled matrix
getVariance <- function(scaledCsvDataMatrix){
  
  numOfNodes <- ncol(scaledCsvDataMatrix)
  numOfSamples <- nrow(scaledCsvDataMatrix)
  
  varArray <- matrix( c(0), nrow = 1, ncol = numOfNodes) 
  
  for(j in 1:numOfNodes){
    sum <- 0
    for(i in 1:numOfSamples){
      sum <- sum + ((scaledCsvDataMatrix[i,j]-getMean(scaledCsvDataMatrix, j))*(scaledCsvDataMatrix[i,j]-getMean(scaledCsvDataMatrix, j)))
    }
    varArray[j] <- sum/(numOfSamples-1)
  }
  return(varArray)
}






#Function to compute density of the Feature matrix
computeDensity <- function(NMIDataMatrix, isSelected, clusterNumber){
  
  print("Compute Density")
  
  numOfNodes <- nrow(NMIDataMatrix)
  numOfCols <- ncol(NMIDataMatrix)
  
  numNodesInducedSet <- 0
  sumWeights <- 0
  
  for(i in 1:numOfNodes){
    if(isSelected[i] == 0){
      numNodesInducedSet <- numNodesInducedSet+1
    }
    for(j in i:numOfCols){
      if(numOfCols==numOfNodes && clusterNumber == 0){ # if condition added to accomodate rectangular mi matrix 
        if(isSelected[i]==0 && isSelected[j]==0 && NMIDataMatrix[i,j]<=1.00){ #if mi not normalized then keep threshold high
          sumWeights <- sumWeights + NMIDataMatrix[i,j]
        }
      }else{
        if(i<6 && j <6)
          if(isSelected[i]==0 && NMIDataMatrix[i,j]<=1.00){ #if mi not normalized then keep threshold high
            sumWeights <- sumWeights + NMIDataMatrix[i,j]
          }  
      }
    }
  }
  
  density <- sumWeights/numNodesInducedSet
  return(density)
}





#Function to compute induced degree
computeInducedDegree <- function(NMIDataMatrix, isSelected, clusterNumber){
  
  print("Compute Induced Degree")
  
  numOfNodes <- nrow(NMIDataMatrix)
  numofCols <- ncol(NMIDataMatrix)
  inducedDegree <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  
  sumDegreeWeights <- 0
  numInducedNodes <- 1
  
  for(i in 1:numOfNodes){
    if (isSelected[i] == 0){
      sumDegreeWeights <- 0
      for(j in 1:numofCols){
        if(numofCols==numOfNodes && clusterNumber == 0){ # if condition added to accomodate rectangular mi matrix 
          if(isSelected[j] == 0 && i!=j && NMIDataMatrix[i,j]<=1.00){ #threshold to be changed as above
            sumDegreeWeights <- sumDegreeWeights + NMIDataMatrix[i,j]
            numInducedNodes <- numInducedNodes+1
            print(paste("sum degree weights = ", sumDegreeWeights))
          }
        }else{
          if(NMIDataMatrix[i,j]<=1.00){ #threshold to be changed as above
            sumDegreeWeights <- sumDegreeWeights + NMIDataMatrix[i,j]
            numInducedNodes <- numInducedNodes+1
            print(paste("sum degree weights = ", sumDegreeWeights))
          }
        }  
      }
      inducedDegree[i] <- sumDegreeWeights
    }else{
      print(0)
      inducedDegree[i] <- 0
    }  
  }
  
  return(inducedDegree)
}




#Function to Compute Ranking (This function needs to be modified Later)
computeRanking <- function(inducedDegree, isShortListed, numOfNodes, rank){
  
  print("Compute Ranking")
  
  isToBeChecked <- matrix( c(0), nrow = 1, ncol = numOfNodes) 
  sizeCheckList <- 0
  r <- 0
  bestIndex <- 0
  
  for(i in 1:numOfNodes){
    if(isShortListed[i] == 1){
      isToBeChecked[i] <- 1
      sizeCheckList <- sizeCheckList +1
    }
  }
  
  for(i in 1:sizeCheckList){
    bestValue <- (-9999)
    for(j in 1:numOfNodes){
      if(isToBeChecked[j]==1){
        if(inducedDegree[j] > bestValue){
          print(paste("Before swapping rank[i] = ", rank[i] ," and rank[j] = ", rank[j]))
          bestIndex <- j
          bestValue <- inducedDegree[j]
          print(paste("After swapping rank[i] = ", rank[i]," and rank[j] = ", rank[j]))
        }
      }
    }
    r <- r+1
    rank[bestIndex] <- r
    isToBeChecked[bestIndex] <- 0
  }
  return(rank)
}



#Function to find the combinations
getComb <- function(n, t){
  val <- factorial(n)/(factorial(n-t)*factorial(t))
  return(val)
}





#function to find the cluster index
findingClusterIndex <- function(NMIDataMatrix, isSelected, clusterIndex, clusterNumber){
  
  print("Compute Cluster Index")
  
  numOfNodes <- nrow(NMIDataMatrix)
  numOfCols <- ncol(NMIDataMatrix)
  
  index <- 0
  for(i in 1:numOfNodes){
    if(isSelected[i] == 0){
      clusterIndex[i] <- index
      index <- index + 1
    }
  }
  
  print(paste("Index = ", index))
  
  for(i in 1:numOfNodes){
    if(isSelected[i] == 1){
      maxSimVal <- 0
      for(j in 1:numOfNodes){
        if(numOfCols == numOfNodes && clusterNumber == 0){ # if condition added to accomodate rectangular mi matrix 
          if(i!=j && isSelected[j] == 0){
            if(NMIDataMatrix[i,j] > maxSimVal){
              maxSimVal <- NMIDataMatrix[i,j]
              clusterIndex[i] <- clusterIndex[j]
            }
          }
        }else{
          correl <- cor(t(NMIDataMatrix[i,]), t(NMIDataMatrix[j,]))[1,1] 
          print(paste("Correlation value = ", correl))
          if(correl > maxSimVal){
            maxSimVal <- cor(t(NMIDataMatrix[i,]), t(NMIDataMatrix[j,])) #assigning set f1 to f1
            clusterIndex[i] <- clusterIndex[j]
          }
        }  
      }
    } 
  }
  return(clusterIndex)
}


#Function to obtain prototyping feature
obtainPrototypeFeature <- function(NMIDataMatrix, scaledCsvDataMatrix, isSelected, clusterIndex, K, clusterNumber){
  
  print("Compute Prototype Feature")
  
  numOfNodes <- nrow(NMIDataMatrix)
  numOfCols <- ncol(NMIDataMatrix)
  numOfSamples <- nrow(scaledCsvDataMatrix)
  
  simDegree <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  clusterPrototypeFeature <- matrix( c(0), nrow = 1, ncol = K )
  numClusterEle <- matrix( c(0), nrow = 1, ncol = K )
  
  for(i in 1:numOfNodes){
    numClusterEle[clusterIndex[i]+1] <- numClusterEle[clusterIndex[i]+1] + 1
  }
  
  for(i in 1:numOfNodes){
    for(j in 1:numOfCols){
      if(numOfCols==numOfNodes && clusterNumber == 0){ # if condition added to accomodate rectangular mi matrix 
        if(i != j && clusterIndex[i] == clusterIndex[j] ){
          simDegree[i] <- simDegree[i] + NMIDataMatrix[i,j]
        }
      }else{
        simDegree[i] <- simDegree[i] + NMIDataMatrix[i,j]
      }  
    }
  }
  
  print("Cluster Degree value of each feature")
  
  for(i in 1:numOfNodes){
    print(simDegree[i])
  }
  
  
  print("Number of cluster elements")
  for(i in 1:K){
    print(numClusterEle[i])
  }
  
  print("Cluster Prototype Features")
  


  for(index in 1:K){ 
    
    maxVal <- 0
    minVal <- 0
    
    for(j in 1:numOfNodes){
      if(clusterIndex[j] == index-1){
        if(simDegree[j] == 1){
          clusterPrototypeFeature[index] <- (j-1)
          break
        }
        if(clusterNumber == 1){
          varArray <- varArray1
        }else if(clusterNumber == 2){
          varArray <- varArray2
        }
        if(varArray[j] > maxVal){
          clusterPrototypeFeature[index] <- (j-1)
          maxVal <- varArray[j]
        }
      }
    }
    print(clusterPrototypeFeature[index])
  }
  

  print(paste("The number of rows in cluster Matrix  = ", numOfNodes))
  clusterNodes <- matrix(c(0), nrow=noOfClusters, ncol=numOfNodes)
  
  for(index in 1:K){
    print(paste("Cluster with index ",index," has the following features"))
    for(j in 1:numOfNodes){
      if(clusterIndex[j]==index-1){
        clusterNodes[index,j] <- j
        print(j)
      }  
    }
    
  }
  
  if(clusterNumber == 0){
    assign("clusterNodes", clusterNodes, envir = .GlobalEnv)
  }else if(clusterNumber == 1){
    assign("clusterNodes1", clusterNodes, envir = .GlobalEnv)
  }else if(clusterNumber == 2){
    assign("clusterNodes2", clusterNodes, envir = .GlobalEnv)
  }
  
  print(paste("Var is ", varArray[1]))
  
  for(j in 1:numOfNodes){
    isSelected[j] <- 1
    for(index in 1:K){
      if(clusterPrototypeFeature[index]==j-1){
        isSelected[j] <- 0
        break
      }  
    }
  }
  return(isSelected)
}




#Function to compute the density variation sequences (Needs some serious debugging)
computeDensityVarSeq <- function(NMIDataMatrix, scaledCsvDataMatrix, varArray, noOfClusters, clusterNumber){
  
  print("Compute Var Seq")
  
  K <- noOfClusters
  numOfNodes <- nrow(NMIDataMatrix)
  numOfCols <- ncol(NMIDataMatrix)
  numOfSamples <- nrow(scaledCsvDataMatrix)
  print(paste(numOfNodes,numOfSamples))
  
  isSelected <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  isCurrentlyDiscarded <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  isShortListed <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  isOptimal <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  inducedDegree <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  rank <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  clusterIndex <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  oldClusterIndex <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  
  
  currentDensity <- computeDensity(NMIDataMatrix, isSelected, clusterNumber)
  optimalDensity <- currentDensity
  
  clusterLoop <- 0
  optimalNodeCount <- 0
  
  print(paste("Optimal Density = ", optimalDensity))
  
  while(1){
    
    nodeCount <- 0
    
    for(i in 1:numOfNodes){
      if(isSelected[i]==0){
        nodeCount <- nodeCount + 1
      }
      rank[i] <- 999999
    }
    
    print("Starting the next iteration")
    print(paste("Total number of non-selected nodes is",nodeCount))
    
    if(nodeCount==0 || nodeCount <= K){
      break
    }
    
   
    inducedDegree <- computeInducedDegree (NMIDataMatrix, isSelected, clusterNumber)
   
    
    noShortListed <- 0
    for(i in 1:numOfNodes){
      if(isSelected[i] == 0){
        if(inducedDegree[i] >= 2*optimalDensity){
          isShortListed[i] <- 1
          noShortListed <- noShortListed + 1
        }else{
          isShortListed[i] <- 0
        }
      }else{
        isShortListed[i] <-0
      }
    }
    
    if(noShortListed > 0){
      print(paste("There is ",noShortListed," shortlisted Candidates " ))
    }
    
    rankLimit <- 2
    
    print(paste("The rank Limit is ", rankLimit))
    
    if(rankLimit > noShortListed){
      rankLimit <- 0.5*noShortListed
    }
    
    if(noShortListed == 1){
      rankLimit <- 2
    }else if(noShortListed == 0){
      break
    }
    print(paste("Second rankLimit is ", rankLimit))
    rank <- computeRanking(inducedDegree, isShortListed, numOfNodes, rank)
    
    for(i in 1:numOfNodes){
      isCurrentlyDiscarded[i] <- 0
      if(isShortListed[i] == 1){
        print("Short Listed")
        if(rank[i] < rankLimit){
          print(paste("Second rankLimit is ", rankLimit))
          isSelected[i] <- 1
          isCurrentlyDiscarded[i] <- 1
        }
      }
    }
    
    currentDensity <- computeDensity(NMIDataMatrix, isSelected, clusterNumber)
    
    density <- 0
    noNodes <- 0
    sumWeight <- 0
    print(paste("Current Density is ", currentDensity))
    
    for(i in 1:numOfNodes){
      if(isSelected[i] == 1 && isCurrentlyDiscarded[i] == 0){
        noNodes <- 0
        sumWeight <- 0
        isSelected[i] <- 0
        
        for(j in 1:numOfCols){
          if(numOfCols == numOfNodes && clusterNumber==0){ #modified for multiview
            if(isSelected[j]==0 && j!=i && (NMIDataMatrix[i,j]<=1.00)){
              noNodes <- noNodes + 1
              sumWeight <- sumWeight + NMIDataMatrix[i,j]
            }
          }else{
            if((NMIDataMatrix[i,j]<=1.00)){
              noNodes <- noNodes + 1
              sumWeight <- sumWeight + NMIDataMatrix[i,j]
            }
          }
        }
        
        density <- (currentDensity*noNodes+sumWeight)/(noNodes+1)
        
        if(density < currentDensity){
          print(paste("One node with id ",i+1," is being readded"))
          currentDensity <- density
        }else{
          isSelected[i] <- 1
        }
      }
    }
    
    if(currentDensity <= optimalDensity){
      optimalDensity <- currentDensity
      optimalNodeCount <- 0
      
      for (i in 1:numOfNodes){
        isOptimal[i]= 1-isSelected[i];
      }
      
      for(i in 1:numOfNodes){
        if(isOptimal[i] == 1){
          optimalNodeCount <- optimalNodeCount + 1
          print(i+1)
        }
      }
      print("")
    }
    print(paste("Optimal Node Count is ", optimalNodeCount))
    print(paste("Optimal Density is ", optimalDensity))
    
    if(optimalNodeCount==K)
      break
  }
  
  
  while(1){
    
    clusterIndex <- findingClusterIndex(NMIDataMatrix, isSelected, clusterIndex, clusterNumber)
    
    print("Cluster Index")
    for(i in 1:numOfNodes){
      print(clusterIndex[i])
    }
    
    continueFlag <- 0
    for(i in 1:numOfNodes){
      if(clusterIndex[i] != oldClusterIndex[i]){
        continueFlag <- 1
        print("First Break")
        break
      }
    }
    
    if(continueFlag == 0){
      print("Second Break")
      break
    }
    
    for(i in 1:numOfNodes){
      oldClusterIndex[i] <- clusterIndex[i]
    }
    isSelected <- obtainPrototypeFeature(NMIDataMatrix, scaledCsvDataMatrix, isSelected, clusterIndex, K, clusterNumber)
    clusterLoop <-clusterLoop +1
    print(paste(" cluster loop = ", clusterLoop))
  }
}


