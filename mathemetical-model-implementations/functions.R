

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
getVariance <- function(scaledCsvDataMatrix, numOfNodes, numOfSamples){
  
  numOfNodes <- ncol(scaledCsvDataMatrix)
  numOfSamples <- nrow(scaledCsvDataMatrix)
  
  varArray <- matrix( c(0), nrow = 1, ncol = numOfNodes) 
  
  for(j in 1:numOfNodes){
    sum <- 0
    for(i in 1:numOfSamples){
      sum <- sum + ((scaledCsvDataMatrix[i,j]-getMean(scaledCsvDataMatrix, j))*(scaledCsvDataMatrix[i,j]-getMean(scaledCsvDataMatrix, j)))
    }
    varArray[j] <- sum/numOfSamples
  }
  return(varArray)
}






#Function to compute density of the Feature matrix
computeDensity <- function(NMIDataMatrix, isSelected){
  
  numOfNodes <- ncol(NMIDataMatrix)
  
  numNodesInducedSet <- 0
  sumWeights <- 0
  
  for(i in 1:numOfNodes){
    if(isSelected[i] == 0){
      numNodesInducedSet <- numNodesInducedSet+1
    }
    for(j in i:numOfNodes){
      if(isSelected[i]==0 && isSelected[j]==0 && NMIDataMatrix[i,j]<=1.00){
        sumWeights <- sumWeights + NMIDataMatrix[i,j]
      }
    }
  }
  
  density <- sumWeights/numNodesInducedSet
  return(density)
}





#Function to compute induced degree
computeInducedDegree <- function(NMIDataMatrix, isSelected){
  
  numOfNodes <- ncol(NMIDataMatrix)
  inducedDegree <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  
  sumDegreeWeights <- 0
  numInducedNodes <- 1
  
  for(i in 1:numOfNodes){
    if (isSelected[i] == 0){
      sumDegreeWeights <- 0
      for(j in 1:numOfNodes){
        if(isSelected[j] == 0 && i!=j && NMIDataMatrix[i,j]<=1.00){
          sumDegreeWeights <- sumDegreeWeights + NMIDataMatrix[i,j]
          numInducedNodes <- numInducedNodes+1
          print(paste("sum degree weights = ", sumDegreeWeights))
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
computeRanking <- function(inducedDegree, isShortListed, numOfNodes){
  
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
          print(paste("\nBefore swapping rank[i] = ", rank[i] ," and rank[j] = ", rank[j]))
          bestIndex <- j
          bestValue <- inducedDegree[j]
          print(paste("\nAfter swapping rank[i] = ", rank[i]," and rank[j] = ", rank[j]))
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
findingClusterIndex <- function(NMIDataMatrix, isSelected, clusterIndex){
  
  numOfNodes <- ncol(NMIDataMatrix)
  
  index <- 0
  for(i in 1:numOfNodes){
    if(isSelected[i] == 0){
      clusterIndex[i] <- index
      index <- index + 1
    }
  }
  
  print(paste("\nIndex = ", index))
  
  for(i in 1:numOfNodes){
    if(isSelected[i] == 1){
      maxSimVal <- 0
      for(j in 1:numOfNodes){
        if(i!=j && isSelected[j] == 0){
          if(NMIDataMatrix[i,j] > maxSimVal){
            maxSimVal <- NMIDataMatrix[i,j]
            clusterIndex[i] <- clusterIndex[j]
          }
        }
      }
    } 
  }
  return(clusterIndex)
}



#Function to obtain prototyping feature
obtainPrototypeFeature <- function(NMIDataMatrix, scaledCsvDataMatrix, isSelected, clusterIndex){
  
  K <- 6
  
  numOfNodes <- ncol(scaledCsvDataMatrix)
  numOfSamples <- nrow(scaledCsvDataMatrix)
  
  simDegree <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  clusterPrototypeFeature <- matrix( c(0), nrow = 1, ncol = K )
  numClusterEle <- matrix( c(0), nrow = 1, ncol = K )
  
  for(i in 1:numOfNodes){
    numClusterEle[clusterIndex[i]] <- numClusterEle[clusterIndex[i]] + 1
  }
  
  for(i in 1:numOfNodes){
    for(j in 1:numOfNodes){
      if(i != j && clusterIndex[i] == clusterIndex[j] ){
        simDegree[i] <- simDegree[i] + NMIDataMatrix[i,j]
      }
    }
  }
  
  print("\nCluster Degree value of each feature\n")
  
  for(i in 1:numOfNodes){
    print(simDegree[i])
  }
  
  for(i in 1:K){
    print(numClusterEle[i])
  }
  
  print("\nCluster Prototype Features")
  
  for(index in 1:K){
    
    maxVal <- 0
    minVal <- 0
    
    for(j in 1:numOfNodes){
      if(clusterIndex[j] == index){
        if(simDegree[j] == 0){
          clusterPrototypeFeature[index] <- j
          break
        }
        
        if(varArray[j] > maxVal){
          clusterPrototypeFeature[index] <- j
          maxVal <- varArray[j]
        }
      }
    }
    print(clusterPrototypeFeature[index])
  }
  
  
  for(index in 1:K){
    print(paste("\nCluster with index ",index," has the following features\n"))
    for(j in 1:numOfNodes){
      if(clusterIndex[j]==index)
        print(j)
    }
    
  }
  print(paste("Var is ", varArray[1]))
  
  for(j in 1:numOfNodes){
    isSelected[j] <- 1
    for(index in 1:K){
      if(clusterPrototypeFeature[index]==j){
        isSelected[j] <- 0
        break
      }  
    }
  }
  return(isSelected)
}




#Function to compute the density variation sequences (Needs some serious debugging)
computeDensityVarSeq <- function(NMIDataMatrix, scaledCsvDataMatrix, maximizer, time, dVal){
  
  K <- 6 #need to know what this means
  
  numOfNodes <- ncol(scaledCsvDataMatrix)
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
  
  currentDensity <- computeDensity(NMIDataMatrix, isSelected)
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
      rank[i] <- 9999
    }
    
    if(nodeCount==0 || nodeCount <= K){
      break
    }
    
   
    inducedDegree <- computeInducedDegree (NMIDataMatrix, isSelected)
   
    
    
    
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
    
    rankLimit <- 1
    
    print(paste("The rank Limit is ", rankLimit))
    
    if(rankLimit > noShortListed){
      rankLimit <- 0.5*noShortListed
    }
    
    if(noShortListed == 1){
      rankLimit <- 1
    }else if(noShortListed == 0){
      break
    }
    
    rank <- computeRanking(inducedDegree, isShortListed, numOfNodes)
    
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
    
    currentDensity <- computeDensity(NMIDataMatrix, isSelected)
    
    density <- 0
    noNodes <- 0
    sumWeight <- 0
    print(paste("Current Density is ", currentDensity))
    
    for(i in 1:numOfNodes){
      if(isSelected[i] == 1 && isCurrentlyDiscarded[i] == 0){
        noNodes <- 0
        sumWeight <- 0
        isSelected[i] <- 0
        
        for(j in 1:numOfNodes){
          if(isSelected[j]==0 && j!=i && (NMIDataMatrix[i,j]<=1.00)){
            noNodes <- noNodes + 1
            sumWeight <- sumWeight + NMIDataMatrix[i,j]
          }
        }
        
        density <- (currentDensity*noNodes+sumWeight)/(noNodes+1)
        
        if(density < currentDensity){
          print(paste("One node with id ",i+1," is being read"))
          currentDensity <- density
        }else{
          isSelected[i] <- 1
        }
      }
    }
    
    if(currentDensity <= optimalDensity){
      optimalDensity <- currentDensity
      optimalNodeCount <- 0
      
      for(i in 1:numOfNodes){
        if(isOptimal[i] == 1){
          optimalNodeCount <- optimalNodeCount + 1
          print(i+1)
        }
      }
      print("\n")
    }
    print(paste("\nOptimal Node Count is ", optimalNodeCount))
    print(paste("\nOptimal Density is ", optimalDensity))
    
    if(optimalNodeCount==K)
      break
  }
  
  
  while(1){
    
    clusterIndex <- findingClusterIndex(NMIDataMatrix, isSelected, clusterIndex)
    continueFlag <- 0
    for(i in 1:numOfNodes){
      if(clusterIndex[i] != oldClusterIndex[i]){
        continueFlag <- 1
        print("\nFirst Break")
        break
      }
    }
    
    if(continueFlag == 0){
      print("\nSecond Break")
      break
    }
    
    for(i in 1:numOfNodes){
      oldClusterIndex[i] <- clusterIndex[i]
    }
    isSelected <- obtainPrototypeFeature(NMIDataMatrix, scaledCsvDataMatrix, isSelected, clusterIndex)
    clusterLoop <-clusterLoop +1
    print(paste("\n cluster loop = ", clusterLoop))
  }
  
  time <- t
}

