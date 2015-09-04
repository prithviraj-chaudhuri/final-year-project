
#This Program is for feature selection
#Author - Prithviraj Chaudhuri

#Loading Libraries
#library("entropy")


#Function to read Data
readData <- function(NMIfilename, isHeader){
  data <- read.table(NMIfilename, header=isHeader)
  return(data)
}


#Function to get MI from CSV Data
#getNMI <- funtion(CSVDataMatrix){
#  NMIDataMatrix <- makemim(CSVDataMatrix)
#  return(NMIDataMatrix)
#}


#function to scale the CSV Data
scaleing <- function(CSVDataMatrix, numOfNodes, numOfSamples){
  
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
getMean <- function(scaledCsvDataMatrix, col, numOfSamples){
  
  sum <- 0
  for(i in 1:numOfSamples){
    sum <- sum + scaledCsvDataMatrix[i,col]
  }
  mean <- sum/numOfSamples
  return(mean)
}


#Function to calculate the variance of the scaled matrix
getVariance <- function(scaledCsvDataMatrix, numOfNodes, numOfSamples){
  
  varArray <- matrix( c(0), nrow = 1, ncol = numOfNodes) 
  
  for(j in 1:numOfNodes){
    sum <- 0
    for(i in 1:numOfSamples){
      sum <- sum + ((scaledCsvDataMatrix[i,j]-getMean(scaledCsvDataMatrix, j, numOfSamples))*(scaledCsvDataMatrix[i,j]-getMean(scaledCsvDataMatrix, j, numOfSamples)))
    }
    varArray[j] <- sum/numOfSamples
  }
  return(varArray)
}


#Function to compute the density variation sequences
computeDensityVarSeq <- function(NMIDataMatrix, scaledCsvDataMatrix, maximizer, time, dVal, numOfNodes){
  
  isSelected <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  isCurrentlyDiscarded <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  isShortlisted <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  isOptimal <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  inducedDegree <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  rank <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  clusterIndex <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  oldClusterIndex <- matrix( c(0), nrow = 1, ncol = numOfNodes )
  
  currentDensity <- computeDensity(NMIDataMatrix, isSelected)
  optimalDensity <- currentDensity
  
  sprintf("Optimal Density = %f", optimalDensity)
  
  while(1){
    
    nodeCount <- 0
    
    for(i in 1:numOfNodes){
      if(isSelected[i]==0){
        nodeCount++
      }
      rank[i]= <- 999999
    }
   
    if(nodeCount==0 || nodeCount <= K){
      break
    }
    
    for(i in 1:numOfNodes){
      if (isSelected[i] == 0){
        inducedDegree[i] <- computeInducedDegree (NMIDataMatrix, isSelected, i)
      }else{
        inducedDegree[i] <- 0
      }  
    }
    
    noShortListed <- 0
    for(i in 1:numOfNodes){
      if(isSelected[i] == 0){
        if(inducedDegree[i] >= 2*optimalDensity){
          isShortListed[i] <- 1
          noShortListed++
        }else{
          isShortListed[i] <- 0
        }
      }else{
        isShortListed[i] <-0
      }
    }
    
    if(noShortListed > 0){
      sprintf("There is %d shortlisted Candidates ", noShortListed)
    }
    
    rankLimit <- 1
    
    sprintf("The rank Limit is %d ", rankLimit)
    
    if(rankLimit > noShortListed){
      rankLimit <- 0.5*noShortListed
    }
      
    if(noShortListed == 1){
      rankLimit <- 1
    }else if(noShortListed == 0){
      break;
    }
    
    computeRanking(inducedDegree, isShortListed, rank)
    
    for(i in 1:numOfNodes){
      isCurrentlyDiscarded[i] <- 0
      if(isShortListed[i] == 1){
        if(rank[i] < rankLimit){
          sprintf("Second rankLimit is %d", rankLimit)
          isSelected[i] <- 1
          isCurrentlyDiscarded[i] <- 1
        }
      }
    }
    
    
    
  }
}



#Exectution of the code

#reading the data
NMIDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test.mi", FALSE)
CSVDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test.feature", TRUE)

#creating the MI
#NMIDataMatrix <- getNMI(CSVDataMatrix)

#setting number of nodes || features
numOfNodes <- ncol(CSVDataMatrix)
numOfSamples <- nrow(CSVDataMatrix)

#scaled csv data matrix
scaledCsvDataMatrix <- scaleing(CSVDataMatrix, numOfNodes, numOfSamples)

#calculating the variance of each of the feature vectors
varArray <- getVariance(scaledCsvDataMatrix, numOfNodes, numOfSamples)
