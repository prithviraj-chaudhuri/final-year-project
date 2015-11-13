#This Program is for feature selection
#Author - Prithviraj Chaudhuri

#including required files

source("functions-multiview.R")


#Exectution of the code

#reading the data
feature1 <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test3.feature", FALSE)
feature2 <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test3-2.feature", FALSE)
NMIDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test-multiview.mi", FALSE)

numOfNodes1 <- nrow(NMIDataMatrix)
numOfNodes2 <- ncol(NMIDataMatrix)

#scaled csv data matrix
scaledFeature1 <- scaleing(feature1)
scaledFeature2 <- scaleing(feature2)

#calculating the variance of each of the feature vectors
varArray1 <- getVariance(scaledFeature1)
varArray2 <- getVariance(scaledFeature2)

#The number of clusters
noOfClusters <- 3

#The arrangement of the features
clusterNodes1 <- matrix(c(0), nrow=noOfClusters, ncol=numOfNodes1)
clusterNodes2 <- matrix(c(0), nrow=noOfClusters, ncol=numOfNodes2)


#assigning ti clusters
computeDensityVarSeq(NMIDataMatrix, scaledFeature1, varArray1, noOfClusters, 1)
NMIDataMatrix2 <- as.data.frame(t(NMIDataMatrix))
row.names(NMIDataMatrix2) <- seq(nrow(NMIDataMatrix2))
computeDensityVarSeq(NMIDataMatrix2, scaledFeature2, varArray2, noOfClusters, 2)