
#This Program is for feature selection
#Author - Prithviraj Chaudhuri

#including required files

source("functions.R")


#Exectution of the code

#reading the data
NMIDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test.mi", FALSE)
CSVDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\test-data\\test.feature", TRUE)

numOfNodes <- ncol(NMIDataMatrix)
maximizerIndex <- matrix( c(0), nrow=1, ncol=numOfNodes)
dVal <- matrix( c(0), nrow=1, ncol=numOfNodes)
#creating the MI
#NMIDataMatrix <- getNMI(CSVDataMatrix)




#scaled csv data matrix
scaledCsvDataMatrix <- scaleing(CSVDataMatrix)

#calculating the variance of each of the feature vectors
varArray <- getVariance(scaledCsvDataMatrix)

#calculating the density of variance sequence
computeDensityVarSeq (NMIDataMatrix, scaledCsvDataMatrix, maximizerIndex, 0, dVal)


