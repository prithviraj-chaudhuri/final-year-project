#Function to read Data
readData <- function(NMIfilename, isHeader){
  data <- read.table(NMIfilename, header=isHeader)
  return(data)
}


#function to scale the CSV Data
scaleing <- function(CSVDataMatrix, numOfNodes, numOfSamples){
  
  minScaledVal = 0;
  maxScaledVal = 1;
  minVal = 0;
  maxVal = 0;
  scaledCsvDataMatrix = 0;
  
  for(j in 0:numOfNodes){
    minVal = CSVDataMatrix[j]
    maxVal = CSVDataMatrix[j]
    
    for(i in 0:numOfSamples){
      scaledCsvDataMatrix[i*numOfNodes+j] <- ((CSVDataMatrix[i*numOfNodes+j]-minVal)/(maxVal-minVal))*(maxScaledVal-minScaledVal)+minScaledVal;
    }
    
  }
  
  return(scaledCsvDataMatrix);
}


#Exectution of the code

#reading the data
NMIDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\facebook\\0.feat", FALSE)
CSVDataMatrix <- readData("C:\\Users\\prithviraj\\Desktop\\final-year-project\\facebook\\facebook_combined.txt", TRUE)

#setting number of nodes || features
numOfNodes <- ncol(NMIData)-1
numOfSamples <- 5673

#scaled csv data matrix
scaledCsvDataMatrix <- scaleing(CSVDataMatrix, numOfNodes, numOfSamples)
