#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define NUMOFSAMPLES 5673 // number of samples
#define NUMOFNODES 2400 // number of features
#define K 6 
//#define NUMOFNODES 4   // number of features
//#define K 2
#define THRESHOLD 1.00
#define EPSILON 0.01
#define ALPHA 0
#define BETA 1

void readNMIData (char *, double *);
void readCSVData (char *, double *);
void scaleing(double *, double *);
double mean(double *,int);
void variance(double *);

void computeDensityVarSeq (double *, double *, int *,int *,double *);
double computeDensity(double *, int *);
double computeInducedDegree (double *, int *, int);
void computeRanking (double *, int *, int *);
void swap (int *, int *);
void findingClusterIndex(double *NMIDataMatrix, int *isSelected, int *clusterIndex);
void obtainingPrototypeFeature(double *NMIDataMatrix, double *scaledcsvDataMatrix, int *isSelected, int *clusterIndex);


double vararray[NUMOFNODES];

int main(int argc, char *argv[])
{
	char NMIFileName[100];
	char csvFileName[100];
	
	double *NMIDataMatrix;
	double *csvDataMatrix;
	double *scaledcsvDataMatrix;
	//double NMIDataMatrix[NUMOFNODES][NUMOFNODES] = {{0,1,1,1,0,0,0,0,0},{1,0,1,1,0,0,0,0,0},{1,1,0,1,0,0,0,0,0},{1,1,1,0,1,0,0,0,0},{0,0,0,1,0,1,1,1,1},{0,0,0,0,1,0,1,1,1},{0,0,0,0,1,1,0,1,1},{0,0,0,0,1,1,1,0,1},{0,0,0,0,1,1,1,1,0}};
	int maximizerIndex[NUMOFNODES];
	double dVal[NUMOFNODES];
	double rateVal[NUMOFNODES];
	int *whetherCoreNode;
	int *coreSetIndex;
	int *clusterCoreIndex;
	int i,j,c,n;
	int time;
	int numCoreNode = 0;
	int numCoreSet;
	int numClusterCore;
	int option;
	
	if(argc!=3)
	{
		fprintf(stderr,"\n%s nmiFile csvFile\n\n",argv[0]);
                exit(1);
        }

	strcpy(NMIFileName,argv[1]);
	strcpy(csvFileName,argv[2]);
	NMIDataMatrix = (double *)calloc(NUMOFNODES*NUMOFNODES,sizeof(double));
	csvDataMatrix = (double *)calloc(NUMOFSAMPLES*NUMOFNODES,sizeof(double));
	scaledcsvDataMatrix = (double *)calloc(NUMOFSAMPLES*NUMOFNODES,sizeof(double));
	whetherCoreNode = (int *)calloc(NUMOFNODES,sizeof(int));
	coreSetIndex = (int *)calloc(NUMOFNODES,sizeof(int));
	clusterCoreIndex = (int *)calloc(NUMOFNODES,sizeof(int));
	
	/*for(i=0;i<NUMOFNODES;i++)
		for(j=0;j<NUMOFNODES;j++)
			if(i!=j)
				NMIDataMatrix[i][j]=(i+1)*(j+1); 	
			else	
				NMIDataMatrix[i][j]=0;*/
	//NMIDataMatrix[NUMOFNODES][NUMOFNODES] = ((1,1),(2,2)); 
	//printf("\n%d %d",NMIDataMatrix[0][0],NMIDataMatrix[8][8]);	

	readNMIData(NMIFileName, NMIDataMatrix);
	readCSVData(csvFileName, csvDataMatrix);
        printf("data read\n");
	scaleing(csvDataMatrix, scaledcsvDataMatrix);
        

	variance(scaledcsvDataMatrix);
	computeDensityVarSeq (NMIDataMatrix, scaledcsvDataMatrix, maximizerIndex, &time, dVal);
	/*identifyCoreNodes (maximizerIndex, dVal, time, rateVal, whetherCoreNode, &numCoreNode);
	makeCoreSet (maximizerIndex, whetherCoreNode, rateVal, time, coreSetIndex, &numCoreSet);
	partitionCoreSetIntoClusterCore (NMIDataMatrix, coreSetIndex, &numCoreSet, clusterCoreIndex, &numClusterCore);
	expandCoreGroupIntoCluster (NMIDataMatrix, whetherCoreNode, clusterCoreIndex, numClusterCore, maximizerIndex, time);*/
	
	while(1)
	{
		printf("\nEnter 1. To show Core Nodes, 2. To Show Cluster Core, 3. To exit");
		scanf("%d",&option);
		switch(option)
		{
			case 1:
				printf("\n\nThe core nodes are :");
				for(i=0;i<NUMOFNODES;i++)
				{
					if(whetherCoreNode[i]==1)
						printf(" %d", i);
				}
				printf("\n");
				printf("\nTotal number of core nodes is %d\n", numCoreNode);
				break;
			case 2:
				for(c=0;c<numClusterCore;c++)
				{
					printf("\nThe nodes in cluster %d are:",c);
					for(n=0;n<NUMOFNODES;n++)
					{
						if(clusterCoreIndex[n]==c)
							printf(" %d",n);
					}
				}
				printf("\n");
				break;
			case 3:
				exit(1);
	
		}
	}

	//for(i=0;i<NUMOFNODES;i++)
	//	printf(" maximizerIndex[%d] = %d", i, maximizerIndex[i]);
	
	printf("\nTime = %d",time);	
	//for(i=0; i<time;i++)
	//	printf(" dVal[%d] = %lf", i, dVal[i]);
	
	//for(i=0;i<NUMOFNODES;i++)
	//	printf("\nCoreSetIndex[%d] = %d", i, coreSetIndex[i]);
		
}

void readNMIData (char *fileName, double *NMIDataMatrix)
{
	FILE *fp;
	int i, j;
	double number;
	
	fp = fopen(fileName,"r");
	if(fp==NULL)
	{
		printf("\nUnable to open the input file.\n");
		exit(1);
	}

	//considering the class label info
	/*for (i=0; i<NUMOFNODES+1; i++)
	{
		for (j=0; j<NUMOFNODES+1; j++)
		{
			if (i!=0 && j!=0)
			{
				fscanf(fp, "%lf", &NMIDataMatrix[(i-1)*NUMOFNODES+(j-1)]);
			}
			else
			{
				fscanf(fp, "%lf", &number);
			}
		}
	}*/
	// considering without the class label info
	for (i=0; i<NUMOFNODES; i++)
	{
		for (j=0; j<NUMOFNODES; j++)
		{
			fscanf(fp, "%lf", &NMIDataMatrix[i*NUMOFNODES+j]);
			//if(i==499&&j==499)
			//printf("i=%d,j=%d,val=%lf,", i,j,NMIDataMatrix[i*NUMOFNODES+j]); 
		}
	}
	fclose(fp);
}

void readCSVData (char *csvFileName, double *csvDataMatrix)
{
	FILE *fp;
	int i, j;
	int class;
	char header[20000];
	
	fp = fopen(csvFileName,"r");
	if(fp==NULL)
	{
		printf("\nUnable to open the input file.\n");
		exit(1);
	}

	//reading the header file
	fscanf(fp, "%s",header);
	printf("\nHeader is %s\n",header);

	//considering the class label info
	for (i=0; i<NUMOFSAMPLES; i++)
	{
		for (j=0; j<(NUMOFNODES+1); j++)
		{
			if(j==0)
			{
				fscanf(fp, "%d",&class);
			}
			else	
				fscanf(fp, ",%lf", &csvDataMatrix[i*NUMOFNODES+(j-1)]);
		}
	}
	fclose(fp);
}

void scaleing(double *csvDataMatrix, double *scaledcsvDataMatrix)
{
	double minVal;
	double maxVal;
	double minScaledVal=0;
	double maxScaledVal=1;	
	
	int i,j;

	for(j=0;j<NUMOFNODES;j++)
	{
		minVal = csvDataMatrix[j];
		maxVal = csvDataMatrix[j];
		for(i=1;i<NUMOFSAMPLES;i++)
		{
			if(csvDataMatrix[i*NUMOFNODES+j] < minVal)
				minVal = csvDataMatrix[i*NUMOFNODES+j];
			else if(csvDataMatrix[i*NUMOFNODES+j] > maxVal)
				maxVal = csvDataMatrix[i*NUMOFNODES+j];
		}
		
		for(i=0;i<NUMOFSAMPLES;i++)
		{
			scaledcsvDataMatrix[i*NUMOFNODES+j]=((csvDataMatrix[i*NUMOFNODES+j]-minVal)/(maxVal-minVal))*(maxScaledVal-minScaledVal)+minScaledVal;
		}
		/*if(j==3)
		{
			for(i=0;i<NUMOFSAMPLES;i++)
				printf("%lf,",scaledcsvDataMatrix[i*NUMOFNODES+j]); 
		}*/
         //  printf("scale \n");
	}
		
}

// This function compute the density variation sequences
void computeDensityVarSeq (double *NMIDataMatrix, double *scaledcsvDataMatrix, int *maximizerIndex, int *time, double *dVal)
{
	int *isSelected;
	int *isCurrentlyDiscarded;
	int *isShortlisted;
	int noShortlisted;
	int *isOptimal;
	int *rank;
	int i,j;
	int t=0;
	int continueFlag;
	double *weight;
	double sumWeight;	
	double *inducedDegree;
	int nodeCount;
	int optimalNodeCount;
	int noNodes;
	int numWeight;
	int ii;
	int rankLimit;
	double density, currentDensity, optimalDensity;
	int *clusterIndex;
	int *oldClusterIndex;
        int clusterloop=0;

	isSelected = (int *)calloc(NUMOFNODES,sizeof(int));
	isCurrentlyDiscarded = (int *)calloc(NUMOFNODES,sizeof(int));
	isShortlisted = (int *)calloc(NUMOFNODES,sizeof(int));
	isOptimal = (int *)calloc(NUMOFNODES,sizeof(int));
	inducedDegree = (double *)calloc(NUMOFNODES,sizeof(double));
	rank = (int *)calloc(NUMOFNODES,sizeof(int));
	clusterIndex= (int *)calloc(NUMOFNODES,sizeof(int));
	oldClusterIndex= (int *)calloc(NUMOFNODES,sizeof(int));

	//for(ii=0;ii<5;ii++)

	printf("1");
	currentDensity = computeDensity(NMIDataMatrix, isSelected);
	//currentDensity = density;
	optimalDensity = currentDensity;
	printf("2");
	printf("\nOptimal Density = %f", optimalDensity);	
	for(;;)
	{
		nodeCount = 0;
		for (i=0; i<NUMOFNODES; i++)
		{
			if (isSelected[i] == 0)
				nodeCount ++;
			rank[i]=999999;
		}
		printf("\nStarting of next iteration");
		printf("\nTotal number of non-selected nodes is %d\n", nodeCount);
		if(nodeCount == 0)
			break;
		//For KM
		if(nodeCount <= K)
			break;
	
		for (i=0; i<NUMOFNODES; i++)
		{
			if (isSelected[i] == 0)
			{
				inducedDegree[i] = computeInducedDegree (NMIDataMatrix, isSelected, i);
			}
			else
				inducedDegree[i]=0;
		}
	
		noShortlisted = 0;	
		for (i=0; i<NUMOFNODES; i++)
		{
			if (isSelected[i] == 0)
			{
				//if (inducedDegree[i] >= 2*(1+EPSILON)*optimalDensity)
				//if (inducedDegree[i] >= (1+EPSILON)*optimalDensity)
				if (inducedDegree[i] >= 2*optimalDensity)
				{
					isShortlisted[i] = 1;
					noShortlisted ++;	
				}
				else
					isShortlisted[i] = 0;
			}
			else
				isShortlisted[i] = 0;
		}
		
		if(noShortlisted > 0)	
			printf("\nThere is %d shortlisted candidate", noShortlisted);
		
		/*rankLimit = (int)((EPSILON/(1+EPSILON))*nodeCount);
		//rankLimit = (rankLimit>0 ? rankLimit: 1);
		rankLimit = (rankLimit>noShortlisted ? 0.5*noShortlisted: rankLimit);
		printf("\nRankLimit is %d", rankLimit);*/
		
		//rankLimit = (int)(0.5*sqrt(nodeCount));
		//rankLimit = pow(pow(100,0.5),0.5);
		rankLimit = 1;
		printf("\nRank Limit is %d",rankLimit);	
		
		//if(nodeCount<=100)
		//	rankLimit = 2;
			
		
		if(rankLimit>noShortlisted)
			rankLimit = 0.5*noShortlisted;
		
		if(noShortlisted==1)
			rankLimit=1;
		else if(noShortlisted==0)
			break;
		
		//if(rankLimit == 0)
		//	break;	
	
		computeRanking (inducedDegree, isShortlisted, rank);	
		
		for (i=0; i<NUMOFNODES; i++)
		{
			isCurrentlyDiscarded[i] = 0;
			if (isShortlisted[i] == 1)
			{
				if (rank[i] < rankLimit)
				{
					printf("\niSecond time RankLimit is %d", rankLimit);
					isSelected[i] = 1;
					isCurrentlyDiscarded[i] = 1;
				}
				
			}
		}
		
		currentDensity = computeDensity(NMIDataMatrix, isSelected);
		
		printf("\nCurrent Density is %lf\n", currentDensity);	
	
		
		//if(rankLimit!=1)		
		for(i=0; i<NUMOFNODES; i++)
		{
			if(isSelected[i] == 1 && isCurrentlyDiscarded[i] == 0)
			{
				noNodes = 0;
				sumWeight = 0;
				isSelected[i] = 0;
				for(j=0; j<NUMOFNODES; j++)
				{
					if(isSelected[j]==0 && j!=i && (NMIDataMatrix[i*NUMOFNODES+j]<=THRESHOLD))
					{
						noNodes++;
						sumWeight += NMIDataMatrix[i*NUMOFNODES+j];
					}
				}
				//version 1
				density = (currentDensity*noNodes+sumWeight)/(noNodes+1);
				//version 2 
				//density = (currentDensity*noNodes+sumWeight)/(comb(noNodes,2)+noNodes);
				
				//printf("\nNumber of nodes is %d",noNodes);
				//printf("\nDensity1 is %lf",density);	
				//density = computeDensity(NMIDataMatrix, isSelected);
				//printf("\nDensity2 is %lf",density);	
				//density = ((currentDensity*)+)/
				if(density<currentDensity)
				{
					//printf("\nDensity is %lf",density);
					//break;
					printf("\n One node with id %d is being readded",i+1);
					currentDensity = density;
					;
				}
				else
					isSelected[i] = 1;	
			}
		}

		if(currentDensity <= optimalDensity)
		{
			optimalDensity = currentDensity;
			optimalNodeCount = 0;
			for (i=0; i<NUMOFNODES; i++)
			{
				isOptimal[i]= 1-isSelected[i];
			}
			
			for (i=0; i<NUMOFNODES; i++)
			{
				if(isOptimal[i]==1)
				{
					optimalNodeCount ++;
					printf("%d,", i+1);
				}
			}	
			printf("\n");
		}

		printf("\nOptimal Node Count is %d\n", optimalNodeCount);
		printf("\nOptimal Density is %f", optimalDensity);
		if(optimalNodeCount==K)
			break;
		//compute the nodes with min weight value 	
		
		/*for(i=0;i<NUMOFNODES;i++)
		{
			if(isSelected[i]==0)
			{
				if(weight[i]==minWeight)
				{
					maximizerIndex[i]=t;	
					isSelected[i]=1;
				}
			}
		}
	
		//compute the d values	
		printf("\nMinWeight = %lf",minWeight);
		//dVal[t++] = (double) minWeight/nodeCount;
		dVal[t++] = (double) minWeight;*/
		
	};
		

	for(;;)
	{	
		findingClusterIndex(NMIDataMatrix, isSelected, clusterIndex);
		continueFlag=0;
		for(i=0;i<NUMOFNODES;i++)
		{
			if(clusterIndex[i]!=oldClusterIndex[i])
			{
				continueFlag=1;
				printf("\nFirst break");
				break;
			}
		}
		if(continueFlag==0)
		{
			printf("\nSecond break");
			break;
		}
		for(i=0;i<NUMOFNODES;i++)
		{
			oldClusterIndex[i]=clusterIndex[i];
		}
		obtainingPrototypeFeature(NMIDataMatrix, scaledcsvDataMatrix, isSelected, clusterIndex);
                clusterloop++;
                printf("%d\n",clusterloop);
                
	}
	*time=t;
	
	free(isSelected);
	free(isOptimal);
	free(inducedDegree);
	free(rank);
	//return 0;
}


void findingClusterIndex(double *NMIDataMatrix, int *isSelected, int *clusterIndex)
{
	int i,j,index=0;
	double maxSimVal;

	//assigning the cluster index for each selected feature 	
	for(i=0;i<NUMOFNODES;i++)
	{
		if(isSelected[i]==0)
			clusterIndex[i]=index++;
	}

	printf("\nIndex = %d",index);	
	//assigning the cluster index for each discarded feature	
	for(i=0;i<NUMOFNODES;i++)
	{
		if(isSelected[i]==1)
		{
			maxSimVal=0;
			for(j=0;j<NUMOFNODES;j++)
			{
				if(i!=j&&isSelected[j]==0)
				{
					if(NMIDataMatrix[i*NUMOFNODES+j]>maxSimVal)
					{		
						maxSimVal=NMIDataMatrix[i*NUMOFNODES+j];
						clusterIndex[i]=clusterIndex[j];
					}
				}
			}
		}	
	}
}
	
void obtainingPrototypeFeature(double *NMIDataMatrix, double *scaledcsvDataMatrix, int *isSelected, int *clusterIndex)
{
	int i,j,jj,index;
	int *clusterPrototypeFeature;
	int *numClusterEle;
	double *simDegree;
	int maxIndex;
	double maxVal;
	double maxVar;

	simDegree=(double *)calloc(NUMOFNODES,sizeof(double));	
	clusterPrototypeFeature=(int *)calloc(K,sizeof(int));
	numClusterEle=(int *)calloc(K,sizeof(int));

	for(i=0;i<NUMOFNODES;i++)
	{
		numClusterEle[clusterIndex[i]]++;
	}

        for(i=0;i<NUMOFNODES;i++)
        {
		//count=1;
                for(j=0;j<NUMOFNODES;j++)
                {
                        if(i!=j&&clusterIndex[i]==clusterIndex[j])
                                simDegree[i]+=NMIDataMatrix[i*NUMOFNODES+j];
                }
        }

        printf("\nCluster Degree value of each feature\n");
        for(i=0;i<NUMOFNODES;i++)
                printf(",%lf",simDegree[i]);
	
	printf("\nNumber of cluster elements\n");
	for(i=0;i<K;i++)
		printf(",%d",numClusterEle[i]);
	
	/*for(i=0;i<K;i++)
	{
		for(j=0;j<NUMOFNODES;j++)
			if(clusterIndex[j]==i)
				printf("\nFeature %d is in %d cluster", j+2, i);
	
	}*/

	//printing the cluster proto features
	printf("\nCluster Prototype features\n");	
	for(index=0;index<K;index++)
	{
		maxVal=0;
		maxVar=0;
		
		for(j=0;j<NUMOFNODES;j++)
		{
			if(clusterIndex[j]==index)
			{
				//if(simDegree[j]+variance(scaledcsvDataMatrix,j)>maxVal)
				if(simDegree[j]==0)
				{
                                	clusterPrototypeFeature[index]=j;
                                	break;
				}
				//1 if((simDegree[j]/numClusterEle[index]>maxVal)&&(variance(scaledcsvDataMatrix,j)>maxVar))
				//2 if((simDegree[j]+variance(scaledcsvDataMatrix,j))>maxVal)
				//3 if(((simDegree[j]/numClusterEle[index])+variance(scaledcsvDataMatrix,j))>maxVal)
				//4 
				if(vararray[j]>maxVal)
				{
					clusterPrototypeFeature[index]=j;
					//1 maxVal=simDegree[j]/numClusterEle[index];
					//1 maxVar=variance(scaledcsvDataMatrix,j);
					//2 maxVal=simDegree[j]+variance(scaledcsvDataMatrix,j);
					//3 maxVal=(simDegree[j]/numClusterEle[index])+variance(scaledcsvDataMatrix,j);
					//4 
					maxVal=vararray[j];
				}
			}
		}
		printf(",%d",clusterPrototypeFeature[index]);
	}
	
	for(index=0;index<K;index++)
	{
		printf("\nCluster with index %d has the following features\n",index);
		for(j=0;j<NUMOFNODES;j++)
		{
			if(clusterIndex[j]==index)
				printf("%d,",j);
		}

	}
	printf("\nVar is %lf", vararray[0]);
	
	for(j=0;j<NUMOFNODES;j++)
	{
		isSelected[j]=1;
		for(index=0;index<K;index++)
		{
			if(clusterPrototypeFeature[index]==j)
			{
				isSelected[j]=0;
				break;
			}	
		}
	}
}

double computeDensity (double *NMIDataMatrix, int *isSelected)
{
	int i, j;
	int numNodesInducedSet = 0;
	double sumWeights = 0;
	double density;
	
	for (i=0; i<NUMOFNODES; i++)
	{
		if(isSelected[i]==0)
			numNodesInducedSet++;
		for (j=i+1; j<NUMOFNODES; j++)
		{
			if(isSelected[i]==0 && isSelected[j]==0 && (NMIDataMatrix[i*NUMOFNODES+j]<=THRESHOLD))
			{
			//	printf("%lf", NMIDataMatrix[i*NUMOFNODES+j]);
				sumWeights += NMIDataMatrix[i*NUMOFNODES+j];
			}
		}
	}

	//printf("\nSumWeights is %lf",sumWeights);
	//printf("\nNumNodesInducedSet is %d",numNodesInducedSet);
	//version 1
	density = sumWeights/numNodesInducedSet;
	//version 2
	//density = sumWeights/comb(numNodesInducedSet,2);
	
	return density;
}

double computeInducedDegree (double *NMIDataMatrix, int *isSelected, int i)
{
	int j;
	double sumDegreeWeights = 0;
	double inducedDegree;
	int numInducedNodes=1;
	
	for (j=0; j<NUMOFNODES; j++)
	{
		if(isSelected[j]==0 && i!=j && NMIDataMatrix[i*NUMOFNODES+j]<=THRESHOLD)
		{
			sumDegreeWeights += NMIDataMatrix[i*NUMOFNODES+j];				
			numInducedNodes++;
		}
	}
	
	//version 1
	inducedDegree = sumDegreeWeights;
	//version 2
	//inducedDegree = sumDegreeWeights/(numInducedNodes-1);
	
	return inducedDegree;
}
		
void computeRanking (double *inducedDegree, int *isShortlisted, int *rank)
{
	int i, j;
	int r = 0;
	int bestIndex;
	int sizeCheckList;
	int *isToBeChecked;
	double bestValue;

	isToBeChecked = (int *)calloc(NUMOFNODES,sizeof(int));
	sizeCheckList = 0;
	
	for (i=0; i< NUMOFNODES; i++)
	{
		if(isShortlisted[i]==1)
		{
			isToBeChecked[i] = 1;
			sizeCheckList++;
		}
	}
	

	for (i=0; i< sizeCheckList; i++)
	{
		bestValue = -999999;

		for (j=0; j< NUMOFNODES; j++)
		{
			if(isToBeChecked[j]==1)	
			{
				if(inducedDegree[j] > bestValue)
				{
					//printf("\nBefore swapping rank[i] = %d and rank[j] = %d", rank[i], rank[j]);
					bestIndex = j;
					bestValue = inducedDegree[j];
					//printf("\nAfter swapping rank[i] = %d and rank[j] = %d", rank[i], rank[j]);
				}
			}
		}
		rank[bestIndex] = r++;
		isToBeChecked[bestIndex] = 0;
	}
	
	//for (i=0; i< NUMOFNODES; i++)
	//	printf("\nFeature %d has incduceed degree of %lf and its Rank is %d", i, inducedDegree[i], rank[i]);
}

/*
void computeRanking (double *inducedDegree, int *isShortlisted, int *rank)
{
	int i, j;
	int r = 0;

	for (i=0; i< NUMOFNODES; i++)
	{
		if(isShortlisted[i]==1)	
			rank[i] = r++;
		else
			rank[i] = 999999;
	}
	
	//for (i=0; i< NUMOFNODES; i++)
	//	printf("\nRank is %d", rank[i]);

	for (i=0; i< NUMOFNODES-1; i++)
	{
		for (j=i+1; j< NUMOFNODES; j++)
		{
			if(isShortlisted[i]==1 && isShortlisted[j]==1)	
			{
				if(inducedDegree[i] < inducedDegree[j])
				{
					//printf("\nBefore swapping rank[i] = %d and rank[j] = %d", rank[i], rank[j]);
					swap (&rank[i], &rank[j]);
					//printf("\nAfter swapping rank[i] = %d and rank[j] = %d", rank[i], rank[j]);
				}
			}
		}
	}
	
}

void swap( int *a, int *b)
{
	int t;
	*a = *a + *b;
	*b = *a - *b;
	*a = *a - *b;
}*/

int comb(int n,int t)
{
        int val;
        val=fact(n)/(fact(n-t)*fact(t));

        return val;
}

int fact(int n)
{
        int f=1;
        int i;
        for(i=1;i<=n;i++)
                f=f*i;
        return f;
}

double mean(double *scaledcsvDataMatrix, int j)
{
	double sum = 0;
	double meanVal;
	int i;

	for(i=0;i<NUMOFSAMPLES;i++)
	{
		sum += scaledcsvDataMatrix[i*NUMOFNODES+j];
	}

	meanVal = sum/NUMOFSAMPLES;
	return meanVal;
}

void variance(double *scaledcsvDataMatrix)
{
	double sum;
	double varVal;
	int i,j;
printf("var function");
	for (j=0;j<NUMOFNODES;j++)
        {
		sum = 0;
		for(i=0;i<NUMOFSAMPLES;i++)
		{
		//sum += powl((scaledcsvDataMatrix[i*NUMOFNODES+j]-mean(scaledcsvDataMatrix,j)),2);
		sum += ((scaledcsvDataMatrix[i*NUMOFNODES+j]-mean(scaledcsvDataMatrix,j))*(scaledcsvDataMatrix[i*NUMOFNODES+j]-mean(scaledcsvDataMatrix,j)));
		}
        
	
	//**varVal = sqrt(sum/(NUMOFSAMPLES-1));
	//varVal = sum/(NUMOFSAMPLES-1);
        	vararray[j]=sum/(NUMOFSAMPLES-1);
	}
	return 0;
}
