#include "MultinomialDataSet.h"
/*
  Methods for the MultinomialDataSet class.
 */
//CONSTRUCTOR
MultinomialDataSet::MultinomialDataSet()
{
  //HARD-WIRE THE GLOBAL HYPERPARAMETER
  globalHyperParameter=2;
}
//*****************************************************************************
//*****************************************************************************
//CONSTRUCTOR - READ IN DATA FROM FILE
MultinomialDataSet::MultinomialDataSet(string dataFile, double hyperParam)
{
  //ASSIGN THE GLOBAL HYPERPARAMETER
  globalHyperParameter=hyperParam;
  //READ IN THE DATA FROM FILE
  ReadInData(dataFile);
}
//*****************************************************************************
//*****************************************************************************
//CONSTRUCTOR - ACCEPT INPUT DATA ARRAY
//need this in order to link to R
MultinomialDataSet::MultinomialDataSet(const vector<vector<int> >& inputData, double hyperParam)
{
  //DECLARATIONS
  int i, j;
  int currentValue;
  int minVal=INT_MAX, maxVal=INT_MIN;//initial values that will be over-written in the first pass
  //ASSIGN THE GLOBAL HYPERPARAMETER
  // JM: TODO: check to make sure private values globalHyperParameter, data,
  // hyperParameter are already declared and can be accessed like this
  globalHyperParameter=hyperParam;
  //COPY THE DATA INTO THE OBJECT
  data = inputData;
  //FIND THE DATA SIZE
  nDataItems = data.size();
  nFeatures  = data[0].size();
  //FIND THE MIN, MAX DATA VALUES
  for (i=0; i<nDataItems; i++)
  {
    for (j=0; j<nFeatures; j++)
    {
      currentValue = data[i][j];
      //keep track of the min and max values
      if (currentValue<minVal) minVal=currentValue;
      if (currentValue>maxVal) maxVal=currentValue;
    }
  }
  //FIND THE NUMBER OF DISCRETE DATA VALUES AND ADJUST THEM TO RUN FROM 0:(n-1)
  nFeatureValues = maxVal - minVal + 1;
  for (i=0; i<nDataItems; i++)
    for (j=0; j<nFeatures; j++)
    {
      data[i][j] -= minVal;//set minimum value to zero
    }
  //DEFINE HYPERPARAMETERS
  ComputeHyperParameters();
}
//*****************************************************************************
//*****************************************************************************
//READ DATA IN FROM A FILE
void MultinomialDataSet::ReadInData(string dataFile)
{
  //DECLARATIONS
  int                  i, j;
  int                  inputValue;
  int                  minVal=INT_MAX;//initial values that will be over-written in the first pass
  int                  maxVal=INT_MIN;//initial values that will be over-written in the first pass
  fstream              file;
  //FIND THE SIZE OF THE INPUT DATA
  FindDataSize(dataFile);
  //OPEN DATA FILE (WITH ERROR-CHECKING)
  file.open (dataFile.c_str(), fstream::in);//dataFile is a string object, so need to convert to a char array here!
  if(file.fail())
  {
    cout<<"Failed to open file "<<dataFile<<"."<<endl;
    if(!system("PAUSE"))/* didn't pause */;
  }
  //READ THE DATA INTO OUR DATA VECTOR
  for (i=0; i<nDataItems; i++)
  {
    data.push_back( vector<int>(nFeatures, 0));//add a vector for this item
    for (j=0; j<nFeatures; j++)
    {
      file >> inputValue;
      data[i][j] = inputValue;
      //also keep track of the min and max values
      if (inputValue<minVal) minVal=inputValue;
      if (inputValue>maxVal) maxVal=inputValue;
    }
  }
  file.close();
  //FIND THE NUMBER OF DISCRETE DATA VALUES AND ADJUST THEM TO RUN FROM 0:(n-1)
  nFeatureValues = maxVal - minVal + 1;
  for (i=0; i<nDataItems; i++)
    for (j=0; j<nFeatures; j++)
    {
      data[i][j] -= minVal;//set minimum value to zero
    }
  //DEFINE HYPERPARAMETERS
  ComputeHyperParameters();
}
//*****************************************************************************
//*****************************************************************************
//COMPUTE THE HYPERPARAMETERS FOR THIS MULTINOMIAL MODEL
//Katherine Heller's BHC code calcualte hyperparameters as being proportional to the population mean proportion of each featuer value
//this is then scaled by the global hyperparameter (cc in her code; here??)
//we want to set this global hyperparameter using Evidence optimisation...
//(by-hand for now.  Automatic at a later date...)
void MultinomialDataSet::ComputeHyperParameters()
{
  //DECLARATIONS
  int            i, j, k;
  vector<double> dataCounter;
  vector<double> newHyperParameter;
  //COMPUTE THE HYPERPARAMETER VALUES
  for (i=0; i<nFeatures; i++)
  {
    //initialise the data counter, vector of new hyperparameters
    dataCounter       = vector<double>(nFeatureValues, 1);
    newHyperParameter = vector<double>(nFeatureValues, 0);
    //count the number of each occurrences of each data value
    for (j=0; j<nDataItems; j++)
      dataCounter[data[j][i]]++;
    //calculate the hyperparameter values
    for (k=0; k<nFeatureValues; k++)
      newHyperParameter[k] = globalHyperParameter * dataCounter[k] / (nDataItems+1); //this is the formula used by Katherine
    //store the hyperparameters
    hyperParameter.push_back(newHyperParameter);
  }
}
//*****************************************************************************
//*****************************************************************************
// COMPUTE THE LOG-EVIDENCE FOR A SINGLE CLUSTER CONTAINING THE DATA ITEMS IDENTIFIED BY 'itemIndex'
double MultinomialDataSet::SingleClusterLogEvidence(vector<int>& itemIndex)
{
  //DECLARATIONS
  int                   i, j;
  double                s1=0, s2=0, s3=0, s4=0;
  int                   index, dataValue;
  int                   nCurrentItems=itemIndex.size();
  int                   counter=0;
  vector< vector<int> > dataCounter;
  vector<double>        sumHyper=vector<double>(nFeatures, 0);
  //INITIALISE THE DATA COUNTER
  for (i=0; i<nFeatures; i++)
    dataCounter.push_back(vector<int>(nFeatureValues, 0));
  //COUNT THE INSTANCES OF EACH VALUE IN EACH FEATURE
  for (i=0; i<nFeatures; i++)
  {
    for (j=0; j<nCurrentItems; j++)
    {
      index     = itemIndex[j];
      dataValue = data[index][i];
      dataCounter[i][dataValue]++;//increment the counter for the current feature and feature-value
      counter++;
    }
  }
  //CONSTRUCT USEFUL INTERMEDIATE VALUES
  //sum of hyperparameter values for each feature
  for (i=0; i<nFeatures; i++)
    for (j=0; j<nFeatureValues; j++) sumHyper[i] += hyperParameter[i][j];
  //intermediates to calculate the log-ev
  for (i=0; i<nFeatures; i++)
  {
    for (j=0; j<nFeatureValues; j++)
    {
      s1 += gammaln(hyperParameter[i][j]);
      s2 += gammaln(hyperParameter[i][j] + dataCounter[i][j]);
    }
    s3 += gammaln(sumHyper[i]);
    s4 += gammaln(sumHyper[i] + nCurrentItems);
  }
  //CALCULATE AND RETURN THE LOG-EVIDENCE
  return(-s1 +s2 +s3 -s4);
}

double MultinomialDataSet::GetClusterNoise(int NodeID)
{ return 0.0; } // No noise
