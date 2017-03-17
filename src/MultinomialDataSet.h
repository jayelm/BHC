#include <iostream>
#include <string>
#include <vector>
#include "DataSet.h"
/*
  This class inherits the superclass DataSet.
  We do this so that BHC can run with different types of data.
  (this is like using a Java interface)
*/
class MultinomialDataSet : public DataSet
{
public:
  //CONSTRUCTORS
  MultinomialDataSet();
  MultinomialDataSet(string dataFile, double hyperParam);
  MultinomialDataSet(const vector<vector<int> >& inputData, double hyperParam);
  //OTHER METHODS
  void   ComputeHyperParameters();
  double SingleClusterLogEvidence(vector<int>& itemIndex);  // returns the log-evidence for a cluster that contains data items defined by itemIndex
  void   ReadInData(string dataFile);                      //read data in from a file
  double GetClusterNoise(int NodeID);
  
  //TAGS
  int    nFeatureValues;                                   //the number of allowed values (1 to n) for the multinomial variables

private:
  double                   globalHyperParameter;           //Not quite sure about this.  Check!!
  vector< vector<int> >    data;                           //array to store the multinomial data
  vector< vector<double> > hyperParameter;                 //array to store hyperparameters
};
