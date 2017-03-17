/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#ifndef DATASET_H
#define DATASET_H

#include "header.h"

/* ----------------------------------------------------------------------
   This is an abstract base class which allows BHC to handle different
   types of data. Note that since the virtual functions have not been
   implemented here, then DataSet cannot be instantiated directly.
---------------------------------------------------------------------- */

class DataSet
{
 public:
  virtual ~DataSet() {}
  // JM NOTE: Multinomial takes a straight-up vector...rather than a pointer
  virtual double SingleClusterLogEvidence(vector<int>& itemIndex) = 0;
  virtual void ReadInData(string dataFile) = 0;
  virtual double GetClusterNoise(int nodeID) = 0;
  // This seems useful...why isn't this implemented?
  // virtual vector<double> GetDataForCluster(vector<int> itemIndex) = 0;

  // Implemented in DataSet.cpp
  void FindDataSize(string dataFile);
  void ReadInDataVector(vector<int> inputData, int nDataItems, int nFeatures);
  int Get_nDataItems();
  int Get_nFeatures();
  
  
  int nDataItems; // the number of data items stored by DataSet
  int nFeatures; // the number of features possessed by each data item
  int noise_mode; // indicates whether there is precalculated fixed noise
  int robust_mode;
  int reps;//the number of replicates per observation
  string dataType;
  vector<double>  noiseData; // noise values (computed from estimators)
};
#endif

