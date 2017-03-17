/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <limits>

#include "multinomial_header.h" // copied from 1.1.0 package

#include "header.h"
#include "MultinomialDataSet.h"
#include "Node.h"
#include "DirichletProcessMixture.h"

/* ----------------------------------------------------------------------
   This is the interface to the time-series and cubic-spline implementations.
   For the multinomial case, see bhcWrapper_multinomial(), below.
   This version supports both the greedy and randomised algorithms; cf. m.
   
   DataTypeID_input should be interpreted as follows:
      0 - multinomial
      1 - time-series (using squared-exponential covariance)
      2 - cubicspline (using cubic spline covariance)
  
   The inputData is ordered as follows:
      item1_feature1  item1_feature2 ... item2_feature1 item2_feature2 ...
   the DataSet classes expect a 2D array of data, (nItems * nFeatures).

   m is the dimension of the subset D_m used in the randomised
   algorithm. If 2<=m<nDataItems, then the randomised algorithm is used,
   else the greedy algorithm.
---------------------------------------------------------------------- */

extern "C" { // required for linking properly to R
  void bhcWrapper(int* dataTypeID_input, double* inputData, double* timePoints,
		  int* nDataItems_input, int* nFeatures_input, double* ghpInput,
		  double* noise_input, int* num_reps, int* set_noise_input,
		  int* robust_input, int* nFeatureValues_input, double* logEvidence,
		  int* node1, int* node2, int* mergeOrder, double* mergeWeight,
		  int* numThreads, int* m)
  {
    // Declarations
    int i, j, counter1=0;
    string dataFile, dataType, outputFile;
    vector<Node> treeNode;
    MultinomialDataSet* dataSet=NULL;
    DirichletProcessMixture bhc;
    
    // Since our inputs are pointers, this saves some effort later
    int nDataItems=*nDataItems_input, 
      nFeatures=*nFeatures_input,
      dataTypeID=*dataTypeID_input,
      num_reps_copy = *num_reps,
      set_noise_copy = *set_noise_input,
      robust_copy = *robust_input;
    double globalHyperParameter= *ghpInput; // unused
    vector<double> noise, timePoints_copy;
    // For ease, we define both these 2D arrays; we only need one
    vector<vector<int> > data_int;
    vector<vector<double> > data_double;

    // Parallelise if required
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(MAX(*numThreads,1));
#endif
    
    if      (dataTypeID==0) dataType = "multinomial";
    else if (dataTypeID==1) { cout << "ERROR: Disabling time course for now" << endl; return; } // dataType = "time-course";
    else if (dataTypeID==2) { cout << "ERROR: Disabling time course for now" << endl; return; }// dataType = "cubicspline";
    else
    {
      cout<<"Oops! invalid dataType! tag one"<<endl;
      return;
    }

    // Copy data into an array of doubles
    if (dataType == "multinomial") {
        // Copy data into an array of ints
        for (i = 0; i < nDataItems; i++) {
            data_int.push_back(vector<int>(nFeatures, 0));
            for (j = 0; j < nFeatures; j++) {
                data_double[i][j] = (int) inputData[counter1++];
            }
        }
    } else {
        // Use doubles
      cout<<"ERROR: Disabling time course for now"<<endl;
      return;
        // for (i=0; i<nDataItems; i++)
          // {
            // data_double.push_back(vector<double>(nFeatures, 0));
            // for (j=0; j<nFeatures; j++)
          // {
            // data_double[i][j] = inputData[counter1++];
          // }
          // }
    }

    //Read in the timePoints
    // for (i=0; i<nFeatures; i++)
      // {
        // timePoints_copy.push_back(timePoints[i]);
      // }

    // Copying over the noise vector
    // if (set_noise_copy ==2)
    // {
      // for (i=0; i<nDataItems; i++)
      // {
        // noise.push_back(noise_input[i]);
      // }
    // }
    // else if (set_noise_copy ==1)
    // {
      // noise.push_back(noise_input[0]);
    // }
    // else if (set_noise_copy ==0)
    // {
      // noise.push_back(0.0);
    // }

    // Instantiate the required type of dataset object
    if (dataType=="time-course")
    {
      cout<<"ERROR: Disabling time course for now"<<endl;
      return;
      // if (robust_copy == 0)
      // {
        // dataSet = new SquaredExponentialTimecourseDataSet(data_double);
        // dataSet->ReadInNoise(noise);
      // }
      // else if (robust_copy == 1)
      // {
        // dataSet = new RobustSquaredExponentialTimecourseDataSet(data_double);
        // dataSet->ReadInNoise(noise);
      // }
    }
    else if (dataType=="cubicspline")
    {
      cout<<"ERROR: Disabling time course for now"<<endl;
      return;
      // if (robust_copy == 0)
      // {
        // dataSet = new CubicSplineTimecourseDataSet(data_double);
        // dataSet->ReadInNoise(noise);
      // }
      // if (robust_copy == 1)
      // {
        // dataSet = new RobustCubicSplineTimecourseDataSet(data_double);
        // dataSet->ReadInNoise(noise);
      // }
    } else if (dataType=="multinomial")
    {
      dataSet = new MultinomialDataSet(data_int, globalHyperParameter);
      // No noise here
    }
    else
    {
      cout<<"Oops! Invalid dataType! tag two"<<endl;
      return;
    }

    // dataSet->SetNoiseMode(set_noise_copy);
    // if (set_noise_copy ==2)
    // {
      // dataSet->SetReps(num_reps_copy);
    // }
    // dataSet->SetRobustMode(robust_copy);
    // dataSet->SetDataType(dataType);
    // dataSet->ReadInTimePoints(timePoints_copy);

    // Run clustering analysis
    if(*m<2 || *m>=nDataItems) {
      treeNode = bhc.GreedyClustering(*dataSet, false);
    }
    else {
      treeNode = bhc.RandomisedClustering(*dataSet, *m, 0, false);
    }

    // Pass the dendrogram data back to R
    for (i=0; i<nDataItems-1; i++)
    {
      node1[i] = treeNode[i+nDataItems].GetLeftChildNodeID()  + 1;
      node2[i] = treeNode[i+nDataItems].GetRightChildNodeID() + 1;
      mergeOrder[i] = i + 1;
      mergeWeight[i] = treeNode[i+nDataItems].GetClusterLogEvidence();
    }
    *logEvidence = treeNode.back().GetLowerBoundLogEvidence();

    // Cleanup
    delete dataSet;
  }
}
