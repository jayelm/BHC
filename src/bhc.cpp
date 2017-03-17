/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <limits>

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
  void bhcWrapper(int* dataTypeID_input, double* inputData,
		  int* nDataItems_input, int* nFeatures_input, double* ghpInput,
          int* num_reps,
          int* nFeatureValues_input, double* logEvidence,
		  int* node1, int* node2, int* mergeOrder, double* mergeWeight,
		  int* numThreads, int* m)
  {
    // Declarations
    // cout << "asdfoijasdfoijasdfoijasdf\noiajsdfoijasdfoiajsdf\noijasdfoijsdf" << endl;
    int i, j, counter1=0;
    string dataFile, dataType, outputFile;
    vector<Node> treeNode;
    MultinomialDataSet* dataSet=NULL;
    DirichletProcessMixture bhc;
    // cout << "Step2" << endl;;
    
    // Since our inputs are pointers, this saves some effort later
    int nDataItems=*nDataItems_input, 
      nFeatures=*nFeatures_input,
      dataTypeID=*dataTypeID_input;
    double globalHyperParameter= *ghpInput; // unused
    vector<double> noise, timePoints_copy;
    // For ease, we define both these 2D arrays; we only need one
    vector<vector<int> > data_int;
    vector<vector<double> > data_double;
    // cout << "Step3" << endl;;

    // Parallelise if required
#ifdef SUPPORT_OPENMP
    omp_set_num_threads(MAX(*numThreads,1));
    cout << "SUPPORT_OPENMP ENABLED; set num threads!!" << endl;
#endif
    
    if (dataTypeID == 0) {
        dataType = "multinomial";
    } else {
      cout << "Oops! invalid dataType! tag one" << endl;
      return;
    }

    // cout << "Step4" << endl;;
    // Copy data into an array of doubles
    if (dataType == "multinomial") {
        // Copy data into an array of ints
        for (i = 0; i < nDataItems; i++) {
            data_int.push_back(vector<int>(nFeatures, 0));
            for (j = 0; j < nFeatures; j++) {
                data_int[i][j] = (int) inputData[counter1++];
            }
        }
    } else {
        // Use doubles
        cout << "ERROR: Disabling time course for now" << endl;
        return;
    }

    // Instantiate the required type of dataset object
    // cout << "Step5" << endl;;
    if (dataType=="multinomial") {
      dataSet = new MultinomialDataSet(data_int, globalHyperParameter);
    } else {
      cout << "Oops! Invalid dataType! tag two" << endl;
      return;
    }

    // Run clustering analysis
    // cout << "Step6 (Before greedyc lustering)" << endl;;
    if(*m<2 || *m>=nDataItems) {
      treeNode = bhc.GreedyClustering(*dataSet, false);
    }
    else {
      treeNode = bhc.RandomisedClustering(*dataSet, *m, 0, false);
    }

    // Pass the dendrogram data back to R
    // cout << "Step7 (After cluster; pass dendro)" << endl;;
    for (i=0; i<nDataItems-1; i++) {
      node1[i] = treeNode[i+nDataItems].GetLeftChildNodeID()  + 1;
      node2[i] = treeNode[i+nDataItems].GetRightChildNodeID() + 1;
      mergeOrder[i] = i + 1;
      mergeWeight[i] = treeNode[i+nDataItems].GetClusterLogEvidence();
    }
    *logEvidence = treeNode.back().GetLowerBoundLogEvidence();

    // Cleanup
    // cout << "Step8 (Done, cleanup)" << endl;;
    delete dataSet;
  }
}
