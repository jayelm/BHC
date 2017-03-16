/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include "multinomial_header.h"

/* ----------------------------------------------------------------------
   Routine to calculate the alpha hyperparameters for the BHC clustering.

   Katherine Heller's BHC code calculates hyperparameters as being
   proportional to the population mean proportion of each feature value
   this is then scaled by the global hyperparameter (cc in her code)
   we want to set this global hyperparameter using evidence optimisation
   (by hand for now; automatic at a later date)
------------------------------------------------------------------------- */

double** CalculateHyperparameters(NODE* tr_node,
				  int dim,
				  int obs,
				  int nFeatureValues,
				  double cc)
{
  //----------------------------------------------------------------------
  // DECLARATIONS --------------------------------------------------------
  //----------------------------------------------------------------------
  int i, j, k;
  //construct dynamically a 2D array 
  double** dataCounter;
  dataCounter = new double* [nFeatureValues];
  for (int i=0; i<nFeatureValues; i++) dataCounter[i] = new double [dim];
  //construct dynamically a 2D array to hold the hyperparameters
  // JM: Rows of the hyperparameter 2D array: number of discrete values used in
  // JM: multinomial. In the binary case, this is 2
  double** hyperParameters;
  hyperParameters = new double* [nFeatureValues];
  // JM: Columns of the hyperparameter 2D array: the number of features in the
  // JM: dataset
  for (int i=0; i<nFeatureValues; i++) hyperParameters[i] = new double [dim];
  //----------------------------------------------------------------------
  // INITIALISE THE DATA COUNTER -----------------------------------------
  //----------------------------------------------------------------------
  //initialise the dataCounter to 1 (i.e. add one to the overall count)
  //this could be viewed as asserting the prior knowledge that each data value is possible. 
  for (i=0;i<dim;i++) 
    for (k=0;k<nFeatureValues;k++) dataCounter[k][i] = 1;
  //----------------------------------------------------------------------
  // CALCULATE THE HYPERPARAMETER VALUES ---------------------------------
  //----------------------------------------------------------------------
  for (i=0;i<dim;i++) {
    //count the number of each data value
    for (j=0;j<obs;j++)
      for (k=0;k<nFeatureValues;k++)
	dataCounter[k][i] += tr_node[j].dat[i][k];
    //calculate the hyperparameter values
    for (k=0;k<nFeatureValues;k++)
      hyperParameters[k][i] = cc * dataCounter[k][i]  / (obs+1); //this is the formula used by Katherine
  }
  //----------------------------------------------------------------------
  // FREE MEMORY ALLOCATION ----------------------------------------------
  //----------------------------------------------------------------------
  for (int i=0; i<nFeatureValues; i++) 
    delete[] dataCounter[i];
  delete []dataCounter;
  //----------------------------------------------------------------------
  // RETURN 2D ARRAY OF HYPERPARAMETERS ----------------------------------
  //----------------------------------------------------------------------
  return hyperParameters;
}
