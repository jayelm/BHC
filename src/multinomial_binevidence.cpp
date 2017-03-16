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
   This calculates the log evidence for data under the Bernoulli model.
   A different routine should be written in place of this for multinomial data.
   
   Inputs:
      tr_node - pointers of nodes in the dendrogram.
      dim - feature size or dimension of the input data.
      alpha - vector (length=dim) of hyperparameters for Beta prior.
      beta - vector (length=dim) of hyperparameters for Beta prior.
      index1 - index of the first input cluster.
      index2 - index of the second input cluster (set to -1 if there is only one input).
   
   Outputs:
      logEvidence - scalar log evidence of data under multivariate Bernoulli model.
---------------------------------------------------------------------- */

double binevidence(NODE* tr_node,
		   int dim,
                   // JM: The number of values X number of features array of
                   // JM: hyperparams
		   double** hyperParameters,
		   int index1,
                   // JM: When used from computing logE of single points,
                   // JM: always -1
		   int index2,
		   int nFeatureValues)
{
  //----------------------------------------------------------------------
  // DECLARATIONS --------------------------------------------------------
  //----------------------------------------------------------------------
  double   logEvidence;
  double*  sumHyper = new double[dim];
  int      i, j, k;
  int      nTotalData;
  double   s2=0., s3=0., s5=0., s6=0.; // unused: s1 and s4
  //construct dynamically a 2D array 
  double** dataCounter;
  dataCounter = new double* [dim];
  for (int i=0; i<dim; i++) dataCounter[i] = new double [nFeatureValues];
  //----------------------------------------------------------------------
  // AGGREGATE VALUES THAT DEPEND ON WHETHER WE'RE MERGING TWO CLUSTERS --
  //----------------------------------------------------------------------
  nTotalData = tr_node[index1].vec_no;  //the total number of data items here
  //count the number of each data value
  for (i=0;i<dim;i++)
    for (k=0;k<nFeatureValues;k++)
      dataCounter[i][k] = tr_node[index1].dat[i][k];
  //if we're merging two clusters., need to add the info for the second cluster
  if (index2!=-1){
    nTotalData += tr_node[index2].vec_no;
    //count the number of each data value
    for (i=0;i<dim;i++)
      for (k=0;k<nFeatureValues;k++)
	dataCounter[i][k] += tr_node[index2].dat[i][k];
  }
  //----------------------------------------------------------------------
  // CONSTRUCT INTERMEDIATE VALUES TO USE IN FINDING LOG-EVIDENCE --------
  //----------------------------------------------------------------------
  //sum of hyperparameter values for each feature
  for (i=0;i<dim;i++){
    sumHyper[i] = 0.;
    for (j=0;j<nFeatureValues;j++) sumHyper[i] += hyperParameters[j][i];
  }
  //intermediates to calculate the log-ev
  for (i=0;i<dim;i++){
    //some require sums over the different allowed feature values
    for (j=0;j<nFeatureValues;j++){
      //if (dataCounter[i][j]!=0) 
      //s1 += fast_gammaln(dataCounter[i][j]); //use lngamma to find log-factorial, except for zero
      s2 += fast_gammaln(hyperParameters[j][i]);
      s3 += fast_gammaln(hyperParameters[j][i] + dataCounter[i][j]);
    }
    //s4 += fast_gammaln(nTotalData); //really only need to do this once!
    s5 += fast_gammaln(sumHyper[i]);
    s6 += fast_gammaln(sumHyper[i] + nTotalData);
  }
  //----------------------------------------------------------------------
  // CALCULATE THE LOG-EVIDENCE ------------------------------------------
  //----------------------------------------------------------------------
  logEvidence =  -s2 +s3 +s5 -s6; 
  //----------------------------------------------------------------------
  // FREE ALLOCATED MEMORY AND RETURN LOG-EVIDENCE -----------------------
  //----------------------------------------------------------------------
  for (i=0; i<dim; i++)
    delete[] dataCounter[i];
  delete []dataCounter;
  delete []sumHyper;
  return logEvidence;
}
