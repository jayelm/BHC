/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include <iostream>
#include <math.h>

#include "multinomial_header.h"

using namespace std;

/* ----------------------------------------------------------------------
   Inputs:
      tr_node - pointers of nodes in the dendrogram.
      dim - feature size or dimension of the input data.
      obs - total number of input data points.
      cc - scalar value of precision of Beta prior (initialized in main.cpp).
      alp - scalar value of Dirichlet Process hyperparameter (initialized in main.cpp).
      min_wt - initial weight value (default=-1000.0 defined in main.cpp).
   
   Output:
      logEvidence - updates NODE information for construction of the dendrogram.
---------------------------------------------------------------------- */

double bayeslink_binf(NODE* tr_node,
              // JM: Number of features (columns)
		      int dim,
              // JM: Number of observations (rows)
		      int obs,
              // JM: Global hyperparameter found by FindOptimalHyperparameter
		      double cc,
              // JM: Hardwired to 0.001
		      double alp,
              // JM: -Inf
		      double min_wt,
              // JM: Number of discrete values used in multinomial
		      int nFeatureValues)
{
  //----------------------------------------------------------------------
  // DECLARATIONS --------------------------------------------------------
  //----------------------------------------------------------------------
  int      i, j;
  // JM: Will be dynamically allocated 2D array of hyperparameters
  double** hyperParameters;
  double   tr1,tr2,a,b,ckt,pk,gell, wt_temp;
  int      node1=-1, node2=-1, merged_node, itr, k;
  double   logEvidence=0.;
  //----------------------------------------------------------------------
  // DEFINE HYPERPARAMETERS ----------------------------------------------
  //----------------------------------------------------------------------
  // JM: These hyperparameters might be important...might want to return these later
  hyperParameters = CalculateHyperparameters(tr_node, dim, obs, nFeatureValues, cc);

  //----------------------------------------------------------------------
  // COMPUTE LOG-EVIDENCE FOR SINGLE DATA POINTS -------------------------
  //----------------------------------------------------------------------
  for (i=0;i<obs;i++){
    // JM: wt gets assigned later as "mergeWeight" in bhc.cpp; it's the ratio
    // JM: of the two hypotheses. Of course with only one node the ratio is
    // JM: just the binevidence
    tr_node[i].wt[i] = binevidence(tr_node,dim,hyperParameters,i,-1, nFeatureValues);
    tr_node[i].ck    = log(alp);
	// JM: nk is the number of values associated with the node!
    tr_node[i].nk    = 1.0;
    // JM: Like I said above, the ratio is just the binevidence.
    // JM: So den is p(D_k \mid H_1)
    tr_node[i].den   = tr_node[i].wt[i];
  }

  //----------------------------------------------------------------------
  // COMPUTE LOG-EVIDENCE OF EVERY PAIR OF POINTS ------------------------
  //----------------------------------------------------------------------
  // JM: This is O(n^2) part
  for (i=0;i<obs;i++)
    for (j=i+1;j<obs;j++){
      tr1                = log(alp) + fast_gammaln(tr_node[i].nk+tr_node[j].nk);
      tr2                = tr_node[i].ck + tr_node[j].ck;
      a                  = max(tr1,tr2);
      b                  = min(tr1,tr2);
      ckt                = a + log(1+exp(b-a));
      // JM: pk is \pi_k. Appears to be set from `alp` (dirichlet process hyperparameter)
      // JM: and some other odd stuff
      pk                 = log(alp) + fast_gammaln(tr_node[i].nk+tr_node[j].nk) - ckt;
      // JM: Appears to be the other component of numerator of r_k eq in Figure 2
      // JM: (binevidence!!)
      gell               = binevidence(tr_node,dim,hyperParameters,i,j, nFeatureValues);
      // JM: num1 has something to do with log P(D_k \mid H_1); product of
      // JM: \pi_k and p(D_k \mid H_1)
      tr_node[i].num1[j] = pk + gell;
      // JM: num2 has something to do with log P(D_k \mid T_k)
      tr_node[i].num2[j] = tr2 - ckt + tr_node[i].wt[i] + tr_node[j].wt[j];
      // JM: So wt is log (r_k / 1 - r_k) ?? maybe?; the log-odds ratio of num1/num2
      tr_node[i].wt[j]   = tr_node[i].num1[j] - tr_node[i].num2[j];
    }
  
  //----------------------------------------------------------------------
  // MERGE AND FORM HIERARCHICAL CLUSTERS --------------------------------
  //----------------------------------------------------------------------
  merged_node = obs-1;
  for (itr=1;itr<obs;itr++){
    //Find clusters to merge.
    wt_temp = min_wt;
    merged_node++;
    for (i=0;i<2*obs;i++){
      if (tr_node[i].flag==0)
	for (j=i+1;j<2*obs;j++)
	  if ((tr_node[j].flag==0)&&(tr_node[i].wt[j]>wt_temp)){
	    wt_temp = tr_node[i].wt[j];
	    node1   = i;
	    node2   = j;
	  }
    }
    //----------------------------------------------------------------------
    // UPDATE NODE INFORMATION ---------------------------------------------
    //----------------------------------------------------------------------
    tr_node[merged_node].pleft           = node1;
    tr_node[merged_node].pright          = node2;
    tr_node[merged_node].wt[merged_node] = tr_node[node1].wt[node2];
    a = max(tr_node[node1].num1[node2],tr_node[node1].num2[node2]);
    b = min(tr_node[node1].num1[node2],tr_node[node1].num2[node2]);
    tr_node[merged_node].den = a + log(1+exp(b-a));
    tr1 = log(alp) + fast_gammaln(tr_node[node1].nk+tr_node[node2].nk);
    tr2 = tr_node[node1].ck + tr_node[node2].ck;
    a   = max(tr1,tr2);
    b   = min(tr1,tr2);
    tr_node[merged_node].ck = a + log(1+exp(b-a));
    tr_node[merged_node].nk = tr_node[node1].nk+tr_node[node2].nk;
      
    //----------------------------------------------------------------------
    // GET NEW DATA POINT --------------------------------------------------
    //----------------------------------------------------------------------
    for (k=0;k<dim;k++) 
      for (j=0;j<nFeatureValues;j++)
	tr_node[merged_node].dat[k][j] = tr_node[node1].dat[k][j] + tr_node[node2].dat[k][j];
    tr_node[merged_node].vec_no = tr_node[node1].vec_no + tr_node[node2].vec_no;
    for (k=0;k<merged_node;k++){
      if ((tr_node[k].flag!=1)&&(k!=node1)&&(k!=node2)){
	tr1  = log(alp) + fast_gammaln(tr_node[merged_node].nk+tr_node[k].nk);
	tr2  = tr_node[merged_node].ck + tr_node[k].ck;
	a    = max(tr1,tr2);
	b    = min(tr1,tr2);
	ckt  = a + log(1+exp(b-a));
	pk   = log(alp) + fast_gammaln(tr_node[merged_node].nk+tr_node[k].nk) - ckt;
	gell = binevidence(tr_node,dim,hyperParameters,merged_node,k, nFeatureValues);
	tr_node[k].num1[merged_node] = pk + gell;
	tr_node[k].num2[merged_node] = tr2 - ckt + tr_node[merged_node].den + tr_node[k].den;
	tr_node[k].wt[merged_node]   = tr_node[k].num1[merged_node] - tr_node[k].num2[merged_node];
      }
    }
    tr_node[node1].flag = 1;
    tr_node[node2].flag = 1;
  }//End merging loop
  
  //----------------------------------------------------------------------
  // FREE ALLOCATED MEMORY -----------------------------------------------
  //--- -------------------------------------------------------------------
  for (i=0; i<nFeatureValues; i++)
    delete[] hyperParameters[i];
  delete[] hyperParameters;
  
  //----------------------------------------------------------------------
  // RETURN THE GLOBAL LOG-EVIDENCE ESTIMATE -----------------------------
  //----------------------------------------------------------------------
  //this should be the final (i.e. root) node; 'den' is the marginal likelihood bound
  logEvidence = tr_node[merged_node].den;

  return logEvidence;
}
