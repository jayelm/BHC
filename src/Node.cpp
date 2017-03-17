/* ----------------------------------------------------------------------
   BHC - Bayesian Hierarchical Clustering
   http://www.bioconductor.org/packages/release/bioc/html/BHC.html
   
   Author: Richard Savage, r.s.savage@warwick.ac.uk
   Contributors: Emma Cooke, Robert Darkins, Yang Xu
   
   This software is distributed under the GNU General Public License.
   
   See the README file.
------------------------------------------------------------------------- */

#include "Node.h"
#include <limits>

/* ----------------------------------------------------------------------
   Hard-coded to the value in Heller/Xu's code
------------------------------------------------------------------------- */

const double Node::dirichletProcessParameter=0.001;

/* ----------------------------------------------------------------------
    Default constructor.
------------------------------------------------------------------------- */

Node::Node(){}

/* ----------------------------------------------------------------------
   Factory method to instantiate this class. It creates a node consisting
   of a single data item (as specified in the argument).
------------------------------------------------------------------------- */

Node Node::CreateDataNode(DataSet& dataSet, const int arg_nodeID)
{
  Node thisNode=Node();
  double lengthScale, noiseFreeScale, noiseSigma, mixtureComponent;

  // Set all the variables to their default values
  thisNode.nodeID=arg_nodeID;
  thisNode.dataID=arg_nodeID;
  thisNode.childDataIDs.push_back(thisNode.dataID);
  thisNode.childNodeIDs.push_back(thisNode.nodeID);
  thisNode.leftChildNodeID=thisNode.rightChildNodeID=-1;
  thisNode.log_d_k=log(dirichletProcessParameter);
  thisNode.mergeFlag=false;
  thisNode.mergePrior=1.0;
  thisNode.clusterLogEvidence=-numeric_limits<double>::infinity();

  // Find the optimised hyperparameters for this node and compute the
  // overall log-evidence estimate
  thisNode.lowerBoundLogEvidence
    = dataSet.SingleClusterLogEvidence(thisNode.childDataIDs);
  return thisNode;
}

/* ----------------------------------------------------------------------
   Factory method to instantiate this class. It takes two nodes and
   merges them.
------------------------------------------------------------------------- */

Node Node::CreateMergerNode(DataSet& dataSet,
			    const Node& node1,
			    const Node& node2,
			    const int arg_nodeID)
{
  Node thisNode=Node();
  double tr1,tr2,a,b,ckt,pk,gell,num1,num2;
  double lengthScale, noiseFreeScale, noiseSigma, mixtureComponent;

  assert(node1.AllowedToMerge() && node2.AllowedToMerge());
  
  // Default values
  thisNode.nodeID=arg_nodeID;
  thisNode.dataID=-1;
  thisNode.childNodeIDs=node1.childNodeIDs;
  thisNode.childNodeIDs.insert(thisNode.childNodeIDs.end(),
			       node2.childNodeIDs.begin(),
			       node2.childNodeIDs.end());
  thisNode.childDataIDs=node1.childDataIDs;
  thisNode.childDataIDs.insert(thisNode.childDataIDs.end(),
			       node2.childDataIDs.begin(),
			       node2.childDataIDs.end());
  thisNode.leftChildNodeID=node1.GetNodeID();
  thisNode.rightChildNodeID=node2.GetNodeID();
  thisNode.mergeFlag=false;

  // Compute the cluster log evidence and lower bound log evidence
  // and update the log_d_k and merge prior
  tr1=log(dirichletProcessParameter)+gammaln(node1.mergePrior+node2.mergePrior);
  tr2=node1.log_d_k+node2.log_d_k;
  a=max(tr1,tr2);
  b=min(tr1,tr2);
  ckt=a+log(1.0+exp(b-a));
  pk=tr1-ckt;
  gell=dataSet.SingleClusterLogEvidence(thisNode.childDataIDs);
  num1=pk+gell;
  num2=tr2-ckt+node1.lowerBoundLogEvidence+node2.lowerBoundLogEvidence;
  thisNode.clusterLogEvidence=num1-num2;
  a=max(num1,num2);
  b=min(num1,num2);
  thisNode.lowerBoundLogEvidence=a+log(1.0+exp(b-a));
  thisNode.log_d_k=ckt;
  thisNode.mergePrior=node1.mergePrior+node2.mergePrior;

  return thisNode;
}

/* ---------------------------------------------------------------------- */

int Node::GetDataID() const
{
  return dataID;
}

/* ---------------------------------------------------------------------- */

int Node::GetNodeID() const
{
  return nodeID;
}

/* ---------------------------------------------------------------------- */

int Node::GetLeftChildNodeID() const
{
  return leftChildNodeID;
}

/* ---------------------------------------------------------------------- */

int Node::GetRightChildNodeID() const
{
  return rightChildNodeID;
}

/* ---------------------------------------------------------------------- */

double Node::GetClusterLogEvidence() const
{
  return clusterLogEvidence;
}

/* ---------------------------------------------------------------------- */

double Node::GetLowerBoundLogEvidence() const
{
  return lowerBoundLogEvidence;
}

/* ---------------------------------------------------------------------- */

int Node::GetNumItems() const
{
  return childNodeIDs.size();
}

/* ---------------------------------------------------------------------- */

bool Node::AllowedToMerge() const
{
  return !mergeFlag;
}

/* ---------------------------------------------------------------------- */

vector<int> Node::GetChildDataIDs() const
{
  return childDataIDs;
}

/* ---------------------------------------------------------------------- */

vector<int> Node::GetChildNodeIDs() const
{
  return childNodeIDs;
}

/* ---------------------------------------------------------------------- */

void Node::SetNodeID(int newID)
{
  nodeID=newID;
}

/* ---------------------------------------------------------------------- */

void Node::SetLeftChildNodeID(int newID)
{
  leftChildNodeID=newID;
}

/* ---------------------------------------------------------------------- */

void Node::SetRightChildNodeID(int newID)
{
  rightChildNodeID=newID;
}

/* ---------------------------------------------------------------------- */

void Node::FlagAsMerged()
{
  mergeFlag=true;
}

/* ---------------------------------------------------------------------- */

void Node::UnflagAsMerged()
{
  mergeFlag=false;
}

/* ----------------------------------------------------------------------
   Write the results from a vector of nodes in a format that is easily
   read into R.
---------------------------------------------------------------------- */

void Node::OutputResultsToFile(const vector<Node>& treeNode, string outputFileName)
{
  int i;
  ofstream outputFile(outputFileName.c_str());
  int nItems = (treeNode.size() + 1)/2;
  
  // Do the writing
  if (outputFile.is_open())
  {
    // There are always (nItems-1) merges in BHC
    for (i=nItems; i<2*nItems-1; i++)
    {
      outputFile << treeNode[i].GetRightChildNodeID()+1 << " "
                 << treeNode[i].GetLeftChildNodeID()+1  << " "
                 << i-nItems+1                         << " "
                 << treeNode[i].GetClusterLogEvidence() << endl;
    }
    outputFile.close();
  }
  else cout << "Unable to open file" << endl;
}

/* ---------------------------------------------------------------------- */
