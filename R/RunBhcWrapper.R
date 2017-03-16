##Wrapper function to run the interfaced C++ code, ensuring correct formating etc
RunBhcWrapper <- function(globalHyperParameter, dataTypeID, data, timePoints, nDataItems, nFeatures, nFeatureValues, noise, numReps=0, noiseMode=0, robust=0, fullOutputSwitch=FALSE, numThreads=1, randomised=FALSE, m=2, verbose=FALSE){
  ##prep
  if(!randomised)
    m = -1 # this is how we tell the library not to use the randomised algorithm
  
  ##generate the output structure
  if(dataTypeID==0){ # Multinomial case: use the 1.1.0 code
    # JM: Note bhcWrapper_multinomial always takes pointers
    # JM: ??? bhcWrapper_multinomial doesn't return anything!
    out <- .C("bhcWrapper_multinomial",
              # JM: Matrix of data?
              as.integer(data),
              # JM: Number of rows in the matrix
              as.integer(nDataItems),
              # JM: Number of columns
              as.integer(nFeatures),
              # JM: Found in FindOptimalHyperparameter.R. I also return this out of bhc.R
              as.double(globalHyperParameter),
              # JM: Number of discrete values used for multinomial
              as.integer(nFeatureValues),
              # JM: Dummy; this value will be changed by `bayeslink_binf`
              logEvidence=as.double(123),
              node1=vector(mode='integer',length=nDataItems-1),
              node2=vector(mode='integer', length=nDataItems-1),
              mergeOrder=vector(mode='integer', length=nDataItems-1),
              mergeWeight=vector(mode='numeric', length=nDataItems-1))
  }else
  {
    out <- .C("bhcWrapper",
              as.integer(dataTypeID),
              as.double(data), ##the data(matrix) is input as a vector read down the matrix columns; only the data is input not the row or column names
              as.double(timePoints),
              as.integer(nDataItems),
              as.integer(nFeatures),
              as.double(globalHyperParameter),
              as.double(noise),
              as.integer(numReps),
              as.integer(noiseMode),
              as.integer(robust),
              as.integer(nFeatureValues),
              logEvidence=as.double(123),
              node1=vector(mode='integer',length=nDataItems-1),
              node2=vector(mode='integer', length=nDataItems-1),
              mergeOrder=vector(mode='integer', length=nDataItems-1),
              mergeWeight=vector(mode='numeric', length=nDataItems-1),
              as.integer(numThreads),
              as.integer(m))
  }
  ##optionally, return the logEvidence so we can use optimisation routines
  if (verbose) print(c(globalHyperParameter, out$logEvidence), quote=FALSE)
  if (fullOutputSwitch) out else out$logEvidence
}
##*****************************************************************************
##*****************************************************************************
