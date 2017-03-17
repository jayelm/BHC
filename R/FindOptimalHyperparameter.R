##Function to run a series of numerical optimsations on the overall log-Evidence, to find
##the optimal value for the global hyperparameter.
##This function assumes that the log-Evidence is well-behaved in the region of param=1
FindOptimalHyperparameter <- function(dataTypeID, data, timePoints, noise, nDataItems, nFeatures, nFeatureValues, numThreads, verbose=FALSE){
  if (verbose) print("Optimising global hyperparameter...", quote=FALSE)
  ##----------------------------------------------------------------------
  ## BASIC SCAN OF (1D) HYPERPARAMETER SPACE TO FIND ROUGH REGION WHERE --
  ## THE OPTIMAL VALUE LIVES ---------------------------------------------
  ##----------------------------------------------------------------------
  param.lower <- 1
  param.upper <- 1
  logEv.lower <- RunBhcWrapper(param.lower, dataTypeID, data, timePoints, nDataItems, nFeatures, nFeatureValues, noise=noise, numThreads = numThreads)  
  logEv.upper <- logEv.lower
  ##find a sensible lower bound >= 2^-maxcounter
  maxcounter <- 10
  counter <- 0
  repeat{
    param.lower <- 0.5 * param.lower
    logEv.new   <- RunBhcWrapper(param.lower, dataTypeID, data, timePoints,
                                 nDataItems, nFeatures, nFeatureValues, noise=noise,
                                 numThreads = numThreads)
    counter <- counter + 1
    if ((logEv.new - logEv.lower) < 0 || counter > maxcounter)
      break
    else
      logEv.lower = logEv.new
  }
  ##find a sensible upper bound <= 2^maxcounter
  counter <- 0
  repeat{ 
    param.upper <- 2 * param.upper
    logEv.new   <- RunBhcWrapper(param.upper, dataTypeID, data, timePoints,
                                 nDataItems, nFeatures, nFeatureValues, noise=noise,
                                 numThreads = numThreads)
    counter <- counter + 1
    if ((logEv.new - logEv.upper) < 0 || counter > maxcounter)
      break
    else
      logEv.upper = logEv.new
  }
  ##----------------------------------------------------------------------
  ## GIVEN A ROUGH RANGE, RUN A SMARTER OPTIMISATION ALGORITHM -----------
  ##----------------------------------------------------------------------
  ##the tolerance here is in log-Evidence units.  This means that 1 is a fairly
  ##sensible value, because values smaller than that descend into the analysis'
  ##uncertainty
  ##(1 -> probability ratio of 'e' -> approx a one-sigma significance (Z-score))
  optimalOutput <- optimise(RunBhcWrapper, interval =c(param.lower, param.upper),
                            dataTypeID, data, timePoints, nDataItems, nFeatures,
                            nFeatureValues, noise=noise, maximum = TRUE, tol=1.,
                            numThreads = numThreads,
                            verbose=verbose)

  ##return the optimal hyperparameter value
  optimalOutput$maximum
}
##*****************************************************************************
##*****************************************************************************
##----------------------------------------------------------------------
## ----------------------------------------
##----------------------------------------------------------------------

