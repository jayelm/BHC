# Function for computing posterior predictions

log_sum_exp = function(u, v) {
  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
}

# Note that the loglikelihood for an entire point is simply the product
# across the likelihood of all its features, where each of its features
# is the Dirichlet-Multinomial model.
# Compute the loglikelihood of the given test point for the given dendrogram
pos_predict = function(dend, y) {
  # Cache feature vector as number of columns of root aprime
  n_features = ncol(attr(dend, "aprime"))
  # Cache indexing vector into aprime matrices
  # Add 1 so that the values are valid indices
  # FIXME: This only works with multinomial values that start from 0!!
  ap_ix = cbind(y + 1, 1:n_features)

  rec_predict = function(dend) {
    # Get this aprime and weight
    aprime = attr(dend, "aprime")
    logwk = attr(dend, "logwk")

    apsum = colSums(aprime)
    apy = aprime[ap_ix]

    # For a single component, the probability is the sum of the selected
    # alphas (normalized to probabilities), plus the weight of the cluster
    # logwk
    val = sum(log(apy / apsum)) + logwk
    if (!is.leaf(dend)) {
      # Get probabilities of children
      lp_1 = rec_predict(dend[[1]])
      lp_2 = rec_predict(dend[[2]])
      # Add these probabilities and try to prevent instability
      val = log_sum_exp(val, log_sum_exp(lp_1, lp_2))
    }
    val
  }

  rec_predict(dend)
}
