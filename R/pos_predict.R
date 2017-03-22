# Function for computing posterior predictions

# Add two numbers
log_sum_exp_single = function(u, v) {
  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
}

# Add multiple numbers
log_sum_exp = function(a) {
  n = length(a)
  if (n == 0) {
    -Inf
  } else if (n == 1) {
    a[1]
  } else {
    Reduce(log_sum_exp_single, a[-1], a[1])
  }
}

# Note that the loglikelihood for an entire point is simply the product
# across the likelihood of all its features, where each of its features
# is the Dirichlet-Multinomial model.
# Compute the loglikelihood of the given test point for the given dendrogram
pos_predict = function(dend, y) {
  # Cache feature vector as number of columns of root logprobs
  n_features = ncol(attr(dend, "logprobs"))
  # Cache indexing vector into logprob matrices
  # Add 1 so that the values are valid indices
  # FIXME: This only works with multinomial values that start from 0!!
  lp_ix = cbind(y + 1, 1:n_features)

  rec_predict = function(dend) {
    # Get these logprobs and weight
    logprobs = attr(dend, "logprobs")
    logwk = attr(dend, "logwk")

    logprobs_y = logprobs[lp_ix]

    # For a single component, the probability is the sum of the selected
    # logprobs, times the weight of the cluster logwk.
    logp = sum(logprobs_y) + logwk
    if (!is.leaf(dend)) {
      # Get probabilities of children
      logp_1 = rec_predict(dend[[1]])
      logp_2 = rec_predict(dend[[2]])
      # Add these probabilities and try to prevent instability
      logp = log_sum_exp(c(logp, logp_1, logp_2))
    }
    logp
  }

  rec_predict(dend)
}
