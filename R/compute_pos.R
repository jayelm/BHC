# Functions for adding posterior predictive distributions to trees from bhc

# PROBABILITY FUNCTIONS ====

# log1mexp(x) = log(1 - exp(-x))
log1mexp = function(x) {
  ifelse(x <= log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

# log1pexp(x) = log(1 + exp(x))
log1pexp = function(x) {
  ifelse(x <= -37, exp(x),
         ifelse(x <= 18, log1p(exp(x)),
                ifelse(x <= 33, x + exp(-x), x)))
}

# If given a log-odds ratio log(p / (1 - p)), this obtains log(p) with
# (hopefully) the least numerical instability issues
loglogistic = function(x) -log1pexp(-x)

# If x = log(p) and y = log(q), this computes log(p + q)
log_add = function(x, y) {
  if (x == -Inf) {
    y
  } else if (y == -Inf) {
    x
  } else {
    max(x, y) + log1pexp(-abs(x - y))
  }
}

# If L = log(p), p = exp(L), so
# log(1 - p) = log(1 - exp(L))
log_complement = function(x) log1mexp(-x)

# ADD WEIGHTS FUNCTIONS ====
# logEvidence is the log-merge ratio i.e. log-odds of r_k.
# This function computes log(r_k) given logEvidence
add_logrk = function(dend) {
  dendrapply(dend, function(n) {
    logEvidencek = attr(n, "logEvidence")
    logrk = loglogistic(logEvidencek)
    attr(n, "logrk") = logrk
    n
  })
}

# This function returns a new dendrogram with logwk assigned.
# lognk_weight is $\prod_{i \in \mathcal{N}_k} (1 - r_i)$ (i.e. product of
# complements of the weights from the parent to the root)
add_logwk = function(dend, lognk_weight = 0) {
  logrk = attr(dend, "logrk")
  attr(dend, "logwk") = logrk + lognk_weight
  # The complement of this node's r_k becomes part of lognk_weight
  new_lognk_weight = log_complement(logrk) + lognk_weight
  if (!is.leaf(dend)) {
    stopifnot(length(dend) == 2) # Enforce bifurcation
    dend[[1]] = add_logwk(dend[[1]], new_lognk_weight)
    dend[[2]] = add_logwk(dend[[2]], new_lognk_weight)
  }
  dend
}

add_weights = function(dend) {

  dend = add_logrk(dend)
  dend = add_logwk(dend)

  dend
}

compute_hyperparameters = function(dend, data) {
  global_hyperparameter = attr(dend, "globalHyperParam")

  n_data_items = nrow(data)
  n_features = ncol(data)
  n_feature_values = length(unique(data[1:length(data)]))

  hyperparameter = sapply(1:n_features, function(i) {
    table_i = table(data[, i])
    # XXX: Basic way of making sure ordering is correct. Probably not ok for
    # more complex categorical data?
    stopifnot(names(table_i)[1] < names(table_i)[2])

    # Convert to data_counter - with as_numeric and also add 1 to everything
    # (for prior)
    data_counter = as.numeric(table_i) + 1

    # Katherine's formula
    global_hyperparameter * data_counter / (n_data_items + 1)
  })

  hyperparameter
}

# FINAL COMPUTING FUNCTION ====
compute_pos = function(dend, data) {
  dend = add_weights(dend)

  hypers = compute_hyperparameters(dend, data)

  dend
}
