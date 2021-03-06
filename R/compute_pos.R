# Functions for adding posterior predictive distributions to trees from bhc

LOG2 = log(2)

# PROBABILITY FUNCTIONS ====
# Note these functions are NOT vectorized with ifelse as they don't work with
# mpfr numbers.

# log1mexp(x) = log(1 - exp(-x))
log1mexp = function(x) {
  if (x <= LOG2) {
    log(-expm1(-x))
  } else {
    log1p(-exp(-x))
  }
}

# log1pexp(x) = log(1 + exp(x))
# TODO: These cutoffs might not be the same based on varying PREC
log1pexp = function(x) {
  if (x <= -37) {
    exp(x)
  } else if (x <= 18) {
    log1p(exp(x))
  } else if (x <= 33) {
    x + exp(-x)
  } else {
    x
  }
}

# If given a log-odds ratio log(p / (1 - p)), this obtains log(p) with
# (hopefully) the least numerical instability issues
loglogistic = function(x) -log1pexp(-x)

# Add two numbers
# TODO: Figure out how to share this with pos_predict
log_sum_exp_single = function(u, v) {
  max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)))
}

# If L = log(p), p = exp(L), so
# log(1 - p) = log(1 - exp(L))
log_complement = function(x) log1mexp(-x)

# ADD WEIGHTS FUNCTIONS ====
# logEvidence is the log-merge ratio i.e. log-odds of r_k.
# This function computes log(r_k) given logEvidence
add_logrk = function(dend) {
  dendrapply(dend, function(node) {
    logEvidencek = attr(node, "logEvidence")
    logrk = loglogistic(logEvidencek)
    attr(node, "logrk") = logrk
    node
  })
}

# This function returns a new dendrogram with logwk assigned.
# lognk_weight is $\prod_{i \in \mathcal{N}_k} (1 - r_i)$ (i.e. product of
# complements of the weights from the parent to the root)
add_logwk = function(dend, lognk_weight = 0, method = 'nk', root_members = 0) {
  if (!(method %in% c('nk', 'eq9', 'tim'))) {
    stop("Weights must be one of nk, eq9, tim")
  }
  if (root_members == 0) {
    # Set # members at root
    root_members = attr(dend, 'members')
  }

  logrk = attr(dend, "logrk")

  if (method == 'tim') {
    # r_k, lognk, AND tim's formula.
    # alpha hardcoded to 0.001
    log_comp_prob = log(attr(dend, 'members') / (root_members + 0.001))
    attr(dend, "logwk") = logrk + lognk_weight + log_comp_prob
  } else {
    # Standard r_k + lognk
    attr(dend, "logwk") = logrk + lognk_weight
  }
  if (!is.leaf(dend)) {
    # The complement of this node's r_k becomes part of lognk_weight
    new_lognk_weight = log_complement(logrk) + lognk_weight
    stopifnot(length(dend) == 2) # Enforce bifurcation

    if (method == 'eq9') {
      new_lognk_weight_left = new_lognk_weight_right = new_lognk_weight
    } else if (method == 'nk') {
      # This is where I deviate from Heller - also encoded in lognk_weight
      # should be the probability of NOT choosing the other child. This is
      # weighted by the size of the subtrees
      n_k = attr(dend, "members")
      n_left = attr(dend[[1]], "members")
      R_left = n_left / n_k

      # For the LEFT child, include in lognk_weight the probability of choosing
      # the left (i.e. log(R_left))
      # for the RIGHT child, include the probability of choosing the right (1 -
      # logp_left)
      new_lognk_weight_left = new_lognk_weight + log(R_left)
      new_lognk_weight_right = new_lognk_weight + log(1 - R_left)
    } else {
      # Tim's formula operates not on children, but on the current node.
      new_lognk_weight_left = new_lognk_weight_right = new_lognk_weight
    }

    dend[[1]] = add_logwk(dend[[1]], new_lognk_weight_left, root_members = root_members, method = method)
    dend[[2]] = add_logwk(dend[[2]], new_lognk_weight_right, root_members = root_members, method = method)
  }
  dend
}

add_weights = function(dend, method = 'nk') {
  dend = add_logrk(dend)
  dend = add_logwk(dend, method = method)

  dend
}

add_data_indices = function(dend, data) {
  # Construct environment mapping dendrogram labels
  # FIXME: This doesn't work when dendrograms have the same label!!
  # Will need a (probably much more drastic) workaround!! E.g.
  # Dive into C++
  dend_labs = labels(dend)
  data_labs = rownames(data)

  # Enforce labels uniqueness
  stopifnot(length(unique(dend_labs)) == length(dend_labs) &&
            length(dend_labs) == length(data_labs) &&
            length(data_labs) == length(unique(data_labs)))
  # Enforce labels matching
  stopifnot(all(sort(dend_labs) == sort(data_labs)))

  data_indices = seq_along(data_labs)
  labs_list = as.list(setNames(data_indices, data_labs))

  labs_map = list2env(labs_list, hash = TRUE)

  # Recursively traverse the tree, keeping track of labels
  add_indices = function(dend) {
    this_lab = attr(dend, "label")
    if (is.leaf(dend)) {
      attr(dend, "ix") = c(labs_map[[this_lab]])
    } else {
      dend[[1]] = add_indices(dend[[1]])
      dend[[2]] = add_indices(dend[[2]])
      d1_ix = attr(dend[[1]], "ix")
      d2_ix = attr(dend[[2]], "ix")
      attr(dend, "ix") = c(d1_ix, d2_ix)
    }
    dend
  }

  dend_ix = add_indices(dend)
  # Just some final assertions on the root of the tree
  ix = attr(dend_ix, "ix")
  stopifnot(length(unique(ix)) == length(ix))
  stopifnot(all(sort(ix) == seq_along(ix)))

  dend_ix
}

# FIXME: Could run faster - might not need inner sapply. Do apply over table,
# then multiply entire thing by scalars
compute_hyperparameters = function(dend, data) {
  global_hyperparameter = attr(dend, "globalHyperParam")
  if (is.null(global_hyperparameter)) {
    # Sensible (inconsequential) default, same as Expt 1.
    global_hyperparameter = 0.46
  }

  n_data_items = nrow(data)
  n_features = ncol(data)

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

  # Transpose since sapply returns columns
  hyperparameter
}

compute_probs = function(dend, data, hypers) {
  feature_values = sort(unique(data[1:length(data)]))

  # First add raw counts as aprime
  add_counts = function(dend) {
    if (is.leaf(dend)) {
      ix = attr(dend, "ix")
      stopifnot(length(ix) == 1)
      # Factor w/ feature_values levels to assign zero counts
      row = factor(data[ix, ], levels = feature_values)
      # Trivially tabulate each column of the row and convert to numeric
      counts = sapply(row, function(t) table(t))
      attr(dend, "aprime") = counts
    } else {
      # Get counts of children
      dend[[1]] = add_counts(dend[[1]])
      dend[[2]] = add_counts(dend[[2]])
      d1_counts = attr(dend[[1]], "aprime")
      d2_counts = attr(dend[[2]], "aprime")

      # Sum up counts of children, including hyperparams
      counts = d1_counts + d2_counts
      attr(dend, "aprime") = counts
    }
    dend
  }

  # Compute probs and drop aprime
  add_probs = function(dend) {
    dendrapply(dend, function(node) {
      # XXX: here I extract raw probabilities. If I wanted to save even more
      # space, I could just store the true probability (for the binary case).
      # But maybe premature optimization
      aprime = attr(node, "aprime") + hypers
      # Normalize by column sums, then take the log
      aprime_scaled = scale(aprime, center = FALSE, scale = colSums(aprime))
      attr(node, "aprime") = NULL
      # FIXME: THIS ONLY WORKS WITH BINARY MODELS - I only keep positive
      # probability (i.e. 1). Row 1 is p(0), row 2 is p(1)
      attr(node, "probs") = aprime_scaled[2, ]
      node
    })
  }

  dend_counts = add_counts(dend)
  dend_probs = add_probs(dend_counts)

  dend_probs
}

extract_mixture_model = function(dend) {
  # To facilitate faster computation with pos_predict, we can now drop the
  # dendrogram structure and return a list with the
  # XXX: This will only work with binaries!
  probs = dendextend::get_nodes_attr(dend, "probs")
  # Aprime is a 3d array, we want to collapse to a
  logwk = dendextend::get_nodes_attr(dend, "logwk")

  list(probs = data.frame(probs), logwk = logwk)
}

# FINAL COMPUTING FUNCTION ====
compute_pos = function(dend, data, verbose = FALSE, method = 'nk') {
  # Make sure data is a matrix
  data = as.matrix(data)
  if (verbose) {
    message("Adding weights (logrk and logwk)")
    message("logwk computed according to ", method)
  }
  dend = add_weights(dend, method = method)
  if (verbose) message("Adding indices of data")
  dend = add_data_indices(dend, data)

  if (verbose) message("Computing hyperparameters")
  hypers = compute_hyperparameters(dend, data)

  if (verbose) message("Computing probs")
  dend = compute_probs(dend, data, hypers)

  if (verbose) message("Extracting mixture model components only")
  mm = extract_mixture_model(dend)

  if (verbose) message("Done, returning")
  mm
}
