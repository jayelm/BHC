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

pos_predict = function(mm, y) {

  predict_component = function(probs, logwk) {
    probs_y = ifelse(y, probs, 1 - probs)
    sum(log(probs_y)) + logwk
  }

  predictions = mapply(predict_component, mm$probs, mm$logwk)

  log_sum_exp(predictions)
}
