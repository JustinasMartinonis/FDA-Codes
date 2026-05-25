trace <- function(x, t.seq) {
  
  # x = matrix/data of functions
  # rows = observations
  # cols = time points
  
  x <- as.matrix(x)
  
  # covariance matrix
  cv <- cov(x)
  
  # step size
  dt <- diff(range(t.seq)) / (length(t.seq) - 1)
  
  # numerical approximation of trace integral
  tr <- sum(diag(cv)) * dt
  
  return(tr)
}