Ztwosample <- function(x, y, t.seq, alpha = 0.05) {
  
  # Convert fd objects to matrices if needed
  if("fd" %in% class(x)) {
    x <- eval.fd(t.seq, x)
  }
  
  if("fd" %in% class(y)) {
    y <- eval.fd(t.seq, y)
  }
  
  x <- t(x)
  y <- t(y)
  
  n1 <- nrow(x)
  n2 <- nrow(y)
  
  p <- ncol(x)
  
  mean1 <- colMeans(x)
  mean2 <- colMeans(y)
  
  var1 <- apply(x, 2, var)
  var2 <- apply(y, 2, var)
  
  zvals <- (mean1 - mean2) /
    sqrt(var1/n1 + var2/n2)
  
  pvals <- 2 * (1 - pnorm(abs(zvals)))
  
  significant <- pvals < alpha
  
  plot(t.seq, pvals,
       type = "l",
       ylim = c(0,1),
       main = "Pointwise Z-test p-values",
       xlab = "Time",
       ylab = "p-value")
  
  abline(h = alpha,
         col = "red",
         lty = 2)
  
  return(list(
    z = zvals,
    pvalue = pvals,
    significant = significant
  ))
}