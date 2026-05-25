L2.stat.twosample <- function(
    x,
    y,
    t.seq,
    method = 1,
    replications = 1000
) {
  
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
  
  mean1 <- colMeans(x)
  mean2 <- colMeans(y)
  
  dt <- diff(range(t.seq)) / (length(t.seq)-1)
  
  stat <- sum((mean1 - mean2)^2) * dt
  
  # Asymptotic approximation
  if(method == 1) {
    
    return(list(
      statistic = stat,
      pvalue = NA
    ))
  }
  
  # Permutation test
  if(method == 2) {
    
    combined <- rbind(x, y)
    
    perm.stats <- numeric(replications)
    
    for(i in 1:replications) {
      
      idx <- sample(1:(n1+n2))
      
      g1 <- combined[idx[1:n1], ]
      g2 <- combined[idx[(n1+1):(n1+n2)], ]
      
      m1 <- colMeans(g1)
      m2 <- colMeans(g2)
      
      perm.stats[i] <- sum((m1 - m2)^2) * dt
    }
    
    pval <- mean(perm.stats >= stat)
    
    hist(perm.stats,
         main = "Permutation distribution",
         xlab = "L2 statistic")
    
    abline(v = stat,
           col = "red",
           lwd = 2)
    
    return(list(
      statistic = stat,
      pvalue = pval,
      perm.stats = perm.stats
    ))
  }
}