F.stat.twosample <- function(
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
  
  var1 <- apply(x, 2, var)
  var2 <- apply(y, 2, var)
  
  dt <- diff(range(t.seq)) / (length(t.seq)-1)
  
  numerator <- sum((mean1 - mean2)^2) * dt
  denominator <- (sum(var1) + sum(var2)) * dt
  
  stat <- numerator / denominator
  
  if(method == 1) {
    
    return(list(
      statistic = stat,
      pvalue = NA
    ))
  }
  
  if(method == 2) {
    
    combined <- rbind(x, y)
    
    perm.stats <- numeric(replications)
    
    for(i in 1:replications) {
      
      idx <- sample(1:(n1+n2))
      
      g1 <- combined[idx[1:n1], ]
      g2 <- combined[idx[(n1+1):(n1+n2)], ]
      
      m1 <- colMeans(g1)
      m2 <- colMeans(g2)
      
      v1 <- apply(g1, 2, var)
      v2 <- apply(g2, 2, var)
      
      num <- sum((m1 - m2)^2) * dt
      den <- (sum(v1) + sum(v2)) * dt
      
      perm.stats[i] <- num / den
    }
    
    pval <- mean(perm.stats >= stat)
    
    hist(perm.stats,
         main = "Permutation distribution",
         xlab = "F statistic")
    
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