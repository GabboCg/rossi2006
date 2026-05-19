gmm_var <- function(y, z, x, heter = 1) {
  
  # Variance of GMM coefficient vector with optimal weighting matrix
  # Returns asymptotic variance matrix
  
  n_obs <- length(y)
  Sxx <- crossprod(x) / n_obs
  Sxz <- crossprod(x, z) / n_obs
  Sxy <- crossprod(x, y) / n_obs
  
  # 1st stage
  W1 <- solve(Sxx)
  b1 <- solve(t(Sxz) %*% W1 %*% Sxz) %*% (t(Sxz) %*% W1 %*% Sxy)
  e <- y - z %*% b1
  
  if (heter == 1) {
    
    S <- crossprod(x * as.vector(e)) / n_obs
    W2 <- solve(S)
    b2 <- solve(t(Sxz) %*% W2 %*% Sxz) %*% (t(Sxz) %*% W2 %*% Sxy)
    avarb <- solve(t(Sxz) %*% W2 %*% Sxz)
    
  } else {
    
    S <- solve(Sxx) * sum(e ^ 2) / n_obs
    W2 <- solve(S)
    b2 <- solve(t(Sxz) %*% W2 %*% Sxz) %*% (t(Sxz) %*% W2 %*% Sxy)
    avarb <- solve(Sxx) * sum((y - z %*% b2) ^ 2) / n_obs
    
  }
  
  return(avarb)
  
}
