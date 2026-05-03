gmm_res <- function(y, z, x, heter = 1) {
  
  # GMM residuals with optimal weighting matrix

  n_obs <- length(y)
  Sxx <- crossprod(x) / n_obs
  Sxz <- crossprod(x, z) / n_obs
  Sxy <- crossprod(x, y) / n_obs

  # 1st stage
  W1 <- solve(Sxx)
  b1 <- solve(t(Sxz) %*% W1 %*% Sxz) %*% (t(Sxz) %*% W1 %*% Sxy)
  e <- y - z %*% b1
  S <- crossprod(x * as.vector(e)) / n_obs
  W2 <- solve(S)

  # 2nd stage
  b2 <- solve(t(Sxz) %*% W2 %*% Sxz) %*% (t(Sxz) %*% W2 %*% Sxy)

  if (heter == 1) {
    
    b <- b2
    
  } else {
    
    b <- b1
    
  }

  res <- as.vector(y - z %*% b)
  
  return(res)
  
}
