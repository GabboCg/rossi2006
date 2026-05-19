gmm_beta <- function(y, z, x, heter = 1) {
  
  # GMM with optimal weighting matrix (Hayashi, ch.3)
  # y: T x 1 dependent variable
  # z: T x l independent variable
  # x: T x k matrix of instruments
  # heter: 1 for heteroskedasticity-robust, 0 for homoskedastic
  
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
  
  return(as.vector(b2))
  
}
