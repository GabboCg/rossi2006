nyblom_star <- function(y, x, z, w, heter = 1, b0) {
  
  # Optimal Nyblom test
  # y: T x 1 dependent variable
  # x: T x p variables whose params are time-varying
  # z: NULL if no subsets, otherwise T x q unrestricted variables
  # w: T x k instruments
  # b0: null hypothesis value for beta
  
  n_obs <- length(y)
  
  if (is.null(z)) {
    
    e <- y - x %*% b0
    k <- ncol(x)
    
    if (heter == 1) {
      
      Sigma <- crossprod(w * as.vector(e)) / n_obs
      
    } else {
      
      s2 <- sum(e ^ 2) / n_obs
      Sigma <- s2 * crossprod(w) / n_obs
      
    }
    
    we <- w * as.vector(e)
    es <- apply(we, 2, cumsum) / sqrt(n_obs)
    Sigma_inv <- solve(Sigma)
    
    sum1 <- 0
    for (i in 1:n_obs) {
      
      sum1 <- sum1 + es[i, ] %*% Sigma_inv %*% es[i, ] / n_obs
    }
    
    return(as.numeric(sum1))
    
  } else {
    
    e_adj <- y - x %*% b0
    e <- gmm_res(e_adj, z, w, heter)
    k <- ncol(z)
    
    if (heter == 1) {
      
      Sigma <- crossprod(w * e) / n_obs
      
    } else {
      
      e2sum <- sum(e^2)
      s2 <- e2sum / (n_obs - k)
      Sigma <- s2 * crossprod(w) / n_obs
      
    }
    
    Sigma_sqrt <- expm::sqrtm(Sigma)
    Sigma_sqrt_inv <- solve(Sigma_sqrt)
    
    Mbetabar <- Sigma_sqrt_inv %*% (-crossprod(w, x) / n_obs)
    Mdeltabar <- Sigma_sqrt_inv %*% (-crossprod(w, z) / n_obs)
    
    Pbardelta <- Mdeltabar %*% solve(crossprod(Mdeltabar)) %*% t(Mdeltabar)
    omegaN <- t(Mbetabar) %*% (diag(nrow(Pbardelta)) - Pbardelta) %*% Mbetabar
    
    GradQ <- matrix(0, n_obs, ncol(x))
    
    for (i in 1:n_obs) {
      
      GradQ[i, ] <- e[i] * w[i, ] %*% Sigma_sqrt_inv %*% Mbetabar
      
    }
    
    GradQpiT <- apply(GradQ, 2, cumsum) / sqrt(n_obs)
    omegaN_inv <- solve(omegaN)
    
    sum1 <- 0
    
    for (i in 1:n_obs) {
      
      sum1 <- sum1 + GradQpiT[i, ] %*% omegaN_inv %*% GradQpiT[i, ] / n_obs
      
    }
    
    return(as.numeric(sum1))
    
  }
  
}
