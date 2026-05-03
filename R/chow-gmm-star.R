chow_gmm_star <- function(y, x, z, w, t_break, heter = 1) {
  
  # Optimal Chow* test statistic by GMM for a fixed time-break
  # y: T x 1 dependent variable
  # x: T x p variables whose params are time-varying (tested)
  # z: NULL if no subsets, otherwise T x q variables not time-varying
  # w: T x k instruments
  # t_break: time of the break

  n_obs <- length(y)
  p <- ncol(x)
  k <- ncol(w)
  pi_val <- t_break / n_obs

  if (is.null(z)) {
    
    y1 <- y[1:t_break]
    x1 <- x[1:t_break, , drop = FALSE]
    w1 <- w[1:t_break, , drop = FALSE]
    y2 <- y[(t_break + 1):n_obs]
    x2 <- x[(t_break + 1):n_obs, , drop = FALSE]
    w2 <- w[(t_break + 1):n_obs, , drop = FALSE]

    b1 <- gmm_beta(y1, x1, w1, heter)
    b2 <- gmm_beta(y2, x2, w2, heter)
    e1 <- gmm_res(y1, x1, w1, heter)
    e2 <- gmm_res(y2, x2, w2, heter)

    if (heter == 1) {
      
      Sigma1 <- crossprod(w1 * e1) / t_break
      Sigma2 <- crossprod(w2 * e2) / (n_obs - t_break)
      
    } else {
      
      s2 <- (sum(e1^2) + sum(e2^2)) / n_obs
      Sigma1 <- s2 * crossprod(w1) / t_break
      Sigma2 <- s2 * crossprod(w2) / (n_obs - t_break)
      
    }

    Gamma <- solve(rbind(
      cbind(pi_val * Sigma1, matrix(0, k, k)),
      cbind(matrix(0, k, k), (1 - pi_val) * Sigma2)
    ))

    Swx1 <- crossprod(w1, x1) / n_obs
    Swx2 <- crossprod(w2, x2) / n_obs
    
    M <- rbind(
      cbind(Swx1, matrix(0, k, p)),
      cbind(matrix(0, k, p), Swx2)
    )

    R <- rbind(
      cbind(diag(p), -diag(p)),
      cbind(pi_val * diag(p), (1 - pi_val) * diag(p))
    )

    VRb <- R %*% solve(t(M) %*% Gamma %*% M) %*% t(R)
    beta_vec <- c(b1 - b2, pi_val * b1 + (1 - pi_val) * b2)
    wald <- n_obs * t(beta_vec) %*% solve(VRb) %*% beta_vec
    
  } else {
    
    q <- ncol(z)
    y1 <- y[1:t_break]
    x1 <- x[1:t_break, , drop = FALSE]
    z1 <- z[1:t_break, , drop = FALSE]
    w1 <- w[1:t_break, , drop = FALSE]
    y2 <- y[(t_break + 1):n_obs]
    x2 <- x[(t_break + 1):n_obs, , drop = FALSE]
    z2 <- z[(t_break + 1):n_obs, , drop = FALSE]
    w2 <- w[(t_break + 1):n_obs, , drop = FALSE]

    X <- rbind(
      cbind(x1, matrix(0, t_break, p), z1),
      cbind(matrix(0, n_obs - t_break, p), x2, z2)
    )
    W <- rbind(
      cbind(w1, matrix(0, t_break, k)),
      cbind(matrix(0, n_obs - t_break, k), w2)
    )

    b <- gmm_beta(y, X, W, heter)
    b1 <- b[1:p]
    b2 <- b[(p + 1):(2 * p)]
    e <- gmm_res(y, X, W, heter)

    if (heter == 1) {
      
      Sigma1 <- crossprod(w1 * e[1:t_break]) / t_break
      Sigma2 <- crossprod(w2 * e[(t_break + 1):n_obs]) / (n_obs - t_break)
      
    } else {
      
      s2 <- sum(e ^ 2) / n_obs
      Sigma1 <- s2 * crossprod(w1) / t_break
      Sigma2 <- s2 * crossprod(w2) / (n_obs - t_break)
      
    }

    Gamma <- solve(rbind(
      cbind(pi_val * Sigma1, matrix(0, k, k)),
      cbind(matrix(0, k, k), (1 - pi_val) * Sigma2)
    ))

    Swx1 <- crossprod(w1, x1) / n_obs
    Swz1 <- crossprod(w1, z1) / n_obs
    Swx2 <- crossprod(w2, x2) / n_obs
    Swz2 <- crossprod(w2, z2) / n_obs
    
    M <- rbind(cbind(Swx1, matrix(0, k, p), Swz1), cbind(matrix(0, k, p), Swx2, Swz2))

    R <- rbind(
      cbind(diag(p), -diag(p), matrix(0, p, q)),
      cbind(pi_val * diag(p), (1 - pi_val) * diag(p), matrix(0, p, q))
    )

    VRb <- R %*% solve(t(M) %*% Gamma %*% M) %*% t(R)
    beta_vec <- c(b1 - b2, pi_val * b1 + (1 - pi_val) * b2)
    wald <- n_obs * t(beta_vec) %*% solve(VRb) %*% beta_vec
    
  }

  return(as.numeric(wald))
  
}
