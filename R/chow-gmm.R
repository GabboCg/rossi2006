chow_gmm <- function(y, x, z, w, t_break, heter = 1) {
  
  # Chow test statistic by GMM for a fixed time-break
  # y: T x 1 dependent variable
  # x: T x p variables whose params are time-varying (tested)
  # z: NULL if no subsets, otherwise T x q variables whose params are not time-varying
  # w: T x k instruments
  # t_break: time of the break (index, not fraction)

  n_obs <- length(y)
  p <- ncol(x)
  k <- ncol(w)
  pi_val <- t_break / n_obs

  if (is.null(z)) {
    
    # No subsets case
    y1 <- y[1:t_break]
    x1 <- x[1:t_break, , drop = FALSE]
    w1 <- w[1:t_break, , drop = FALSE]
    y2 <- y[(t_break + 1):n_obs]
    x2 <- x[(t_break + 1):n_obs, , drop = FALSE]
    w2 <- w[(t_break + 1):n_obs, , drop = FALSE]

    X <- rbind(
      cbind(x1, matrix(0, t_break, p)),
      cbind(matrix(0, n_obs - t_break, p), x2)
    )
    W <- rbind(
      cbind(w1, matrix(0, t_break, k)),
      cbind(matrix(0, n_obs - t_break, k), w2)
    )

    b <- gmm_beta(y, X, W, heter)
    b1 <- b[1:p]
    b2 <- b[(p + 1):(2 * p)]
    e <- gmm_res(y, X, W, heter)

    M1 <- crossprod(w1, x1) / t_break
    M2 <- crossprod(w2, x2) / (n_obs - t_break)

    if (heter == 1) {
      
      S1 <- crossprod(w1 * e[1:t_break]) / t_break
      S2 <- crossprod(w2 * e[(t_break + 1):n_obs]) / (n_obs - t_break)
      
    } else {
      
      s2 <- sum(e ^ 2) / n_obs
      S1 <- s2 * crossprod(w1) / t_break
      S2 <- s2 * crossprod(w2) / (n_obs - t_break)
      
    }

    V1 <- solve(t(M1) %*% solve(S1) %*% M1)
    V2 <- solve(t(M2) %*% solve(S2) %*% M2)
    wald <- n_obs * t(b1 - b2) %*% solve(V1 / pi_val + V2 / (1 - pi_val)) %*% (b1 - b2)
    
  } else {
    
    # Subsets case
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

    M1 <- crossprod(w1, x1) / t_break
    M2 <- crossprod(w2, x2) / (n_obs - t_break)

    if (heter == 1) {
      
      S1 <- crossprod(w1 * e[1:t_break]) / t_break
      S2 <- crossprod(w2 * e[(t_break + 1):n_obs]) / (n_obs - t_break)
      
    } else {
      
      s2 <- sum(e^2) / n_obs
      S1 <- s2 * crossprod(w1) / t_break
      S2 <- s2 * crossprod(w2) / (n_obs - t_break)
      
    }

    V1 <- solve(t(M1) %*% solve(S1) %*% M1)
    V2 <- solve(t(M2) %*% solve(S2) %*% M2)
    wald <- n_obs * t(b1 - b2) %*% solve(V1 / pi_val + V2 / (1 - pi_val)) %*% (b1 - b2)
    
  }

  return(as.numeric(wald))
  
}
