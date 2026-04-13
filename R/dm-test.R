dm_test <- function(actual, forecast_1, forecast_2) {
  
  e_1 <- (actual - forecast_1) ^ 2
  e_2 <- (actual - forecast_2) ^ 2
  
  f_hat <- e_2 - e_1
  
  var_f <- var(f_hat) / nrow(f_hat)
  
  dm <- mean(f_hat) / sqrt(var_f)
  
  p_value <- 1 - pchisq(dm ^ 2, 1)
  
  output <- list(dm, p_value)
  names(output) <- c('dm', 'p-value')
  
  return(output)
  
}