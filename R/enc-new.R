enc_new <- function(actual, forecast_1, forecast_2) {
  
  e_1 <- (actual - forecast_1) ^ 2
  e_2 <- (actual - forecast_2) ^ 2
  
  enc <- nrow(e_1) * (mean(e_1 ^ 2  - e_2 * e_1) / mean((e_2 - mean(e_2)) ^ 2))
  
  output <- list(enc)
  names(output) <- c("enc-new")
  
  return(output)
  
}