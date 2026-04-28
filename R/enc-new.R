enc_new <- function(actual, forecast_1, forecast_2) {

  u_1 <- actual - forecast_1
  u_2 <- actual - forecast_2

  enc <- nrow(u_1) * (mean(u_1 ^ 2 - u_1 * u_2) / mean(u_2 ^ 2))

  output <- list(enc)
  names(output) <- c("enc-new")

  return(output)

}