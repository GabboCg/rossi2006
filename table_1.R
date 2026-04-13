#'---
#' author: Gabriel E. Cabrera
#' title: Comparing Out-Of-Sample Forecasts
#' date: Last updated `r Sys.Date()`
#'---

# packages
library("readxl")
library("dplyr")

start_time <- Sys.time()

# auxiliary functions
source("Functions/dm-test.R")
source("Functions/enc-new.R")

# Read Rossi (2005) dataset
predictor_data <- read_xls("Dataset/PredictorData1998.xls")

# Out-Of-Sample Forecasting

ys <- predictor_data %>% 
  janitor::clean_names() %>%
  select(date, 9, 14, 19, 24, 29) %>% 
  mutate_at(vars(2:6), list(~ log(. / lag(., 1)))) %>%  
  rename("can" = 2, "fr" = 3, "ger" = 4, "it" = 5, "jap" = 6) 

n_row <- nrow(ys) - 3
n_col <- ncol(ys)

r <- round(0.5 * n_row)

p <- n_row - r

fc_exchange <- matrix(NA, p, 5)
colnames(fc_exchange) <- c("can", "fr", "ger", "it", "jp")

actual <- matrix(NA, p, 5)
colnames(actual) <- c("can", "fr", "ger", "it", "jp")

fc_rw <- matrix(NA, p, 1)
colnames(fc_rw) <- c("rw")

for (i in 1:(n_col - 1)) {
  
  # generate dataset with lag
  xs_i <- ys[,(i + 1)] %>%
    rename("y_i" = 1) %>%  
    mutate(x_i = lag(y_i, 1)) %>% 
    na.omit() %>% 
    slice(-1)
    
  for (j in 1:p) {
    
    cat("Iteration:", j, "of", p, "\n")
    
    # fit ar1
    ar1_i_j <- lm(y_i ~ ., data = xs_i[1:(r + j - 1),]) 
    
    # out-of-sample 
    xnew <- xs_i[(r + j), 2] 
    
    # out-of-sample prediction
    fc_exchange[j, i] <- predict(ar1_i_j, xnew)
    
    # actual
    actual[j, i] <- xs_i[(r + j), 1] %>% 
      pull() 
    
    # random walk
    fc_rw[j, 1] <- 0
    
  }

}

# Compare Out-Of-Sample Forecast 

oosres <- matrix(NA, 5, 3)

colnames(oosres) <- c("dm", "p-value", "enc-new")
rownames(oosres) <- colnames(ys)[-1]

for (i in 1:5) {
  
  dm_res <- dm_test(actual[, i, drop = FALSE], fc_rw[,1], fc_exchange[, i, drop = FALSE])
  
  oosres[i, 1] <- dm_res$dm
  oosres[i, 2] <- dm_res$`p-value`
  
  enc_new_res <- enc_new(actual[, i, drop = FALSE], fc_rw[,1], fc_exchange[, i, drop = FALSE])
  
  oosres[i, 3] <- enc_new_res$`enc-new`
  
}

# Replicate part of Table 1 > AR(1) and AR(2) models
round(oosres, 2)

end_time <- Sys.time()

elasped_time <- end_time - start_time

cat('\n')
cat(sprintf('%s%#.2f\n', 'elapsed time = ', elasped_time))
