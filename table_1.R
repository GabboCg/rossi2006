#!/usr/bin/env Rscript
# ======================================================== #
#
#                  Replication of Table 1
#         AR(1) and AR(2) models (Rossi, 2006)
#
#                 Gabriel E. Cabrera-Guzman
#                The University of Manchester
#
#                        Spring, 2026
#
#                https://gcabrerag.rbind.io
#
# ------------------------------ #
# email: gabriel.cabreraguzman@postgrad.manchester.ac.uk
# ======================================================== #

library("expm")

start_time <- Sys.time()

# MATLAB-style rounding (round half up, not banker's rounding)
mround <- function(x) floor(x + 0.5)

# Source helper functions
source("R/gmm-beta.R")
source("R/gmm-var.R")
source("R/gmm-res.R")
source("R/chow-gmm.R")
source("R/chow-gmm-star.R")
source("R/nyblom-star.R")
source("R/pv-calc.R")

# Load p-value tables
pvqlrsb  <- as.matrix(read.table("data-raw/pvqlrsb.txt"))
pvapisb  <- as.matrix(read.table("data-raw/pvapisb.txt"))
pvnybsb  <- as.matrix(read.table("data-raw/pvnybsb.txt"))
pvqlropt <- as.matrix(read.table("data-raw/pvqlropt.txt"))
pvapiopt <- as.matrix(read.table("data-raw/pvapiopt.txt"))
pvap0opt <- as.matrix(read.table("data-raw/pvap0opt.txt"))
pvnybopt <- as.matrix(read.table("data-raw/pvnybopt.txt"))

# Load data
data_raw <- as.matrix(read.table("data/PredictorData1998.txt"))
heter <- 1  # robust to conditional heteroskedasticity

countries <- c("Can.", "Fr.", "Ger.", "It.", "Jap.")

# ==========================================
#     Run tests for both AR(1) and AR(2)
# ------------------------------------------

run_tests <- function(ntable, data_full) {

  # Full results: 25 rows
  full_results <- matrix(NA, 25, 5)

  for (ii in 1:5) {

    data <- data_full

    # No sample truncation for AR models (Table 1):
    # AR models only use exchange rates, not fundamentals with missing data

    # Extract exchange rate (log level)
    e <- log(data[, 4 + 5 * (ii - 1) + 4])
    T_full <- length(e)

    # Build regressors based on model type
    if (ntable == 71) {
      
      # AR(1): e1 = diff(e, t vs t-1), p1 = lag of diff(e)
      e1 <- e[4:T_full] - e[3:(T_full - 1)]
      p1 <- matrix(e[3:(T_full - 1)] - e[2:(T_full - 2)], ncol = 1)
      listrestr <- c(1, 2)
      listunr <- NULL
      
    } else if (ntable == 72) {
      
      # AR(2): e1 = diff(e), p1 = [lag1, lag2] of diff(e)
      e1 <- e[4:T_full] - e[3:(T_full - 1)]
      p1 <- cbind(e[3:(T_full - 1)] - e[2:(T_full - 2)], e[2:(T_full - 2)] - e[1:(T_full - 3)])
      listrestr <- c(1, 2, 3)
      listunr <- NULL
      
    }

    n_obs <- length(e1)
    kk <- ncol(p1)
    p1plusc <- cbind(1, p1)
    p1restr <- p1plusc[, listrestr, drop = FALSE]
    
    if (is.null(listunr)) {
      
      p1unr <- NULL
      
    } else {
      
      p1unr <- p1plusc[, listunr, drop = FALSE]
      
    }

    restr <- length(listrestr)

    # ==========================================
    #                   LR Test
    # ------------------------------------------
    z_full <- cbind(1, p1)
    bhat <- gmm_beta(e1, z_full, z_full, heter)
    varb <- gmm_var(e1, z_full, z_full, heter)
    iden <- diag(ncol(p1) + 1)
    RR <- iden[listrestr, , drop = FALSE]
    LLR <- n_obs * t(RR %*% bhat) %*% solve(RR %*% varb %*% t(RR)) %*% (RR %*% bhat)
    pvLLR <- cdf_chisq(as.numeric(LLR), restr)
    
    # ==========================================
    #             TVP and Optimal tests
    # ------------------------------------------
    t2v <- mround(n_obs * 0.15):mround(n_obs * 0.85)
    Andrews <- numeric(length(t2v))
    AP <- 0
    Nyb <- 0
    LLR7v <- numeric(length(t2v))
    AP0 <- 0
    AP00 <- 0
    Nyb0 <- 0

    cat(sprintf("Country %d: Computing TVP tests (%d break points)...\n", ii, length(t2v)))

    for (idx in seq_along(t2v)) {
      
      t2 <- t2v[idx]

      chow <- chow_gmm(e1, p1restr, p1unr, p1plusc, t2, heter)
      Andrews[idx] <- chow
      AP <- AP + exp(0.5 * chow)
      Nyb <- Nyb + chow * (t2 / n_obs) * (1 - t2 / n_obs)

      LLR7 <- chow_gmm_star(e1, p1restr, p1unr, p1plusc, t2, heter)
      LLR7v[idx] <- LLR7
      AP0 <- AP0 + exp(0.5 * LLR7)
      AP00 <- AP00 + LLR7
      Nyb0 <- Nyb0 + LLR7 * (t2 / n_obs) * (1 - t2 / n_obs)
      
    }

    # TVP test statistics
    SupLR <- max(Andrews)
    ExpW <- log((1 / (0.85 - 0.15)) * AP / n_obs)
    Nyblom_val <- (1 / (0.85 - 0.15)) * Nyb / n_obs

    # Optimal test statistics
    SupLRopt <- max(LLR7v)
    ExpWopt <- log((1 / (0.85 - 0.15)) * AP0 / n_obs)
    MeanWopt <- (1 / (0.85 - 0.15)) * AP00 / n_obs
    Nyblomopt <- nyblom_star(e1, p1restr, p1unr, p1plusc, heter, rep(0, ncol(p1restr)))

    # P-values for TVP tests
    pvSupLR <- pv_calc(SupLR, pvqlrsb, restr)
    pvExpW <- pv_calc(ExpW, pvapisb, restr)
    pvNyblom <- pv_calc(Nyblom_val, pvnybsb, restr)

    # P-values for optimal tests
    pvSupLRopt <- pv_calc(SupLRopt, pvqlropt, restr)
    pvExpWopt <- pv_calc(ExpWopt, pvapiopt, restr)
    pvMeanWopt <- pv_calc(MeanWopt, pvap0opt, restr)
    pvNyblomopt <- pv_calc(Nyblomopt, pvnybopt, restr)
    
    # ==========================================
    #              Out-of-sample tests
    # ------------------------------------------
    R <- mround(0.5 * n_obs)
    Pred <- n_obs - R

    yo_recur <- numeric(n_obs)
    yo_split <- numeric(n_obs)
    yo_roll <- numeric(n_obs)
    rw <- numeric(n_obs)
    true_val <- numeric(n_obs)

    cat(sprintf("Country %d: Computing OOS forecasts (%d periods)...\n", ii, Pred))

    for (j in 1:Pred) {
      
      # Recursive (expanding window)
      x_recur <- cbind(1, p1[1:(R + j - 1), , drop = FALSE])
      b_recur <- solve(crossprod(x_recur)) %*% crossprod(x_recur, e1[1:(R + j - 1)])
      yo_recur[R + j] <- c(1, p1[R + j,]) %*% b_recur

      # Split (fixed window)
      x_split <- cbind(1, p1[1:R, , drop = FALSE])
      b_split <- solve(crossprod(x_split)) %*% crossprod(x_split, e1[1:R])
      yo_split[R + j] <- c(1, p1[R + j,]) %*% b_split

      # Rolling window
      x_roll <- cbind(1, p1[j:(R + j - 1), , drop = FALSE])
      b_roll <- solve(crossprod(x_roll)) %*% crossprod(x_roll, e1[j:(R + j - 1)])
      yo_roll[R + j] <- c(1, p1[R + j,]) %*% b_roll

      rw[R + j] <- 0
      true_val[R + j] <- e1[R + j]
      
    }

    # Extract OOS period
    yo_r <- yo_recur[(R + 1):n_obs]
    yo_s <- yo_split[(R + 1):n_obs]
    yo_rl <- yo_roll[(R + 1):n_obs]
    rw_oos <- rw[(R + 1):n_obs]
    true_oos <- true_val[(R + 1):n_obs]

    # --- DM and ENC-NEW for RECURSIVE ---
    u1_r <- (true_oos - yo_r) ^ 2
    u2 <- (true_oos - rw_oos) ^ 2
    f_r <- u1_r - u2
    varf_r <- var(f_r) / Pred
    DMrecurs <- mean(f_r) / sqrt(varf_r)
    pvDMrecurs <- cdf_chisq(DMrecurs ^ 2, 1)
    ENCrecurs <- Pred * mean(u2 ^ 2 - u1_r * u2) / mean((u1_r - mean(u1_r)) ^ 2)
    pvENCrecurs <- 1
    
    if (ENCrecurs > oos_cv(restr, 0.10, 2, 1)) pvENCrecurs <- 0.10
    if (ENCrecurs > oos_cv(restr, 0.05, 2, 1)) pvENCrecurs <- 0.05
    if (ENCrecurs > oos_cv(restr, 0.01, 2, 1)) pvENCrecurs <- 0.01

    # --- DM and ENC-NEW for ROLLING ---
    u1_rl <- (true_oos - yo_rl) ^ 2
    f_rl <- u1_rl - u2
    varf_rl <- var(f_rl) / Pred
    DMroll <- mean(f_rl) / sqrt(varf_rl)
    pvDMroll <- cdf_chisq(DMroll ^ 2, 1)
    ENCroll <- Pred * mean(u2 ^ 2 - u1_rl * u2) / mean((u1_rl - mean(u1_rl)) ^ 2)
    pvENCroll <- 1
    
    if (ENCroll > oos_cv(restr, 0.10, 2, 2)) pvENCroll <- 0.10
    if (ENCroll > oos_cv(restr, 0.05, 2, 2)) pvENCroll <- 0.05
    if (ENCroll > oos_cv(restr, 0.01, 2, 2)) pvENCroll <- 0.01

    # --- DM and ENC-NEW for SPLIT ---
    u1_s <- (true_oos - yo_s) ^ 2
    f_s <- u1_s - u2
    varf_s <- var(f_s) / Pred
    DMsplit <- mean(f_s) / sqrt(varf_s)
    pvDMsplit <- cdf_chisq(DMsplit ^ 2, 1)
    ENCsplit <- Pred * mean(u2 ^ 2 - u1_s * u2) / mean((u1_s - mean(u1_s)) ^ 2)
    pvENCsplit <- 1
    
    if (ENCsplit > oos_cv(restr, 0.10, 2, 3)) pvENCsplit <- 0.10
    if (ENCsplit > oos_cv(restr, 0.05, 2, 3)) pvENCsplit <- 0.05
    if (ENCsplit > oos_cv(restr, 0.01, 2, 3)) pvENCsplit <- 0.01

    # Store results
    full_results[1, ii]  <- round(as.numeric(LLR), 2)
    full_results[2, ii]  <- round(pvLLR, 2)
    full_results[3, ii]  <- round(SupLR, 2)
    full_results[4, ii]  <- round(pvSupLR, 2)
    full_results[5, ii]  <- round(ExpW, 2)
    full_results[6, ii]  <- round(pvExpW, 2)
    full_results[7, ii]  <- round(Nyblom_val, 2)
    full_results[8, ii]  <- round(pvNyblom, 2)
    full_results[9, ii]  <- round(ExpWopt, 2)
    full_results[10, ii] <- round(pvExpWopt, 2)
    full_results[11, ii] <- round(MeanWopt, 2)
    full_results[12, ii] <- round(pvMeanWopt, 2)
    full_results[13, ii] <- round(Nyblomopt, 2)
    full_results[14, ii] <- round(pvNyblomopt, 2)
    full_results[15, ii] <- round(SupLRopt, 2)
    full_results[16, ii] <- round(pvSupLRopt, 2)
    full_results[17, ii] <- round(DMsplit, 2)
    full_results[18, ii] <- round(pvDMsplit, 2)
    full_results[19, ii] <- round(DMrecurs, 2)
    full_results[20, ii] <- round(pvDMrecurs, 2)
    full_results[21, ii] <- round(DMroll, 2)
    full_results[22, ii] <- round(pvDMroll, 2)
    full_results[23, ii] <- ENCsplit
    full_results[24, ii] <- ENCrecurs
    full_results[25, ii] <- ENCroll

    cat(sprintf("Country %d done.\n\n", ii))
    
  }

  colnames(full_results) <- countries

  # Format ENC-NEW significance markers
  enc_format <- function(val, pv) {
    
    if (is.na(val)) return(NA_character_)
    marker <- ""
    
    if (!is.na(pv) && pv <= 0.01) marker <- "^1"
    
    else if (!is.na(pv) && pv <= 0.05) marker <- "^5"
    else if (!is.na(pv) && pv <= 0.10) marker <- "^{10}"
    
    sprintf("%.2f%s", val, marker)
    
  }

  # Build row names
  rn <- c(
    "LR_T", "(p-value)",
    "QLR_T", "(p-value)", "Exp-W_T", "(p-value)", "Nyblom_T", "(p-value)",
    "Exp-W*_T", "(p-value)", "Mean-W*_T", "(p-value)",
    "Nyblom*_T", "(p-value)", "QLR*_T", "(p-value)",
    "DM_T split", "(p-value)", "DM_T recur", "(p-value)",
    "DM_T roll", "(p-value)",
    "ENC_T split", "ENC_T recur", "ENC_T roll"
  )
  rownames(full_results) <- rn

  return(full_results)
  
}

# ==========================================
#             Run for AR(1), p = 2
# ------------------------------------------
cat("========================================\n")
cat("         AR(1) model, p = 2\n")
cat("========================================\n\n")
res_ar1 <- run_tests(71, data_raw)

# ==========================================
#             Run for AR(1), p = 3
# ------------------------------------------
cat("========================================\n")
cat("         AR(2) model, p = 3\n")
cat("========================================\n\n")
res_ar2 <- run_tests(72, data_raw)

# ==========================================
#                Print Table 1
# ------------------------------------------
cat("\n")
cat("================================================================\n")
cat("          TABLE 1. AR(1) and AR(2) models\n")
cat("================================================================\n\n")

# Helper: format ENC-NEW with significance superscripts
format_enc <- function(val, restr, m2_type) {
  
  if (is.na(val)) return("    NA")
  
  pv <- 1
  
  if (val > oos_cv(restr, 0.10, 2, m2_type)) pv <- 0.10
  if (val > oos_cv(restr, 0.05, 2, m2_type)) pv <- 0.05
  if (val > oos_cv(restr, 0.01, 2, m2_type)) pv <- 0.01

  marker <- ""
  
  if (pv <= 0.01) marker <- "***"
  
  else if (pv <= 0.05) marker <- "** "
  else if (pv <= 0.10) marker <- "*  "
  else marker <- "   "

  sprintf("%7.2f%s", val, marker)
  
}

print_table <- function(res, restr, label) {
  
  cat(sprintf("                     %s\n", label))
  cat(sprintf("%-16s %8s %8s %8s %8s %8s\n",
              "", "Can.", "Fr.", "Ger.", "It.", "Jap."))
  cat(paste(rep("-", 60), collapse = ""), "\n")

  # LR_T
  cat(sprintf("%-16s", "LR_T"))
  cat(sprintf(" %8.2f", res[1, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[2, ])))
  cat("\n")

  # TVP tests header
  cat(paste(rep(" ", 25), collapse = ""), "TVP tests\n")

  # QLR_T
  cat(sprintf("%-16s", "QLR_T"))
  cat(sprintf(" %8.2f", res[3, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[4, ])))
  cat("\n")

  # Exp-W_T
  cat(sprintf("%-16s", "Exp-W_T"))
  cat(sprintf(" %8.2f", res[5, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[6, ])))
  cat("\n")

  # Nyblom_T
  cat(sprintf("%-16s", "Nyblom_T"))
  cat(sprintf(" %8.2f", res[7, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[8, ])))
  cat("\n")

  # Optimal tests header
  cat(paste(rep(" ", 22), collapse = ""), "Optimal tests\n")

  # Exp-W*_T
  cat(sprintf("%-16s", "Exp-W*_T"))
  cat(sprintf(" %8.2f", res[9, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[10, ])))
  cat("\n")

  # Mean-W*_T
  cat(sprintf("%-16s", "Mean-W*_T"))
  cat(sprintf(" %8.2f", res[11, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[12, ])))
  cat("\n")

  # Nyblom*_T
  cat(sprintf("%-16s", "Nyblom*_T"))
  cat(sprintf(" %8.2f", res[13, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[14, ])))
  cat("\n")

  # QLR*_T
  cat(sprintf("%-16s", "QLR*_T"))
  cat(sprintf(" %8.2f", res[15, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[16, ])))
  cat("\n")

  # OOS tests header
  cat(paste(rep(" ", 24), collapse = ""), "OOS tests\n")

  # DM_T split
  cat(sprintf("%-16s", "DM_T split"))
  cat(sprintf(" %8.2f", res[17, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[18, ])))
  cat("\n")

  # DM_T recur
  cat(sprintf("%-16s", "DM_T recur"))
  cat(sprintf(" %8.2f", res[19, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[20, ])))
  cat("\n")

  # DM_T roll
  cat(sprintf("%-16s", "DM_T roll"))
  cat(sprintf(" %8.2f", res[21, ]))
  cat("\n")
  cat(sprintf("%-16s", ""))
  cat(sprintf(" %8s", sprintf("(%.2f)", res[22, ])))
  cat("\n")

  # ENC_T split (with significance markers)
  cat(sprintf("%-16s", "ENC_T split"))
  for (j in 1:5) cat(format_enc(res[23, j], restr, 3))
  cat("\n")

  # ENC_T recur
  cat(sprintf("%-16s", "ENC_T recur"))
  for (j in 1:5) cat(format_enc(res[24, j], restr, 1))
  cat("\n")

  # ENC_T roll
  cat(sprintf("%-16s", "ENC_T roll"))
  for (j in 1:5) cat(format_enc(res[25, j], restr, 2))
  cat("\n")
  
}

print_table(res_ar1, 2, "AR(1), p = 2")
cat("\n")
print_table(res_ar2, 3, "AR(2), p = 3")

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "secs")
cat(sprintf("\nElapsed time: %.1f seconds\n", as.numeric(elapsed_time)))
