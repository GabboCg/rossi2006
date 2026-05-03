pv_calc <- function(osser, tavola, kk) {
  
  # Calculate p-values of TVP and Optimal tests via interpolation
  # osser: observed value of the statistic
  # tavola: table with critical values (34 x (1 + num_restr_cols))
  # kk: number of restrictions (column index)

  col_idx <- kk + 2  # column 1 is p-value grid, column kk+2 is the critical value
  numero <- tavola[, c(1, col_idx)]

  if (osser <= tavola[1, col_idx]) {
    
    pv <- 1
    
  } else if (osser >= tavola[34, col_idx]) {
    
    pv <- 0
    
  } else {
    
    rigal <- which(numero[, 2] <= osser)
    riga <- c(max(rigal), max(rigal) + 1)
    sel <- tavola[riga, c(1, col_idx)]
    pv <- sel[2, 1] + (sel[2, 2] - osser) * (sel[1, 1] - sel[2, 1]) / (sel[2, 2] - sel[1, 2])
    
  }

  return(pv)
  
}

cdf_chisq <- function(x, k) {
  
  # p-value from chi-squared distribution
  1 - pchisq(x, k)
  
}

oos_cv <- function(r, a, m1, m2) {
  
  # Critical values for OOS tests
  # r: number of restrictions
  # a: significance level (0.01, 0.05, 0.10)
  # m1: 1 = DM, 2 = ENC-NEW
  # m2: 1 = recursive, 2 = rolling, 3 = fixed/split

  ooscvdat <- as.matrix(read.table("data-raw/ooscvdat.txt"))

  if (m1 == 1) data <- ooscvdat[,1:9]
  if (m1 == 2) data <- ooscvdat[,10:18]

  if (m2 == 1) data1 <- data[,1:3]
  if (m2 == 2) data1 <- data[,4:6]
  if (m2 == 3) data1 <- data[,7:9]

  row <- data1[r,]

  if (a == 0.01) col <- 1
  if (a == 0.05) col <- 2
  if (a == 0.10) col <- 3

  return(row[col])
  
}
