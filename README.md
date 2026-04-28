# rossi2006

Replication of out-of-sample forecast comparison tests from **Table 1** (AR(1) model, recursive scheme) in:

> Rossi, B. (2006). Are Exchange Rates Really Random Walks? Some Evidence Robust to Parameter Instability. *Macroeconomic Dynamics*, 10(1), 20-38. DOI: [10.1017/S1365100506050085](https://doi.org/10.1017/S1365100506050085)

## Overview

The paper tests whether macroeconomic fundamentals can forecast nominal exchange rates better than a random walk, using tests that are robust to parameter instability. This replication focuses on the **out-of-sample tests** (Diebold-Mariano and ENC-NEW) for the AR(1) model using a recursive (expanding window) estimation scheme.

Five bilateral exchange rates against the U.S. dollar are considered: Canada, France, Germany, Italy, and Japan. The exchange rate growth is modeled as an AR(1) process and compared against a driftless random walk (forecast = 0).

## Results

Replication of the DM<sub>T</sub> recursive row from Table 1 (AR(1), *p* = 2):

|         | Can. | Fr.  | Ger. | It.  | Jap. |
|---------|------|------|------|------|------|
| DM      | 1.16 | 2.24 | 0.14 | 0.77 | 0.00 |
| p-value | 0.25 | 0.03 | 0.89 | 0.44 | 1.00 |
| ENC-NEW | -0.45| -1.19| 0.26 | 0.17 | 1.19 |

The Diebold-Mariano test statistics match the paper exactly. The ENC-NEW values differ because the paper uses heteroskedasticity-robust covariance estimates, while this implementation uses the standard Clark and McCracken (2001) formula.

## Methodology

- **Data**: Monthly Datastream series, March 1973 to December 1998 (310 observations)
- **Sample split**: First 50% of usable observations for initial estimation, remainder for out-of-sample evaluation (recursive/expanding window)
- **Restricted model (H₀)**: Driftless random walk — forecast of log exchange rate change is zero
- **Unrestricted model (H₁)**: AR(1) — `x_{1,t} = beta_1 + beta_2 * x_{1,t-1} + epsilon_t`
- **DM test**: Diebold and Mariano (1995) test for equal predictive accuracy, p-values from chi-squared distribution
- **ENC-NEW test**: Clark and McCracken (2001) encompassing test for nested models

## Repository Structure

```
data/             Goyal-Welch predictor dataset (XLS format)
R/                Helper functions
  dm-test.R         Diebold-Mariano test statistic
  enc-new.R         Clark-McCracken ENC-NEW test statistic
refs/             Paper (PDF)
table_1.R         Main replication script
```

## How to Run

1. Open `rossi2006.Rproj` in RStudio (or set working directory to the repo root)
2. Run `table_1.R`

### Dependencies

- `readxl`
- `dplyr`
- `janitor`

## References

- Clark, T. E., & McCracken, M. W. (2001). Tests of Equal Forecast Accuracy and Encompassing for Nested Models. *Journal of Econometrics*, 105(1), 85-110.
- Diebold, F. X., & Mariano, R. S. (1995). Comparing Predictive Accuracy. *Journal of Business & Economic Statistics*, 13(3), 253-263.
- Meese, R. A., & Rogoff, K. (1983a). Empirical Exchange Rate Models of the Seventies: Do They Fit Out of Sample? *Journal of International Economics*, 14(1-2), 3-24.
- Rossi, B. (2006). Are Exchange Rates Really Random Walks? Some Evidence Robust to Parameter Instability. *Macroeconomic Dynamics*, 10(1), 20-38.
