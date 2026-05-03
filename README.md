# Are Exchange Rates Really Random Walks? Some Evidence Robust to Parameter Instability

Replication of Rossi (2006), "Are exchange rates really random walks? Some evidence robust to parameter instability," *Macroeconomic Dynamics*, Vol. 10, pp. 20-38.

## What is replicated

| Output | Script | Description |
|--------|--------|-------------|
| Table 1 | `table_1.R` | LR, TVP, optimal, and OOS tests for AR(1) and AR(2) models across 5 exchange rates |

Five bilateral nominal exchange rates versus the U.S. dollar: Canada, France, Germany, Italy, and Japan. AR(1) with *p* = 2 restrictions and AR(2) with *p* = 3 restrictions, tested against a driftless random walk.

## How to run

All scripts are R and must be run from the **repo root**.

```r
source("table_1.R")   # Table 1: AR(1) and AR(2) models (~2 sec)
```

### Dependencies

- `expm`

## Methods

### Likelihood ratio test (LR<sub>T</sub>)

Wald-type test for &beta; = 0 using GMM with heteroskedasticity-robust (HC) covariance matrices. *p*-values from the &chi;<sup>2</sup> distribution with *p* degrees of freedom.

### TVP tests (QLR<sub>T</sub>, Exp-W<sub>T</sub>, Nyblom<sub>T</sub>)

Quandt (1960) likelihood ratio, Andrews and Ploberger (1994) exponential Wald, and Nyblom (1989) tests for time-varying parameters. Structural break scanned over [0.15*T*, 0.85*T*]. *p*-values from simulated critical value tables (5,000 Monte Carlo replications).

### Optimal tests (Exp-W<sup>\*</sup><sub>T</sub>, Mean-W<sup>\*</sup><sub>T</sub>, Nyblom<sup>\*</sup><sub>T</sub>, QLR<sup>\*</sup><sub>T</sub>)

Rossi (2005) joint tests for both parameter instability and &beta; = 0. *p*-values from simulated critical value tables.

### Out-of-sample tests (DM<sub>T</sub>, ENC<sub>T</sub>)

Diebold and Mariano (1995) equal predictive accuracy test and Clark and McCracken (2001) ENC-NEW encompassing test for nested models, under three estimation schemes: split (fixed window), recursive (expanding window), and rolling window.

| Metric | Description |
|--------|-------------|
| DM<sub>T</sub> | Diebold-Mariano test; *p*-values from &chi;<sup>2</sup>(1) |
| ENC<sub>T</sub> | Clark-McCracken ENC-NEW; significance from nonstandard critical value tables |

### Key parameters

| Parameter | Value | Meaning |
|-----------|-------|---------|
| Sample | 1973:03-1998:12 | 310 monthly Datastream observations |
| Sample split | *R* = *T*/2 | First 50% for estimation |
| Covariance | HC-robust | Conditionally heteroskedastic |
| Break trimming | [0.15, 0.85] | Andrews (1993) trimming |
| AR(1) restrictions | *p* = 2 | Constant + 1 lag |
| AR(2) restrictions | *p* = 3 | Constant + 2 lags |

## Repository structure

```
data/
  PredictorData1998.txt             <- Datastream monthly series (310 obs x 29 vars)
data-raw/
  ooscvdat.txt                      <- Clark-McCracken ENC-NEW critical values
  pv*.txt                           <- simulated p-value tables for TVP and optimal tests
R/
  gmm-beta.R                        <- gmm_beta(): GMM coefficient estimation (Hayashi, ch. 3)
  gmm-var.R                         <- gmm_var(): GMM asymptotic variance
  gmm-res.R                         <- gmm_res(): GMM residuals
  chow-gmm.R                        <- chow_gmm(): Chow structural break test via GMM
  chow-gmm-star.R                   <- chow_gmm_star(): optimal Chow* test (Rossi, 2005)
  nyblom-star.R                     <- nyblom_star(): optimal Nyblom test
  pv-calc.R                         <- pv_calc(), cdf_chisq(), oos_cv(): p-value interpolation
  dm-test.R                         <- dm_test(): Diebold-Mariano test
  enc-new.R                         <- enc_new(): Clark-McCracken ENC-NEW test
refs/
  *.pdf                             <- paper (Rossi, 2006)
table_1.R                           <- main replication script
```

## References

- Andrews, D. W. K. (1993). Tests for parameter instability and structural change with unknown change point. *Econometrica*, 61(4), 821-856.
- Andrews, D. W. K., and Ploberger, W. (1994). Optimal tests when a nuisance parameter is present only under the alternative. *Econometrica*, 62(6), 1383-1414.
- Clark, T. E., and McCracken, M. W. (2001). Tests of equal forecast accuracy and encompassing for nested models. *Journal of Econometrics*, 105(1), 85-110.
- Diebold, F. X., and Mariano, R. S. (1995). Comparing predictive accuracy. *Journal of Business & Economic Statistics*, 13(3), 253-263.
- Meese, R. A., and Rogoff, K. (1983). Empirical exchange rate models of the seventies: Do they fit out of sample? *Journal of International Economics*, 14(1-2), 3-24.
- Nyblom, J. (1989). Testing for the constancy of parameters over time. *Journal of the American Statistical Association*, 84(405), 223-230.
- Rossi, B. (2005). Optimal tests for nested model selection with underlying parameter instability. *Econometric Theory*, 21(5), 962-990.
- Rossi, B. (2006). Are exchange rates really random walks? Some evidence robust to parameter instability. *Macroeconomic Dynamics*, 10(1), 20-38.
