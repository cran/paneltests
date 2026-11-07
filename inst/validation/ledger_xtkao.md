# Derivation ledger: xtkao() (paneltests 1.0.9)

Origin: audit ledger for 1.0.7 (2026-10-30), updated to the 1.0.8 code lines and verification status. Source: R/xtkao.R. Paper: Kao (1999) J. Econometrics 90, eqs (1)-(5) pp. 3-4, section 3.1 pp. 7-9, section 3.2 pp. 9-11, eqs (10)-(12) p. 11.

## Input handling (new in 1.0.8)

| Item | Code file:line | Convention | Verification |
|---|---|---|---|
| formula terms | xtkao.R:103 (.coint_frame in xtpedroni.R) | model.frame + model.matrix | test "formula terms log() and I()" |
| missing values | .coint_frame | stop with count and row positions | test |
| k > 7, N < 2, T0 < 10 | xtkao.R:109-111 | warnings (Kao tabulates finite-sample sizes for one regressor only) | tests |

| Paper equation / page | Mathematical expression | Code file:line | Convention adopted | Source | Verification |
|---|---|---|---|---|---|
| eq. (4), p. 3 | LSDV / within regression, common slope | xtkao.R:130-134 | on the t >= 2 sample (one period lost to the AR(1) step) | Stata, stated in the man page | audit 1e-10; unchanged |
| eqs (7)-(8), p. 7 | pooled DF regression without constant, rho_hat | xtkao.R:137 | T = T0 - 2 periods per panel | paper (form), Stata (sample) | audit 1e-10; unchanged |
| pp. 7-8 | t_rho, s_e^2 = (1/NT) sum (e_it - rho e_{i,t-1})^2 | xtkao.R:138 | RSS / n with n = N T | paper | audit 1e-10; unchanged |
| Remark 6, p. 8 | DF_rho*, DF_t* | xtkao.R:159-161 | as written | paper | audit 1e-10; Stata 4.6e-7 |
| Corollary 2, p. 9 | unadjusted DF_rho, DF_t | xtkao.R:163-164 | as written | paper | audit 1e-10; Stata 4.6e-7 |
| eq. (9), p. 9 | pooled ADF regression with p lags, no constant | xtkao.R:141-147 | common lag; sample (p + 2):T per panel | paper (form), Stata (common lag) | audit 1e-10; unchanged |
| pp. 9-10 | s_v^2 = (1/NT) sum sum v_itp^2; t_ADF = (rho - 1)[sum e_i'Q_i e_i]^{1/2}/s_v | xtkao.R:149-151 | s_v^2 = RSS / (N T), T = number of DF periods (reported as T) | paper | FIXED in 1.0.8 (1.0.7 divided by the number of ADF observations N(T - p)); test rebuilds t_ADF with lm.fit(). Stata gap on ch6 (standardized ADF, R minus Stata): -0.0088 (set 1, bw 0/1/2), -0.0099 (set 2), -0.185 to -0.188 (set 3), -0.157 to -0.159 (set 4); 1.0.7 gaps were +0.14, +0.14, +0.02, +0.04. Stata's convention is not the OLS t of eq. (9) with any fixed divisor on this sample; recorded as a known gap in the man page |
| Remark 8, p. 10 | standardized ADF | xtkao.R:162 | same add and den as DF_t* | paper | audit 1e-10 given t_adf |
| eqs (1)-(3), (10)-(12) | Sigma_hat, Omega_hat from w_it = (dy, dx) | xtkao.R:114-128 | residuals of pooled AR(1) regressions, within-demeaned per panel, averaged over panels | Stata, stated in the man page as the package's choice | Stata 4.6e-7 on ch6 confirms the convention; paper form differs (recorded) |
| Theorem 1 notation | conditional variances sigma_v^2, sigma_0v^2 | xtkao.R:127-128 | cond() | paper | audit 1e-10 |
| kernel | Bartlett, bandwidth l | xtkao.R:121-126 | one bandwidth for all panels | package choice | unchanged |
| p-values | all five statistics N(0,1) under H0 and diverge to minus infinity (Theorem 2, Corollary 2, Remarks 6, 8) | xtkao.R:167 | Phi(stat) | paper | FIXED in 1.0.8 (was 1 - Phi(abs stat)); test "p-values ... left tail for Kao" |

## Stata reference (ch6 panel, N = 6, T0 = 17, bandwidths 0, 1, 2, lags = 1)

Modified DF t, DF t, unadjusted modified DF t, unadjusted DF t: max |R - Stata| = 4.5e-7 (absolute) in 1.0.8, unchanged. ADF t: see the row above; the gated test allows the documented gap (0.02 on set 1, bandwidth 2).

## 1.0.9: convention = "kao" (default)

Oracle: Kao (1999, Table 4), the sizes at 5 percent produced by the author for the design of section 5.2 (pp. 15-16): a1 = 1, a2 = -1, beta = 2, alpha_i ~ U[0, 10], rho = 1, sigma = 1, pi = 0, theta = 0, two ADF lags, Bartlett window 5, critical value -1.645. Script: inst/validation/kao_table4.R.

| Paper equation / page | Mathematical expression | Convention adopted with "kao" | How it was settled |
|---|---|---|---|
| eq. (4), p. 3 | LSDV on t = 1, ..., T | all T0 periods | DF_rho column of Table 4 (see the T0 row below) |
| eqs (10)-(12), p. 11 | Sigma_hat, Omega_hat from w_it = (u_it, eps_it')' | first differences (dy_it, dx_it')', not demeaned | DF_rho* column: mean absolute difference from Table 4 of 0.004 without demeaning, 0.038 with it (15 cells, 5000 replications) |
| p. 7 | s_e^2 = (1/NT) sum_i sum_{t=2}^T (e_it - rho e_i,t-1)^2 | RSS / (N (T0 - 1)) | unadjusted DF_t depends on this divisor only: 0.006 with N (T0 - 1), 0.054 with N T0 (15 cells, 5000 replications) |
| p. 9 | s_v^2 = (1/NT) sum sum v_itp^2 | RSS / (N (T0 - 1)) | ADF column with common lag coefficients: 0.016 with N (T0 - 1), 0.069 with N T0, 0.093 with the number of ADF residuals (9 cells, 3000 replications) |
| eq. (9), p. 9 (Q_i notation) | lagged differences in the ADF regression | coefficients common to all panels | ADF column with panel-specific coefficients: 0.041 with N (T0 - 1), 0.091 with N T0, 0.068 with the number of ADF residuals (same 9 cells) |
| pp. 7-8, Remark 6 | sqrt(N) T (rho_hat - 1) | T = T0 | DF_rho column: 0.006 with T0, 0.053 with T0 - 1 (9 cells, 3000 replications); DF_rho*: 0.006 and 0.020 |
| Remark 8, p. 10 | numerator sqrt(6N) sigma_v / (2 sigma_0u) | sigma_0v, as in Remark 6 | under otherwise equal settings the ADF column is 0.094 from Table 4 with sigma_0v and 0.181 with sigma_0u (15 cells, 5000 replications); read as a typographical error |

Numerical checks of the 1.0.9 code:

* convention = "kao" against an independent implementation of the steps above: max absolute difference 2.3e-14 over 200 random panels (N = 2 to 12, T0 = 10 to 40, p = 1 to 3, bandwidth 0 to 5).
* convention = "stata" against 1.0.8 on the same 200 panels: identical (difference 0).
* convention = "kao" against Kao (1999, Table 4), 15 cells (T = 10, 25, 50; N = 2, 5, 10, 25, 50), 5000 replications each, mean absolute difference: DF_rho 0.006, DF_t 0.007, DF_rho* 0.002, DF_t* 0.009, ADF 0.009 (largest single-cell difference: ADF at T = 10, N = 2, 0.174 against 0.138).
* Stata reference panel (N = 6, T0 = 17): the agreement test runs with convention = "stata" and passes as in 1.0.8.

## 1.0.9: further oracles and stress tests (2026-11-04)

Kao (1999, Tables 5 and 6), script inst/validation/kao_tables_5_6.R, two ADF lags, Bartlett window 5. Mean absolute difference from the printed rejection rates (columns DF_rho, DF_t, DF_rho*, DF_t*, ADF):

| Table | Cells, replications | convention = "kao" | convention = "stata" |
|---|---|---|---|
| 5, power, rho = 0.9 | 8 cells (row T = 50, N = 2 excluded: it repeats the row T = 25, N = 10 digit for digit), 3000 | 0.054, 0.025, 0.023, 0.032, 0.019 | larger (all 9 cells: 0.149, 0.027, 0.151, 0.087, 0.148) |
| 6, size, N = T = 50, 27 designs (sigma, pi, theta) | 27, 1500 | 0.004, 0.005, 0.004, 0.005, 0.023 | 0.015, 0.004, 0.055, 0.034, 0.075 |

Largest differences with "kao": Table 5, DF_rho at T = 10, N = 50 (0.459 against 0.626); Table 6, ADF with pi = -0.8 and sigma = 0.25 (0.172 against 0.038 for theta = 0, 0.265 against 0.062 for theta = 0.5). For pi = 0 and pi = 0.8 all five columns of Table 6 are within 0.008 on average.

Stress grid without a printed oracle (lags chosen by AIC up to 2, default bandwidth, 1000 replications per design, 76 designs): N = 3, 6, 10; T0 = 15, 17, 20, 30; K = 1, 2, 4 regressors; endogeneity (correlation 0.5) and MA(1) errors (+0.5, -0.5) in the regressors; power at rho = 0.6, 0.8, 0.9. Results:

* No errors and no non-finite statistics in any run.
* Size at the nominal 5 percent level, N = 6, T0 = 17, K = 1: 0.13 to 0.19 ("kao"), 0.05 to 0.17 ("stata"); K = 4: 0.35 to 0.50 and 0.21 to 0.43. From T0 = 15 to T0 = 30 (N = 3 to 10, K = 1) the rejection rates of the unadjusted DF t, DF_t* and ADF fall and those of DF_rho and DF_rho* do not.
* Power at rho = 0.8 rises with N and with T0 for both conventions (no decrease larger than 0.02 along either dimension); it is at least as high with "kao" as with "stata" in 139 of 140 design-statistic cells (one cell lower by 0.001). The rejection rates are not size-adjusted.

Functional checks (inputs): duplicated (panel, time) pairs, gaps in a numeric or Date time index, a response without within-panel variation, collinear or within-panel constant regressors, kernel_lags negative, fractional or not smaller than the number of innovations, lags 0 or infeasible, and an exact ADF fit stop with a message. Invariance (to 1e-10): rescaling and sign changes of the variables, panel-specific constants ("kao"; not "stata"), row order, panel labels and order, time labels.
