# Derivation ledger: xtkao() (paneltests 1.0.8)

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
