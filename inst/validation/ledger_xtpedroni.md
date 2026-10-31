# Derivation ledger: xtpedroni() (paneltests 1.0.8)

Origin: audit ledger for 1.0.7 (2026-10-30), updated to the 1.0.8 code lines and verification status. Source: R/xtpedroni.R. Paper texts: Pedroni (1999) OBES 61(S1), Table 1 and steps 1-5 (pp. 658-662), Table 2 (p. 666), p. 668 (tails); Pedroni (2004) ET 20(3), Corollary 1 and Table 1 (pp. 606-607). "Hand rebuild" means an independent lm()/lm.fit() computation in the audit scripts or in tests/testthat/test-xtpedroni-xtkao.R.

## Input handling (new in 1.0.8)

| Item | Code file:line | Convention | Verification |
|---|---|---|---|
| formula terms | xtpedroni.R:128 (.coint_frame) | model.frame + model.matrix; each column one regressor; log(), I() work | test "formula terms log() and I()" (equal to pre-transformed columns) |
| missing values | .coint_frame | stop with count and row positions | test "missing values are refused" |
| k > 7 | xtpedroni.R:130 | stop, no extrapolation of Table 2 | test "at most 7" |
| N < 2, T0 < 10 | xtpedroni.R:135-136 | warning | tests |

## Residuals and nuisance parameters

| Paper equation / page | Mathematical expression | Code file:line | Convention adopted | Source | Verification |
|---|---|---|---|---|---|
| Pedroni 1999 step 1, p. 659 | e_it from the panel-specific level regression with intercept (and trend), all T0 periods | xtpedroni.R:148-150 | lm.fit with constant, optional trend seq_len(T0) | paper | audit hand rebuild 1e-10; unchanged in 1.0.8 |
| step 4(a), p. 659 | DF regression without constant; s_i^2, sigma_i^2, lambda_i | xtpedroni.R:155-156 | residuals demeaned before the variances (Stata; stated in the man page as the package's choice); divisor T = T0 - 1 | paper (form), Stata (demeaning) | audit hand rebuild 1e-10; unchanged |
| Table 1 footnote | Bartlett kernel, autocovariances / T | xtpedroni.R:139-141 | bart() | paper | audit 1e-10; unchanged |
| steps 2-3, pp. 658-659 | L11i^2 conditional long-run variance of dy given dx | xtpedroni.R:157-159 | Bartlett long-run covariance of demeaned first differences | paper (projection), Stata (demeaning) | audit 1e-10; unchanged |
| step 4(b), p. 662 | ADF regression with p lags, no constant; s_i*^2 | xtpedroni.R:173-177 | common lag, reduced sample (p + 2):T0, demeaned residuals | Stata | audit 1e-10; unchanged |
| bandwidth | per-panel k_i allowed by the paper | xtpedroni.R:137-138 | one integer for all panels; default floor(4 (T/100)^(2/9)) | package choice | unchanged |
| ADF lag selection | AIC/BIC on the pooled regression | xtpedroni.R:162-170 | common lag 1..lags | Stata | test "lag selection"; unchanged |

## The seven statistics (Pedroni 1999, Table 1)

| Paper | Expression | Code file:line | Convention | Source | Verification |
|---|---|---|---|---|---|
| no. 1 | T^2 N^{3/2} Z_v | xtpedroni.R:181 | as written | paper | audit 1e-10; Stata 1e-6 on ch6 (unchanged) |
| no. 2 | T sqrt(N) Z_rho | xtpedroni.R:182 | lambda summed over T | paper | audit 1e-10; Stata 1e-6 |
| no. 3 | Z_t, sigma~^2 = N^-1 sum L11i^-2 sigma_i^2 | xtpedroni.R:183 | as written | paper | audit 1e-10; Stata 1e-6 |
| no. 4 | Z_t*, s~^2 unweighted in the paper | xtpedroni.R:184 | s*^2 weighted by L11i^-2 (Stata) | Stata, stated in the man page | Stata 1e-6 (ch6 test); paper reading differs, recorded |
| no. 5-7 | group rho, group PP t, group ADF t | xtpedroni.R:185-187 | as written | paper | audit 1e-10; Stata 1e-6 |
| eq. (2), p. 665 | z = (chi - mu sqrt N)/sqrt(nu) | xtpedroni.R:189 | as written | paper | test with case 3 moments |
| p. 668 | panel v: right tail; other six: left tail | xtpedroni.R:196-197 | p = 1 - Phi(z) for panel v, Phi(z) otherwise; element `tail` | paper | FIXED in 1.0.8 (was 1 - Phi(abs z)); test "p-values: right tail ... left tail" |

## Moment tables

| Paper table | Code file:line | Status |
|---|---|---|
| Pedroni 2004 Corollary 1 case 1 (standard) | xtpedroni.R:220, 224 | FIXED in 1.0.8: (4.00, -2.77, -1.01, -6.84, -1.39), (27.81, 24.91, 1.50, 26.78, 0.78); unreachable from the exported function; test pins the values |
| Pedroni 2004 Corollary 1 case 2 (fixed effects) | xtpedroni.R:219, 223 | matched to printed precision; test pins the values |
| Pedroni 2004 Corollary 1 case 3 (fixed effects and trends) | xtpedroni.R:221, 225 | FIXED in 1.0.8: (17.86, -10.54, -2.29, -13.65, -2.53), (101.68, 39.52, 0.66, 50.91, 0.56); the audit's T = 1000 simulation gave (18.49, -10.90, -2.34, -13.62, -2.54) with s.e. (0.30, 0.18, 0.02, 0.20, 0.02), 0.2 to 2.2 s.e. from these values and 9 to 19 s.e. from the 1.0.7 row; test pins the values and the standardization |
| Pedroni 1999 Table 2, standard case m = 2..7 | xtpedroni.R:232-246 col 1 | 60 of 60 matched |
| Pedroni 1999 Table 2, heterogeneous intercepts m = 2..7 | col 2 | FIXED in 1.0.8: m = 3 panel t mean -2.576 (was -2.567); 60 of 60 matched; test pins the m = 3 row |
| Pedroni 1999 Table 2, trends and intercepts m = 2..7 | col 3 | 60 of 60 matched |

## Stata reference (ch6 panel, N = 6, T0 = 17, fixed bandwidths 0, 1, 2, ADF lag 1)

All seven standardized statistics, four specifications with k = 4, 6, 5, 4: max |R - Stata| = 1.2e-6 in 1.0.8 (unchanged; tests/testthat/test-xtpedroni-xtkao.R, gated by PANELTESTS_CH6_DATA). The p-values are not in the reference file and follow the paper.
