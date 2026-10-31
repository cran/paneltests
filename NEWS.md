# paneltests 1.0.8

## Breaking changes

* `xtpedroni()`: p-values are one-sided in the tail stated by Pedroni (1999, p. 668): `1 - pnorm(z)` for the panel variance ratio and `pnorm(z)` for the other six statistics. Up to 1.0.7 all seven were `1 - pnorm(abs(z))`. The new element `tail` records the tail of each statistic and the print method shows it.
* `xtkao()`: p-values are `pnorm(stat)` (left tail) for all five statistics. Up to 1.0.7 they were `1 - pnorm(abs(stat))`.
* `xtpedroni()` with one regressor and `trend = TRUE`: the standardization moments are Pedroni (2004, Corollary 1) case 3, mu = (17.86, -10.54, -2.29, -13.65, -2.53), nu = (101.68, 39.52, 0.66, 50.91, 0.56). Up to 1.0.7 the Pedroni (1999, Table 2) row for two regressors with trend was used. The internal "none" case holds Corollary 1 case 1. The warning attached to these cases is removed.
* `xtpedroni()`: Pedroni (1999, Table 2) cell for three regressors, heterogeneous intercepts, panel t mean is -2.576 (was -2.567); it enters the panel PP t and panel ADF t for `k = 3`, `trend = FALSE`.
* `xtkao()`: the ADF t statistic uses `s_v^2 = RSS / (N T)` with `T` the number of periods of the DF step, as in Kao (1999, p. 9-10). Up to 1.0.7 the residual sum of squares was divided by the number of ADF observations `N (T - p)`. On the Stata reference panel the standardized ADF t now differs from `xtcointtest kao` by -0.009 to -0.19 (was +0.02 to +0.14); the other four Kao statistics and all seven Pedroni statistics are unchanged and agree with Stata within 1e-6.
* `xtpretest()`, module HETEROGENEITY: the Swamy statistic is centred on the `1/s_ii` weighted pooled estimator of Swamy (1970, eq. 6.2); up to 1.0.7 it was centred on the unweighted fixed-effects estimator. The returned list gains `beta_hat`.
* `xtpretest()`: rows with missing values in the response or regressors stop the function with an error giving their number and positions; up to 1.0.7 the F statistics were `NA` but printed as `F = 0.0000, p = 1.0000 [accept]` with a recommendation.
* `xtpretest()`: regressor columns exactly collinear with the others are removed with a warning naming them and `K` is the rank of the pooled design (returned as `dropped`); up to 1.0.7 `K` counted the aliased column and the degrees of freedom and F statistics were wrong.
* `xtpretest()`: fewer than two panels, or a panel with fewer than `K + 2` observations, stop the function with an error naming the panels; up to 1.0.7 short panels were skipped in `S1` while the degrees of freedom kept `N`.
* `xtpretest()`, modules HSIAO and HETEROGENEITY: an individual regression that is rank deficient (a regressor constant within one panel) stops the module with a message naming the panel, recorded in `failed`; up to 1.0.7 `lm.fit` fitted the reduced-rank design silently (HSIAO) or the panel contributed 0 to the Swamy statistic.

## Other changes

* `xtpretest()`, `xtpedroni()`, `xtkao()`: formula terms are evaluated with `model.frame()` and `model.matrix()`, so `log(x)`, `I(x^2)` and interactions work; each model-matrix column is one regressor.
* `xtpedroni()`, `xtkao()`: missing values are reported with their count and row positions; a warning is issued when `N < 2` or when the panel has fewer than 10 periods; `xtkao()` warns with more than seven regressors.
* `xtpretest()`, module ROBUST: the printed output and the man page state the simulated size of the HC1 Wald slope test under a true null (33.4 percent for N = 10, T = 20, from 8000 replications; 64.8 percent for N = 25, T = 20; 15.4 percent for N = 10, T = 50; 11.8 percent for N = 25, T = 100, the last three from 500 replications each; nominal 5 percent) and refer to the classical F test or a bootstrap. The estimator is unchanged.
* Man pages: the labels F1, F2, F3 of `xtpretest()` are mapped to Hsiao (2014) eqs. (2.2.14), (2.2.15), (2.2.16), there named F3, F1, F4; the CD module is documented as using fixed-effects residuals with pairwise-complete overlaps; the Stata conventions kept in `xtpedroni()` (demeaned DF residuals, weighted s*^2 in the panel ADF t) and `xtkao()` (AR(1) residual innovations, T0 - 1 fixed-effects sample) are stated as package choices.
* Header of `xtpretest()` prints the installed package version.
* Printed reference lines of `xtcsdq()` and `xtqsh()` use "and" instead of "&" between author names.
* Arrows removed from printed messages and comments; the 1.0.5 NEWS entry rewritten without the HTML-escaped arrow.
* Tests: Pedroni (2004) case 1, 2, 3 moments and the Table 2 m = 3 cell pinned to the printed values; p-value tails checked against `pnorm()`; Kao ADF t rebuilt with `lm.fit()`; Swamy statistic rebuilt by hand with `lm()`; Hsiao F2 with a collinear column rebuilt with `lm()`; errors for missing values, N = 1, short panels and rank-deficient panels; formula terms. The test `.xmp_locf fills forward correctly` referenced the retired namespace `xtmispanel:::` and now uses `paneltests:::`. The Stata agreement test for the Kao ADF t allows the documented gap (0.02 on the tested specification).
* `inst/validation/`: derivation ledgers for `xtpretest()`, `xtpedroni()` and `xtkao()`.
* `inst/CITATION` version updated to 1.0.8.

# paneltests 1.0.7

* New functions `xtpedroni()` and `xtkao()`: the seven Pedroni (1999, 2004) and five Kao (1999) residual-based panel cointegration statistics, following the formulas of the Stata manual entry for `xtcointtest` and its sample conventions (DF and ADF regressions without constant, demeaned residuals for the contemporaneous and Bartlett long-run variances, sums over the differenced sample; Pedroni moments from Pedroni (2004) for one regressor and Pedroni (1999, Table 2) for two to seven regressors; one-sided p-values from the absolute value of the standardized statistic). With the same kernel bandwidth and ADF lag, all seven Pedroni statistics and four of the five Kao statistics reproduce Stata 19.5 to at least 1e-6 on a six-country panel (four specifications, bandwidths 0, 1 and 2). The Kao ADF t is larger in absolute value in Stata by 0.02 to 0.25; its internal convention could not be identified. Stata's default Newey-West (1994) per-panel bandwidth is not reproduced; pass `kernel_lags`.
* `xtpretest()`, module ROBUST: the interacted design matrix included the regressors together with the interactions of all N panels, which are exactly collinear (the N interactions of a regressor sum to the regressor itself). `solve(crossprod(X))` therefore failed with "system is computationally singular" on every data set, and because the error was not caught it aborted the remaining modules when `tests = "ALL"`. The module now uses deviations from the first panel (N-1 dummies and (N-1)K interactions), so the slope test has (N-1)K restrictions as in Hsiao's F2 and the intercept test N-1 restrictions. Aliased columns (for example a regressor constant within a panel) are dropped and reported in `robust$aliased`.
* Module ROBUST: the HC1 "meat" was computed as `crossprod(X * e^2)`, which is `sum e^4 x x'`; it is now `sum e^2 x x'`. The inverse uses a QR decomposition.
* Each module of `xtpretest()` runs inside `tryCatch()`; a failing module is skipped with a message and recorded in the new element `failed` instead of stopping the function.
* Documentation now states that the ROBUST module is an extension of the package and not part of Hsiao (2014), who only notes that the F-tests assume i.i.d. errors (3rd ed., p. 24).
* New tests for the ROBUST module (clean panel, agreement with a direct sandwich computation, aliased regressor, no abort of other modules).

# paneltests 1.0.6

* Corrected the DOIs of Galvao, Juhl, Montes-Rojas and Olmo (2017) (10.1093/jjfinec/nbx016), Little (1988) (10.1080/01621459.1988.10478722) and Swamy (1970) (10.2307/1913012).
* Demetrescu, Hosseinkouchack and Rodrigues (2023): title corrected to the published working-paper title and DOI corrected to 10.4419/96973210 (Ruhr Economic Papers 1041); the journal DOI previously given in DESCRIPTION pointed to an unrelated article and was removed.
* No changes to code.

# paneltests 1.0.5

- Resubmission addressing CRAN reviewer (Benjamin Altmann) feedback on 1.0.4:
  - Added DOI / ISBN links to all references in the Description field (Hsiao 2014, Swamy 1970, Pesaran 2004).
  - Removed the surrounding single quotes around function names in the Description text and added trailing parentheses, e.g. `'xtpretest'` became `xtpretest()`.
- Synchronised CITATION version with DESCRIPTION (1.0.5).
- No code or behavior changes.

# paneltests 1.0.4

- Resubmission to CRAN.
- Removed UTF-8 BOM (Byte Order Mark) from `inst/CITATION` and from `R/xtcsdq.R`, `R/xtmispanel.R`, `R/xtqsh.R`. The BOM in `inst/CITATION`was causing the CRAN incoming auto-check to fail with: `Invalid citation information in 'inst/CITATION': 1:1: unexpected input`.
- Synchronised the version string in `inst/CITATION` with `DESCRIPTION`.
- No code or behavior changes.

# paneltests 1.0.3

- Resubmission. Removed redundant "a suite of tools for" from the Description start. All references use `doi:10.xxx` format.

# paneltests 1.0.0

- Initial CRAN release.
- Consolidates four previously separate packages (xtpretest, xtmispanel, xtcsdq, xtqsh) into a single panel-data pre-testing suite.
- Implements:
  - Hsiao (2014) homogeneity F-tests (F1/F2/F3),
  - Swamy (1970) parameter heterogeneity test,
  - Pesaran (2004) cross-sectional dependence test (`xtpretest`),
  - Missing-data detection, mechanism testing, and imputation for unbalanced panels (`xtmispanel`),
  - Demetrescu, Hosseinkouchack and Rodrigues (2023) quantile-regression cross-sectional dependence tests (`xtcsdq`),
  - Galvao, Juhl, Montes-Rojas and Olmo (2017) panel quantile-regression slope homogeneity statistics (`xtqsh`).
