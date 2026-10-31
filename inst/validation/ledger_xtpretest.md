# Derivation ledger: xtpretest() (paneltests 1.0.8)

Origin: audit ledger for 1.0.7 (2026-10-30), updated to the 1.0.8 code lines and verification status. Source: R/xtpretest.R. Paper texts: Hsiao (2014) 3rd ed., section 2.2 (pp. 20-22); Swamy (1970) Econometrica 38(2), eq. (6.2) p. 319 and eq. (2.11) p. 315; Pesaran (2004): no full text available.

## Input handling (new in 1.0.8)

| Item | Code file:line | Convention | Verification |
|---|---|---|---|
| formula terms | xtpretest.R:149-157 | model.frame + model.matrix (intercept column removed); each column one regressor | test "formula terms log() and I() are evaluated" |
| missing values | xtpretest.R:162-166 | stop with count and row positions (the tests define no treatment) | test "missing values are refused" |
| exactly collinear regressors | xtpretest.R:170-179 | QR rank of cbind(1, X); aliased columns dropped with a warning; K = rank; returned as `dropped` | test "exactly collinear regressors are dropped": F2 equals the F2 without the column and the lm() hand value (df (10, 72)) |
| N = 1 | xtpretest.R:189 | stop | test |
| panel with T_i < K + 2 | xtpretest.R:190-194 | stop naming the panels | test |
| regressor constant within one panel | xtpretest.R:306-308 (HSIAO), 486-488 (HETEROGENEITY) | module stops with the panel name, recorded in `failed`; ROBUST drops the aliased columns and reports them (xtpretest.R:415-421) | test "a regressor constant within one panel stops HSIAO" |

## Module HSIAO

| Paper equation / page | Expression | Code file:line | Convention | Source | Verification |
|---|---|---|---|---|---|
| Table 2.1, p. 22 | S1 = RSS of N individual OLS regressions, df N(T - K - 1) | xtpretest.R:302-311 | lm.fit per panel on cbind(1, X_i); rank checked | paper | audit 1e-10; collinearity test rebuilds S1 with lm() |
| (2.2.11), p. 20 | S2 = RSS of the within regression | xtpretest.R:314-318 | lm.fit on N dummies plus X | paper | audit 1e-10; test with lm() |
| (2.2.13), p. 20 | S3 = RSS of pooled OLS | xtpretest.R:321-325 | lm.fit on cbind(1, X) | paper | audit 1e-10 |
| (2.2.14), p. 21, Hsiao's F3 | overall homogeneity | xtpretest.R:329-331 | reported as F1; df (N - 1)(K + 1), df1 | paper (formula); label mapping now in the man page and printed output | audit 1e-10 |
| (2.2.15), p. 21, Hsiao's F1 | slope homogeneity | xtpretest.R:333-335 | reported as F2; df (N - 1)K, df1 | paper | audit 1e-10; test |
| (2.2.16), p. 21, Hsiao's F4 | intercept homogeneity | xtpretest.R:337-339 | reported as F3; df N - 1, N(T - 1) - K | paper | audit 1e-10 |
| p. 20 | exact F under i.i.d. normal errors | pf() upper tail | paper | audit 1e-10 |
| printing | NA statistics | xtpretest.R:362-365 | printed as "not computed (NA)" (1.0.7 printed F = 0, p = 1, accept) | package | fixed; unreachable after the input checks except through module failure |

## Module ROBUST (package extension)

| Item | Code file:line | Convention | Verification |
|---|---|---|---|
| interacted design, first panel as reference | xtpretest.R:387-401 | N - 1 dummies, (N - 1)K interactions | audit 1e-10; test vs direct sandwich |
| HC1 sandwich | xtpretest.R:446-457 | n/(n - k) sum e^2 x x', QR inverse | audit matched sandwich::vcovHC HC1 at 4e-16 |
| Wald/q to F(q, n - k) | xtpretest.R:460-469 | q = (N - 1)K and N - 1 | audit 1e-10 |
| size | man page and printed note (xtpretest.R:437-439) | simulated rejection of a true null at nominal 5 percent: 33.4 percent (N = 10, T = 20; R = 8000, s.e. 0.5 points; the audit's R = 500 value of 34.4 percent is within 2 s.e.), 64.8 (N = 25, T = 20), 15.4 (N = 10, T = 50), 11.8 (N = 25, T = 100) (last three from the audit's R = 500 grid); classical F 3.4 to 6.4 percent | documented in 1.0.8; estimator unchanged |

## Module HETEROGENEITY (Swamy 1970)

| Paper equation / page | Expression | Code file:line | Convention | Source | Verification |
|---|---|---|---|---|---|
| (6.2), p. 319 | H = sum_i (b_i - beta_hat)' X_i'X_i (b_i - beta_hat) / s_ii, beta_hat = [sum X_i'X_i/s_ii]^-1 sum (X_i'X_i/s_ii) b_i | xtpretest.R:483-500 | X_i within-panel demeaned (Frisch-Waugh: equals the inverse slope block of s_ii (X_i'X_i)^-1 with intercept); weighted centring | paper | FIXED in 1.0.8 (1.0.7 centred on the unweighted FE estimator; audit: 10.3546 vs 10.2835 on its test data; size under heteroskedastic panels 0.162 vs 0.062 at N = 5, T = 100); test rebuilds H and beta_hat with lm() at 1e-10 and checks W_i = solve(vcov(lm)[-1, -1]) |
| (2.11), p. 315 | s_ii | xtpretest.R:489 | RSS_i / (T_i - K - 1), intercept counted | package choice (Pesaran and Yamagata 2008 use the same); Swamy divides by T - Lambda | recorded; stated in the man page |
| p. 319 | chi-square with Lambda(N - 1) df | xtpretest.R:501-502 | K(N - 1) | paper | test |
| singular V_i | xtpretest.R:486-488 | stop (1.0.7 added 0 silently) | package | test (constant regressor in one panel) |

## Module CSD (Pesaran 2004; no full text)

| Item | Code file:line | Convention | Verification |
|---|---|---|---|
| CD = sqrt(2/(N(N - 1))) sum_{i<j} sqrt(T_ij) rho_ij | xtpretest.R:540-553 | pairwise-complete T_ij; residuals of the FE regression with common slopes (stated in the man page in 1.0.8) | balanced form rebuilt in the audit at 1e-10; not verifiable against the paper |
| p-value | xtpretest.R:554 | two-sided normal | package choice |
