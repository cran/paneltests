## Reference values: Stata 19.5 `xtcointtest pedroni` and `xtcointtest kao`
## on a six-country panel (2008-2024, 17 periods) with fixed Bartlett
## bandwidths. The data set is internal to the project and is not
## distributed with the package, so these tests run only when the file is
## available (environment variable PANELTESTS_CH6_DATA). The simulated-data
## tests below always run.

ch6 <- Sys.getenv("PANELTESTS_CH6_DATA", unset = "")

sim_panel <- function(N = 6, TT = 20, k = 2, seed = 1) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(N), function(i) {
    X <- matrix(cumsum(rnorm(TT * k)), TT, k); y <- i + X %*% rep(0.5, k) + rnorm(TT)
    cbind(data.frame(id = i, t = seq_len(TT), y = as.numeric(y)), setNames(as.data.frame(X), paste0("x", seq_len(k)))) }))
}

test_that("xtpedroni returns seven standardized statistics and sensible structure", {
  d <- sim_panel()
  r <- xtpedroni(d, y ~ x1 + x2, index = c("id", "t"), kernel_lags = 1)
  expect_s3_class(r, "xtpedroni")
  expect_length(r$stats, 7); expect_length(r$p, 7)
  expect_true(all(is.finite(r$stats)))
  expect_true(all(r$p >= 0 & r$p <= 1))
  expect_equal(r$k, 2L); expect_equal(r$N, 6L); expect_equal(r$T, 19L)
  ## tails of Pedroni (1999, p. 668): right for the variance ratio, left otherwise
  expect_equal(unname(r$p[1]), 1 - pnorm(unname(r$stats[1])))
  expect_equal(unname(r$p[2:7]), pnorm(unname(r$stats[2:7])))
  expect_equal(unname(r$tail), c("right", rep("left", 6)))
  expect_output(print(r), "Pedroni")
})

test_that("xtpedroni: moments depend on the number of regressors and the deterministic case", {
  m1 <- .pedroni_moments(1L, "intercept"); m2 <- .pedroni_moments(2L, "intercept"); m2t <- .pedroni_moments(2L, "trend")
  expect_equal(unname(m1$mu["panel_v"]), 8.62); expect_equal(unname(m1$v["group_t"]), 0.66)
  expect_equal(unname(m2$mu["panel_v"]), 11.754); expect_equal(unname(m2t$mu["panel_v"]), 21.162)
  expect_error(xtpedroni(sim_panel(k = 8), y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, index = c("id", "t")), "at most 7")
})

test_that("xtpedroni: k = 1 moments are Pedroni (2004, Corollary 1) cases 1, 2 and 3", {
  ## Pedroni (2004, Econometric Theory 20, p. 606), Corollary 1, printed
  ## values for (panel v, panel rho, panel t, group rho, group t); the ADF t
  ## statistics share the moments of the t statistics.
  m_none <- .pedroni_moments(1L, "none"); m_fe <- .pedroni_moments(1L, "intercept"); m_tr <- .pedroni_moments(1L, "trend")
  five <- c("panel_v", "panel_rho", "panel_t", "group_rho", "group_t")
  expect_equal(unname(m_none$mu[five]), c(4.00, -2.77, -1.01, -6.84, -1.39))
  expect_equal(unname(m_none$v[five]),  c(27.81, 24.91, 1.50, 26.78, 0.78))
  expect_equal(unname(m_fe$mu[five]),   c(8.62, -6.02, -1.73, -9.05, -2.03))
  expect_equal(unname(m_fe$v[five]),    c(60.75, 31.27, 0.93, 35.98, 0.66))
  expect_equal(unname(m_tr$mu[five]),   c(17.86, -10.54, -2.29, -13.65, -2.53))
  expect_equal(unname(m_tr$v[five]),    c(101.68, 39.52, 0.66, 50.91, 0.56))
  expect_equal(unname(m_tr$mu["panel_adf"]), -2.29); expect_equal(unname(m_tr$v["group_adf"]), 0.56)
  ## trend = TRUE with one regressor uses case 3, without a warning
  d <- sim_panel(k = 1)
  expect_silent(r <- xtpedroni(d, y ~ x1, index = c("id", "t"), trend = TRUE, kernel_lags = 1))
  expect_equal(unname(r$moments$mu["panel_v"]), 17.86)
  expect_equal(unname(r$z["panel_t"]), unname((r$raw["panel_t"] - (-2.29) * sqrt(6)) / sqrt(0.66)))
})

test_that("xtpedroni: Pedroni (1999) Table 2, intercept case, m = 3, panel t mean is -2.576", {
  m3 <- .pedroni_moments(3L, "intercept")
  expect_equal(unname(m3$mu), c(15.197, -13.256, -2.576, -2.576, -16.888, -2.827, -2.827))
  expect_equal(unname(m3$v),  c(151.094, 81.772, 0.923, 0.923, 67.123, 0.585, 0.585))
})

test_that("p-values: right tail for the Pedroni variance ratio, left tail otherwise, left tail for Kao", {
  d <- sim_panel(k = 1, TT = 25, seed = 9)
  r <- xtpedroni(d, y ~ x1, index = c("id", "t"), kernel_lags = 1)
  z <- unname(r$stats)
  expect_equal(unname(r$p), c(pnorm(z[1], lower.tail = FALSE), pnorm(z[2:7])))
  k <- xtkao(d, y ~ x1, index = c("id", "t"), kernel_lags = 1)
  expect_equal(unname(k$p), pnorm(unname(k$stats)))
  ## a wrong-signed statistic has a large p-value
  expect_true(all(k$p[k$stats > 0] > 0.5)); expect_true(all(r$p[-1][r$stats[-1] > 0] > 0.5))
})

test_that("xtkao: ADF t uses s_v^2 = RSS / (N T) with T the number of DF periods (Kao 1999, p. 9-10)", {
  d <- sim_panel(k = 1, TT = 20, seed = 4)
  k <- xtkao(d, y ~ x1, index = c("id", "t"), lags = 1, kernel_lags = 0)
  ## rebuild: AR(1) residual step drops t = 1; FE on t >= 2; pooled DF / ADF
  N <- 6; T0 <- 20
  dd <- d[d$t >= 2, ]
  dm <- function(v) v - ave(v, dd$id)
  e <- dm(dd$y) - dm(dd$x1) * sum(dm(dd$x1) * dm(dd$y)) / sum(dm(dd$x1)^2)
  E <- split(e, dd$id)
  ys <- c(); Xs <- NULL
  for (i in seq_len(N)) { z <- E[[i]]; dz <- diff(z); idx <- 3:length(z)
    Xs <- rbind(Xs, cbind(z[idx - 1], dz[idx - 2])); ys <- c(ys, z[idx]) }
  fit <- lm.fit(Xs, ys); rss <- sum(fit$residuals^2)
  TT <- T0 - 2                                   # periods of the DF step, reported as T
  t_hand <- (fit$coefficients[1] - 1) / sqrt(rss / (N * TT) * solve(crossprod(Xs))[1, 1])
  expect_equal(k$T, TT)
  expect_equal(k$t_adf, unname(t_hand), tolerance = 1e-10)
  ## not the OLS t with the number of ADF observations
  t_obs <- (fit$coefficients[1] - 1) / sqrt(rss / length(ys) * solve(crossprod(Xs))[1, 1])
  expect_false(isTRUE(all.equal(k$t_adf, unname(t_obs))))
})

test_that("formula terms log() and I() are evaluated in xtpedroni and xtkao", {
  d <- sim_panel(k = 1, seed = 2); d$x1 <- exp(d$x1 / 5)
  d2 <- d; d2$lx <- log(d$x1); d2$sx <- d$x1^2
  r1 <- xtpedroni(d, y ~ log(x1) + I(x1^2), index = c("id", "t"), kernel_lags = 1)
  r2 <- xtpedroni(d2, y ~ lx + sx, index = c("id", "t"), kernel_lags = 1)
  expect_equal(r1$stats, r2$stats); expect_equal(r1$k, 2L)
  k1 <- xtkao(d, y ~ log(x1), index = c("id", "t"), kernel_lags = 1)
  k2 <- xtkao(d2, y ~ lx, index = c("id", "t"), kernel_lags = 1)
  expect_equal(k1$stats, k2$stats)
})

test_that("missing values are refused with count and positions; small N, short T and many regressors warn", {
  d <- sim_panel(k = 1); d$x1[c(4, 30)] <- NA
  expect_error(xtpedroni(d, y ~ x1, index = c("id", "t")), "2 row\\(s\\).*rows 4, 30")
  expect_error(xtkao(d, y ~ x1, index = c("id", "t")), "2 row\\(s\\).*rows 4, 30")
  d1 <- sim_panel(k = 1, N = 1)
  expect_warning(xtpedroni(d1, y ~ x1, index = c("id", "t"), kernel_lags = 1), "N = 1")
  expect_warning(xtkao(d1, y ~ x1, index = c("id", "t"), kernel_lags = 1), "N = 1")
  ds <- sim_panel(k = 1, TT = 8)
  expect_warning(xtpedroni(ds, y ~ x1, index = c("id", "t"), kernel_lags = 0), "8 periods")
  expect_warning(xtkao(ds, y ~ x1, index = c("id", "t"), kernel_lags = 0), "8 periods")
  d8 <- sim_panel(k = 8, TT = 40)
  expect_warning(xtkao(d8, y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8, index = c("id", "t"), kernel_lags = 1), "8 regressors")
})

test_that("xtpedroni: lag selection chooses between 1 and the maximum", {
  d <- sim_panel(TT = 30, seed = 3)
  r <- xtpedroni(d, y ~ x1 + x2, index = c("id", "t"), lags = 3, lag_select = "aic", kernel_lags = 2)
  expect_true(r$adf_lags %in% 1:3)
  r2 <- xtpedroni(d, y ~ x1 + x2, index = c("id", "t"), lags = r$adf_lags, kernel_lags = 2)
  expect_equal(r$stats, r2$stats)
})

test_that("xtkao returns five statistics, and the unadjusted statistics do not depend on the kernel", {
  d <- sim_panel()
  k0 <- xtkao(d, y ~ x1 + x2, index = c("id", "t"), kernel_lags = 0)
  k2 <- xtkao(d, y ~ x1 + x2, index = c("id", "t"), kernel_lags = 2)
  expect_s3_class(k0, "xtkao"); expect_length(k0$stats, 5)
  expect_equal(k0$stats[4:5], k2$stats[4:5])
  expect_equal(k0$sigma2_nu, k0$omega2_nu)          # bandwidth 0: long-run = contemporaneous
  expect_equal(unname(k0$stats[1]), unname(k0$stats[4]))  # hence modified DF t = unadjusted version
  expect_equal(unname(k0$stats[2]), unname(k0$stats[5]))
  expect_output(print(k0), "Kao")
})

test_that("input checks", {
  d <- sim_panel(); d2 <- d[-1, ]
  expect_error(xtpedroni(d2, y ~ x1, index = c("id", "t")), "balanced")
  expect_error(xtkao(d2, y ~ x1, index = c("id", "t")), "balanced")
  expect_error(xtkao(d, y ~ 1, index = c("id", "t")), "regressor")
  expect_error(xtpedroni(d, y ~ x1, index = "id"), "length 2")
})

test_that("agreement with Stata xtcointtest on the chapter-6 panel (fixed bandwidths)", {
  skip_if(ch6 == "" || !file.exists(ch6), "chapter-6 data not available")
  d <- read.csv(ch6)
  ## set 1: Z ~ LEV LDR RR lnGPR, lags(1) kernel(bartlett 0 / 2)
  p0 <- xtpedroni(d, Z ~ LEV + LDR + RR + lnGPR, index = c("iso", "year"), lags = 1, kernel_lags = 0)
  expect_equal(unname(p0$z[c("group_rho", "group_t", "group_adf")]), c(2.025521, -1.6838803, -0.81624361), tolerance = 1e-6)
  expect_equal(unname(p0$z[c("panel_v", "panel_rho", "panel_t", "panel_adf")]), c(-1.7036008, 0.87688916, -1.7514026, -1.7734208), tolerance = 1e-6)
  p2 <- xtpedroni(d, Z ~ LEV + LDR + RR + lnGPR, index = c("iso", "year"), lags = 1, kernel_lags = 2)
  expect_equal(unname(p2$z[c("group_rho", "group_t")]), c(1.9051592, -1.6520418), tolerance = 1e-6)
  expect_equal(unname(p2$z[c("panel_v", "panel_rho", "panel_t", "panel_adf")]), c(-2.0616892, 0.81700553, -1.7428972, -1.5996663), tolerance = 1e-6)
  ## set 4: AFSI ~ MaPP_cum lnWUI FED lnOIL, lags(1) kernel(bartlett 1)
  p4 <- xtpedroni(d, AFSI ~ MaPP_cum + lnWUI + FED + lnOIL, index = c("iso", "year"), lags = 1, kernel_lags = 1)
  expect_equal(unname(p4$z[c("group_rho", "group_t", "group_adf")]), c(1.5534802, -3.8168484, -2.3863107), tolerance = 1e-6)
  ## Kao, set 1, lags(1) kernel(bartlett 2): four statistics at 1e-6
  k <- xtkao(d, Z ~ LEV + LDR + RR + lnGPR, index = c("iso", "year"), lags = 1, kernel_lags = 2)
  expect_lt(max(abs(unname(k$stats[c(1, 2, 4, 5)]) - c(0.06566431, -0.2500834, 0.17883823, -0.16659824))), 1e-6)
  ## ADF t follows Kao (1999) s_v^2 = RSS/(NT); Stata's ADF t (-1.1560 here)
  ## is a known gap: within 0.02 on this specification (absolute), up to
  ## 0.19 on the other reference specifications (see the man page)
  expect_lt(abs(unname(k$stats[3]) - (-1.1560)), 0.02)
  expect_false(isTRUE(all.equal(unname(k$stats[3]), -1.1560, tolerance = 1e-4)))
})
