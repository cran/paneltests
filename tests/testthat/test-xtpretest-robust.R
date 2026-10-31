## Regression tests for 1.0.7: the ROBUST module of xtpretest.
## In 1.0.6 the interacted design matrix contained X together with the
## interactions of ALL panels, which are exactly collinear (the N
## interactions of a regressor sum to the regressor), so solve() failed with
## "system is computationally singular" on every data set, and the error
## aborted the other modules when tests = "ALL".

make_panel2 <- function(n = 6, t = 17, seed = 1) {
  set.seed(seed)
  data.frame(
    id   = rep(seq_len(n), each = t),
    time = rep(seq_len(t), times = n),
    y    = rnorm(n * t),
    x1   = rnorm(n * t),
    x2   = rnorm(n * t)
  )
}

test_that("ROBUST module runs on a clean panel and returns finite statistics", {
  df  <- make_panel2()
  res <- suppressMessages(xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = "ROBUST"))
  expect_false(is.null(res$robust))
  expect_true(is.finite(res$robust$rF2$F))
  expect_true(is.finite(res$robust$rF3$F))
  expect_equal(res$robust$rF2$df1, (6 - 1) * 2)   # (N-1)K restrictions
  expect_equal(res$robust$rF3$df1, 6 - 1)         # N-1 restrictions
  expect_equal(res$robust$rF2$df2, 6 * 17 - (1 + 2 + 5 + 10))
  expect_true(res$robust$rF2$p >= 0 && res$robust$rF2$p <= 1)
  expect_length(res$robust$aliased, 0)
})

test_that("ROBUST HC1 matches a direct sandwich computation", {
  df <- make_panel2(seed = 3)
  res <- suppressMessages(xtpretest(df, y ~ x1, index = c("id", "time"), tests = "ROBUST"))
  pf <- factor(df$id)
  D  <- model.matrix(~ pf - 1)[, -1]
  Xf <- cbind(1, df$x1, D, D * df$x1)
  fit <- lm.fit(Xf, df$y)
  n <- nrow(Xf); k <- ncol(Xf)
  XtXi <- solve(crossprod(Xf))
  meat <- crossprod(Xf * (fit$residuals^2 * n / (n - k)), Xf)
  V <- XtXi %*% meat %*% XtXi
  idx <- (2 + 5 + 1):k
  b <- fit$coefficients[idx]
  F_direct <- as.numeric(t(b) %*% solve(V[idx, idx]) %*% b) / length(idx)
  expect_equal(res$robust$rF2$F, F_direct, tolerance = 1e-10)
})

test_that("time-invariant regressor is dropped as aliased, not fatal", {
  df <- make_panel2(seed = 4)
  df$z <- rep(c(1, 2, 3, 4, 5, 6), each = 17)   # constant within panel
  res <- suppressMessages(xtpretest(df, y ~ x1 + z, index = c("id", "time"), tests = "ROBUST"))
  expect_false(is.null(res$robust))
  expect_true(length(res$robust$aliased) >= 1)
  expect_true(is.finite(res$robust$rF2$F))
})

test_that("a failing module does not abort the others", {
  df <- make_panel2(seed = 5)
  res <- suppressMessages(xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = "ALL"))
  expect_false(is.null(res$hsiao))
  expect_false(is.null(res$csd))
  expect_false(is.null(res$robust))
  expect_null(res$failed)
})
