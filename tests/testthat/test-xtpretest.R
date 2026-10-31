make_panel <- function(n = 5, t = 10, seed = 1) {
  set.seed(seed)
  data.frame(
    id   = rep(seq_len(n), each = t),
    time = rep(seq_len(t), times = n),
    y    = rnorm(n * t),
    x1   = rnorm(n * t),
    x2   = rnorm(n * t)
  )
}

test_that("xtpretest runs Hsiao tests", {
  df  <- make_panel()
  res <- xtpretest(df, y ~ x1, index = c("id", "time"), tests = "HSIAO")
  expect_type(res, "list")
  expect_named(res$hsiao, c("F1", "F1p", "F1_nd", "F1_dd",
                             "F2", "F2p", "F2_nd", "F2_dd",
                             "F3", "F3p", "F3_nd", "F3_dd",
                             "S1", "S2", "S3"), ignore.order = TRUE)
  expect_true(is.numeric(res$hsiao$F1))
  expect_true(res$hsiao$F1p >= 0 && res$hsiao$F1p <= 1)
})

test_that("xtpretest runs CSD test", {
  df  <- make_panel(seed = 5)
  res <- xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = "CSD")
  expect_named(res$csd, c("cd_stat", "pval", "avg_rho"), ignore.order = TRUE)
  expect_true(is.numeric(res$csd$cd_stat))
})

test_that("xtpretest returns a recommendation", {
  df  <- make_panel()
  res <- xtpretest(df, y ~ x1, index = c("id", "time"),
                   tests = c("HSIAO", "CSD"))
  expect_true(is.character(res$recommendation))
  expect_true(nchar(res$recommendation) > 0L)
})

test_that("xtpretest errors on bad index", {
  df <- make_panel()
  expect_error(xtpretest(df, y ~ x1, index = c("bad_id", "time")))
})

test_that("xtpretest summary module works", {
  df  <- make_panel()
  res <- xtpretest(df, y ~ x1, index = c("id", "time"), tests = "SUMMARY")
  expect_type(res$summary, "list")
})

## ---- 1.0.8 checks -----------------------------------------------------

test_that("Swamy statistic centres on the 1/s_ii weighted estimator of Swamy (1970, eq. 6.2)", {
  set.seed(20261030)
  N <- 6; TT <- 15; K <- 2
  df <- data.frame(id = rep(1:N, each = TT), time = rep(1:TT, N))
  df$x1 <- rnorm(N * TT); df$x2 <- rnorm(N * TT)
  df$y <- 1 + 0.5 * df$x1 - 0.3 * df$x2 + rnorm(N * TT, sd = rep(c(0.5, 1, 2, 0.7, 1.5, 3), each = TT))
  res <- suppressMessages(xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = "HETEROGENEITY"))
  ## hand computation with lm(): b_i, s_ii = RSS_i/(T - K - 1), W_i = X_i'X_i / s_ii
  ## (X_i demeaned), beta_hat = (sum W_i)^-1 sum W_i b_i, H = sum (b_i - beta_hat)' W_i (b_i - beta_hat)
  B <- NULL; W <- list()
  for (i in 1:N) {
    di <- df[df$id == i, ]
    fit <- lm(y ~ x1 + x2, data = di)
    s_ii <- sum(residuals(fit)^2) / (TT - K - 1)
    Xd <- scale(cbind(di$x1, di$x2), scale = FALSE)
    B <- rbind(B, coef(fit)[-1]); W[[i]] <- crossprod(Xd) / s_ii
  }
  Wsum <- Reduce(`+`, W)
  beta_hat <- solve(Wsum, Reduce(`+`, lapply(1:N, function(i) W[[i]] %*% B[i, ])))
  H <- sum(sapply(1:N, function(i) { d <- B[i, ] - beta_hat; t(d) %*% W[[i]] %*% d }))
  expect_equal(res$swamy$chi2, H, tolerance = 1e-10)
  expect_equal(res$swamy$beta_hat, as.numeric(beta_hat), tolerance = 1e-10)
  expect_equal(res$swamy$df, K * (N - 1))
  expect_equal(res$swamy$pval, pchisq(H, K * (N - 1), lower.tail = FALSE))
  ## W_i equals the inverse of the slope block of vcov(lm) (Frisch-Waugh)
  fit1 <- lm(y ~ x1 + x2, data = df[df$id == 1, ])
  expect_equal(W[[1]], solve(vcov(fit1)[-1, -1]), tolerance = 1e-10, ignore_attr = TRUE)
  ## and differs from the statistic centred on the unweighted FE estimator
  b_fe <- coef(lm(y ~ factor(id) + x1 + x2, data = df))[c("x1", "x2")]
  H_fe <- sum(sapply(1:N, function(i) { d <- B[i, ] - b_fe; t(d) %*% W[[i]] %*% d }))
  expect_false(isTRUE(all.equal(res$swamy$chi2, H_fe)))
})

test_that("missing values are refused with their count and positions", {
  df <- make_panel(seed = 7)
  df$x1[c(3, 17)] <- NA
  expect_error(xtpretest(df, y ~ x1 + x2, index = c("id", "time")), "2 row\\(s\\).*rows 3, 17")
  df2 <- make_panel(seed = 7); df2$y[8] <- NA
  expect_error(xtpretest(df2, y ~ x1, index = c("id", "time")), "1 row\\(s\\).*rows 8")
})

test_that("exactly collinear regressors are dropped with a warning and K is the rank", {
  df <- make_panel(n = 6, t = 15, seed = 11)
  df$x3 <- 2 * df$x1
  expect_warning(res <- suppressMessages(xtpretest(df, y ~ x1 + x2 + x3, index = c("id", "time"), tests = "HSIAO")), "x3")
  ref <- suppressMessages(xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = "HSIAO"))
  expect_equal(res$dropped, "x3")
  expect_equal(res$hsiao$F2, ref$hsiao$F2); expect_equal(res$hsiao$F2_nd, (6 - 1) * 2); expect_equal(res$hsiao$F2_dd, 6 * (15 - 2 - 1))
  ## hand values from lm(): S1, S2, S3 and F2 = [(S2 - S1)/((N-1)K)] / [S1/(N(T-K-1))]
  S1 <- sum(sapply(1:6, function(i) sum(residuals(lm(y ~ x1 + x2, data = df[df$id == i, ]))^2)))
  S2 <- sum(residuals(lm(y ~ factor(id) + x1 + x2, data = df))^2)
  expect_equal(res$hsiao$F2, ((S2 - S1) / 10) / (S1 / 72), tolerance = 1e-10)
})

test_that("N = 1 and panels shorter than K + 2 are refused", {
  df <- make_panel(seed = 2)
  expect_error(xtpretest(df[df$id == 1, ], y ~ x1, index = c("id", "time")), "at least two panels")
  short <- df[!(df$id == 2 & df$time > 3), ]           # panel 2 keeps 3 obs, K + 2 = 4
  expect_error(xtpretest(short, y ~ x1 + x2, index = c("id", "time")), "fewer than K \\+ 2 = 4.*2")
})

test_that("formula terms log() and I() are evaluated", {
  df <- make_panel(seed = 12); df$x1 <- exp(df$x1)
  df2 <- df; df2$lx <- log(df$x1); df2$sx <- df$x2^2
  r1 <- suppressMessages(xtpretest(df, y ~ log(x1) + I(x2^2), index = c("id", "time"), tests = "HSIAO"))
  r2 <- suppressMessages(xtpretest(df2, y ~ lx + sx, index = c("id", "time"), tests = "HSIAO"))
  expect_equal(r1$hsiao, r2$hsiao)
})

test_that("a regressor constant within one panel stops HSIAO with a message instead of a silent fit", {
  df <- make_panel(n = 5, t = 10, seed = 13)
  df$x2[df$id == 3] <- 1
  res <- suppressMessages(xtpretest(df, y ~ x1 + x2, index = c("id", "time"), tests = c("HSIAO", "HETEROGENEITY")))
  expect_null(res$hsiao); expect_null(res$swamy)
  expect_match(res$failed$HSIAO, "panel 3")
  expect_match(res$failed$HETEROGENEITY, "panel 3")
})

test_that("ROBUST output carries the small-T size note", {
  df <- make_panel(n = 5, t = 10, seed = 14)
  expect_message(xtpretest(df, y ~ x1, index = c("id", "time"), tests = "ROBUST"), "over-rejects at small T")
})
