#' Kao (1999) Panel Cointegration Tests
#'
#' Residual-based panel cointegration tests of Kao (1999) with a common
#' cointegrating vector and panel fixed effects: the modified DF t, DF t and
#' ADF t (which allow endogenous regressors through the long-run variance
#' adjustment) and the unadjusted modified DF t and DF t (strict exogeneity,
#' no serial correlation). The default \code{convention = "kao"} follows
#' Kao (1999) and reproduces the finite-sample sizes of Kao (1999, Table 4);
#' \code{convention = "stata"} follows the Stata command
#' \code{xtcointtest kao} (see Details). The p-values are left-tail normal
#' probabilities, the direction of rejection stated by Kao (1999).
#'
#' @inheritParams xtpedroni
#' @param formula \code{y ~ x1 + x2 + ...}; terms are evaluated with
#'   \code{model.frame} and \code{model.matrix}, so \code{log(x)} or
#'   \code{I(x^2)} are allowed.
#' @param lags Integer. Number of lagged differences in the pooled ADF
#'   regression (default 1); with \code{lag_select} it is the maximum and
#'   the lag is chosen by AIC or BIC over 1, ..., \code{lags} on the common
#'   sample of the largest regression.
#' @param kernel_lags Integer bandwidth of the Bartlett kernel for the
#'   long-run covariance matrices, applied to every panel, or \code{NULL}
#'   (default) for \code{floor(4 (T/100)^(2/9))} with \eqn{T = T_0 - 1}.
#'   Stata selects a bandwidth per panel with the Newey-West automatic rule
#'   and reports the average; that rule is not reproduced here.
#' @param convention \code{"kao"} (default) or \code{"stata"}: the sample
#'   and innovation conventions described in Details.
#'
#' @return An object of class \code{"xtkao"}: a list with \code{stats},
#'   \code{p} (one-sided left-tail p-values \code{pnorm(stat)}), \code{rho},
#'   \code{t_rho}, \code{t_adf}, \code{adf_lags}, \code{kernel_lags},
#'   \code{sigma2_nu}, \code{omega2_nu}, \code{N}, \code{T} (periods of
#'   the DF regression), \code{T0} (periods per panel in the data),
#'   \code{convention} and \code{call}.
#'
#' @details
#' With \eqn{T_0} periods per panel, the two conventions are:
#'
#' \strong{\code{convention = "kao"}.} The fixed-effects regression of
#' \eqn{y} on \eqn{x} (Kao 1999, eq. 4) uses all \eqn{T_0} periods. The
#' innovations entering the long-run covariance matrices are the first
#' differences \eqn{\hat w_{it} = (\Delta y_{it}, \Delta x_{it}')'}, not
#' demeaned (Kao 1999, eqs. 10-12). The pooled DF regression
#' \eqn{\hat e_{it} = \rho \hat e_{i,t-1} + \nu_{it}} (eqs. 7-8) and the
#' pooled ADF regression with \eqn{p} common lagged differences (eq. 9) use
#' the residual variances \eqn{s_e^2} and \eqn{s_\nu^2} computed as the
#' residual sum of squares divided by \eqn{N (T_0 - 1)}, the number of
#' observations of the DF regression, and the modified DF statistics use
#' \eqn{\sqrt{N} T_0 (\hat\rho - 1)}. With these choices the five
#' statistics reproduce the sizes of Kao (1999, Table 4) in a replication of
#' its design (two ADF lags, Bartlett window 5, \eqn{T = 10, 25, 50},
#' \eqn{N = 2} to 50, 5000 replications per cell): mean absolute
#' differences of 0.002 to 0.009 across the five columns
#' (\code{inst/validation/kao_table4.R}).
#'
#' \strong{Size in short panels.} Both conventions use the asymptotic
#' N(0, 1) critical values and over-reject in short panels, as Kao (1999,
#' Table 4) reports. In a simulation with \eqn{N = 6}, \eqn{T_0 = 17}, one
#' regressor, no cointegration, lags chosen by AIC up to 2 and the default
#' bandwidth (1000 replications), the rejection rates at the nominal 5
#' percent level were 0.13 to 0.19 with \code{"kao"} and 0.05 to 0.17 with
#' \code{"stata"} across the five statistics, and 0.35 to 0.50 and 0.21 to
#' 0.43 with four regressors. Rejections in such samples are weak evidence
#' of cointegration.
#'
#' \strong{\code{convention = "stata"}.} As in \code{xtcointtest kao}, the
#' sample first loses one period to the pooled AR(1) regressions (without
#' constant) of each variable on its own lag; their residuals, demeaned
#' within panel, are the innovations. The fixed-effects regression uses the
#' remaining \eqn{T_0 - 1} periods, the DF regression \eqn{T_0 - 2}, and
#' the residual variances and the modified DF statistics use
#' \eqn{T = T_0 - 2}. This is the computation of versions 1.0.8 and earlier.
#' The unadjusted statistics reproduce \code{xtcointtest kao} to at least
#' 1e-6; the modified DF t and DF t do so when the same kernel bandwidth is
#' used for all panels (verified on a six-country panel, four
#' specifications, bandwidths 0, 1 and 2). Stata's ADF t differs by up to
#' 0.19 (standardized units) on the reference data; its internal convention
#' could not be identified from the documentation. The AR(1) regressions of
#' this convention have no constant, so its statistics change when a
#' panel-specific constant is added to a variable; the statistics of
#' \code{convention = "kao"} do not.
#'
#' \strong{Input checks.} The function stops when a (panel, time) pair is
#' duplicated, when a numeric or \code{Date} time index is not equally
#' spaced within a panel, when the response has no within-panel variation,
#' when a regressor is constant within panels or collinear with the others
#' after removing the fixed effects, when \code{kernel_lags} is not a
#' non-negative integer smaller than the number of innovations per panel,
#' and when \code{lags} leaves too few observations for the pooled ADF
#' regression or makes it fit exactly.
#'
#' In both conventions the long-run covariance matrices are computed per
#' panel with a Bartlett kernel and averaged over panels, and the
#' conditional variances
#' \eqn{\hat\sigma_\nu^2 = \hat\sigma_u^2 - \hat\Sigma_{u\epsilon}
#' \hat\Sigma_{\epsilon\epsilon}^{-1} \hat\Sigma_{\epsilon u}} and
#' \eqn{\hat\omega_\nu^2} (same with the long-run matrix) are taken from
#' the averaged matrices. The ADF t statistic is
#' \eqn{t_{ADF} = (\hat\rho - 1)[\sum_i \hat e_i' Q_i \hat e_i]^{1/2} /
#' s_\nu} (Kao 1999, p. 9).
#'
#' All five statistics are standard normal under the null and diverge to
#' minus infinity under the alternative (Kao 1999, Theorem 2, Corollary 2
#' and Remarks 6 and 8), so the p-values are \eqn{\Phi(z)}, \eqn{z} being
#' the standardized statistic. Kao (1999) tabulates finite-sample sizes for
#' one regressor only (Table 4); the function warns when there are more than
#' seven regressors (no finite-sample evidence at all), when \eqn{N < 2}
#' (the standardization is asymptotic in \eqn{N}) or when the panel has
#' fewer than 10 periods.
#'
#' @references
#' Kao, C. (1999). Spurious regression and residual-based tests for
#' cointegration in panel data. \emph{Journal of Econometrics}, 90(1), 1-44.
#' \doi{10.1016/S0304-4076(98)00023-2}
#'
#' @examples
#' set.seed(2)
#' N <- 6; TT <- 20
#' d <- do.call(rbind, lapply(1:N, function(i) {
#'   x <- cumsum(rnorm(TT)); y <- i + 0.5 * x + rnorm(TT)
#'   data.frame(id = i, t = 1:TT, y = y, x = x) }))
#' xtkao(d, y ~ x, index = c("id", "t"), kernel_lags = 1)
#' @export
xtkao <- function(data, formula, index, lags = 1L, lag_select = c("none", "aic", "bic"),
                  kernel_lags = NULL, convention = c("kao", "stata")) {
  lag_select <- match.arg(lag_select)
  convention <- match.arg(convention)
  if (!inherits(data, "data.frame")) stop("'data' must be a data.frame.")
  if (length(index) != 2L) stop("'index' must be a character vector of length 2.")
  pid_col <- index[1L]; tid_col <- index[2L]
  if (!all(c(pid_col, tid_col) %in% names(data))) stop("index columns not found in 'data'.")
  data <- .coint_frame(formula, data, pid_col, tid_col, "xtkao")
  y <- "y"; x <- setdiff(names(data), c(pid_col, tid_col, "y"))
  ids <- unique(data[[pid_col]]); N <- length(ids)
  Ts <- table(data[[pid_col]])
  if (length(unique(Ts)) != 1L) stop("xtkao requires a balanced panel")
  T0 <- as.integer(Ts[1])
  if (length(x) > 7L) warning(sprintf("%d regressors: Kao (1999) reports finite-sample evidence for one regressor only and no evidence exists for more than seven; the asymptotic standardization is unreliable", length(x)))
  if (N < 2L) warning("N = 1 panel: the standardization is asymptotic in N and the p-values are not meaningful for a single panel")
  if (T0 < 10L) warning(sprintf("only %d periods per panel: the standardized statistics are unreliable for short panels", T0))
  rows_i <- split(seq_len(nrow(data)), data[[pid_col]])[as.character(ids)]
  ## within-panel variation of the response and rank of the regressors
  dm0 <- function(v) stats::ave(v, data[[pid_col]], FUN = function(z) z - mean(z))
  if (sum(dm0(data[[y]])^2) <= 1e-12 * max(sum(data[[y]]^2), .Machine$double.xmin))
    stop("xtkao: the response has no variation within panels")
  qx <- qr(vapply(x, function(v) dm0(data[[v]]), numeric(nrow(data))))
  if (qx$rank < length(x))
    stop(sprintf("xtkao: regressor(s) %s are constant within panels or collinear with the others after removing the fixed effects",
                 paste(x[qx$pivot[(qx$rank + 1):length(x)]], collapse = ", ")))
  TT1 <- T0 - 1L
  if (is.null(kernel_lags)) kernel_lags <- floor(4 * (TT1 / 100)^(2 / 9))
  kernel_lags <- .check_kernel_lags(kernel_lags, TT1, "xtkao")
  lags <- .check_lags(lags, "xtkao")
  bart_mat <- function(W, m) { n <- nrow(W); S <- crossprod(W) / n; O <- S
    if (m > 0) for (j in seq_len(m)) { G <- crossprod(W[(j + 1):n, , drop = FALSE], W[1:(n - j), , drop = FALSE]) / n; O <- O + (1 - j / (m + 1)) * (G + t(G)) }
    list(S = S, O = O) }
  if (convention == "kao") {
    ## 1. innovations w_it = (dy_it, dx_it')' (Kao 1999, eqs. 10-12), not demeaned
    W_i <- lapply(rows_i, function(r) apply(as.matrix(data[r, c(y, x), drop = FALSE]), 2, diff))
    keep <- rep(TRUE, nrow(data))
  } else {
    ## 1. pooled AR(1) residuals (no constant) of each variable, t >= 2,
    ##    demeaned within panel (xtcointtest kao)
    ar_res <- function(v) { out <- rep(NA_real_, nrow(data)); yy <- c(); xx <- c(); rows <- c()
      for (i in seq_len(N)) { r <- rows_i[[i]]; z <- data[[v]][r]; yy <- c(yy, z[-1]); xx <- c(xx, z[-length(z)]); rows <- c(rows, r[-1]) }
      out[rows] <- yy - sum(xx * yy) / sum(xx^2) * xx; out }
    U <- cbind(ar_res(y), sapply(x, ar_res)); keep <- stats::complete.cases(U)
    W_i <- lapply(seq_len(N), function(i) scale(U[keep & data[[pid_col]] == ids[i], , drop = FALSE], scale = FALSE))
  }
  S <- 0; O <- 0
  for (i in seq_len(N)) { cv <- bart_mat(W_i[[i]], kernel_lags); S <- S + cv$S / N; O <- O + cv$O / N }
  cond <- function(M) as.numeric(M[1, 1] - M[1, -1, drop = FALSE] %*% solve(M[-1, -1, drop = FALSE]) %*% M[-1, 1, drop = FALSE])
  sigma2 <- cond(S); omega2 <- cond(O)
  ## 2. fixed-effects regression: t = 1, ..., T0 (kao) or t >= 2 (stata)
  dd <- data[keep, ]
  dm <- function(v) stats::ave(v, dd[[pid_col]], FUN = function(z) z - mean(z))
  Yd <- dm(dd[[y]]); Xd <- sapply(x, function(v) dm(dd[[v]])); Xd <- matrix(Xd, ncol = length(x))
  b <- solve(crossprod(Xd), crossprod(Xd, Yd)); e <- as.numeric(Yd - Xd %*% b)
  E <- split(e, dd[[pid_col]])[as.character(ids)]
  ## 3. pooled DF regression
  el <- unlist(lapply(E, function(z) z[-length(z)])); et <- unlist(lapply(E, function(z) z[-1]))
  rho <- sum(el * et) / sum(el^2); res <- et - rho * el; n <- length(et)
  t_rho <- (rho - 1) / sqrt(sum(res^2) / n / sum(el^2))
  TT <- n / N
  ## T multiplying (rho - 1) in the modified DF statistics: the series length
  ## T0 (kao; Kao 1999, pp. 7-8 and Table 4) or the DF periods (stata)
  Tm <- if (convention == "kao") T0 else TT
  ## 4. pooled ADF regression, common lag
  adf_fit <- function(p, pmax) { ys <- c(); Xs <- NULL
    if (length(E[[1L]]) < pmax + 2L)
      stop(sprintf("xtkao: lags = %d leaves no observation for the pooled ADF regression (%d periods after the fixed-effects step); use fewer lags",
                   pmax, length(E[[1L]])))
    for (i in seq_len(N)) { z <- E[[i]]; Tz <- length(z); dz <- diff(z); idx <- (pmax + 2):Tz
      Xi <- cbind(z[idx - 1]); for (j in seq_len(p)) Xi <- cbind(Xi, dz[idx - 1 - j])
      ys <- c(ys, z[idx]); Xs <- rbind(Xs, Xi) }
    kk <- ncol(Xs); nn <- length(ys)
    if (nn <= kk || qr(Xs)$rank < kk)
      stop(sprintf("xtkao: lags = %d leaves %d observation(s) for the %d coefficients of the pooled ADF regression; use fewer lags",
                   pmax, nn, kk))
    fit <- stats::lm.fit(Xs, ys); rss <- sum(fit$residuals^2)
    if (rss <= 1e-12 * max(sum(ys^2), .Machine$double.xmin))
      stop(sprintf("xtkao: the pooled ADF regression with %d lag(s) fits exactly (zero residual variance); use fewer lags or a longer panel", p))
    ll <- -nn / 2 * (log(2 * pi) + log(rss / nn) + 1)
    ## s_v^2 = RSS / (N TT), TT the number of DF periods: T0 - 1 (kao) or
    ## T0 - 2 (stata); with TT = T0 - 1 this reproduces Kao (1999, Table 4)
    s2v <- rss / (N * TT)
    list(p = p, aic = -2 * ll + 2 * kk, bic = -2 * ll + log(nn) * kk,
         t = (fit$coefficients[1] - 1) / sqrt(s2v * solve(crossprod(Xs))[1, 1])) }
  if (lag_select == "none") { p <- as.integer(lags); af <- adf_fit(p, p) } else {
    ## lag chosen on the common sample of the largest regression, then the
    ## chosen regression is re-estimated on its own (larger) sample, as Stata does
    cands <- lapply(seq_len(lags), function(p) adf_fit(p, lags)); crit <- sapply(cands, `[[`, lag_select)
    p <- cands[[which.min(crit)]]$p; af <- adf_fit(p, p) }
  t_adf <- as.numeric(af$t)
  ## 5. statistics
  den <- sqrt(omega2 / (2 * sigma2) + 3 * sigma2 / (10 * omega2)); add <- sqrt(6 * N * sigma2) / (2 * sqrt(omega2))
  st <- c("Modified DF t" = (sqrt(N) * Tm * (rho - 1) + 3 * sqrt(N) * sigma2 / omega2) / sqrt(3 + 36 * sigma2^2 / (5 * omega2^2)),
          "DF t" = (t_rho + add) / den,
          "ADF t" = (t_adf + add) / den,
          "Unadjusted modified DF t" = (sqrt(N) * Tm * (rho - 1) + 3 * sqrt(N)) / sqrt(51 / 5),
          "Unadjusted DF t" = sqrt(5 / 4) * t_rho + sqrt(15 * N / 8))
  ## left tail: all five statistics diverge to minus infinity under the
  ## alternative (Kao 1999, Theorem 2, Corollary 2, Remarks 6 and 8)
  out <- list(stats = st, p = stats::pnorm(st), rho = rho, t_rho = t_rho, t_adf = t_adf,
              adf_lags = p, kernel_lags = kernel_lags, sigma2_nu = sigma2, omega2_nu = omega2,
              N = N, T = TT, T0 = T0, convention = convention, call = match.call())
  class(out) <- "xtkao"
  out
}

#' @export
#' @method print xtkao
print.xtkao <- function(x, ...) {
  cat("\nKao (1999) panel cointegration tests (common cointegrating vector, fixed effects)\n")
  cat(sprintf("N = %d panels, %d periods per panel (%g in the DF regression) | Bartlett kernel lags: %d | ADF lags: %d\n", x$N, x$T0, x$T, x$kernel_lags, x$adf_lags))
  cat(sprintf("Convention: %s\n", if (identical(x$convention, "stata")) "stata (xtcointtest kao)" else "kao (Kao 1999)"))
  cat("H0: no cointegration\n\n")
  tab <- data.frame(Statistic = sprintf("%9.4f", x$stats), `p-value` = sprintf("%7.4f", x$p), row.names = names(x$stats), check.names = FALSE)
  print(tab)
  cat("\np-values are Phi(stat): all five statistics reject for large negative values (Kao 1999).\n")
  invisible(x)
}
