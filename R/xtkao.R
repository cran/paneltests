#' Kao (1999) Panel Cointegration Tests
#'
#' Residual-based panel cointegration tests of Kao (1999) with a common
#' cointegrating vector and panel fixed effects: the modified DF t, DF t and
#' ADF t (which allow endogenous regressors through the long-run variance
#' adjustment) and the unadjusted modified DF t and DF t (strict exogeneity,
#' no serial correlation). The statistics are those of Kao (1999); the sample
#' conventions are those of the Stata command \code{xtcointtest kao}, which
#' four of the five statistics reproduce numerically (see Details). The
#' p-values are left-tail normal probabilities, the direction of rejection
#' stated by Kao (1999).
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
#'   (default) for \code{floor(4 (T/100)^(2/9))}. Stata selects a bandwidth
#'   per panel with the Newey-West automatic rule and reports the average;
#'   that rule is not reproduced here.
#'
#' @return An object of class \code{"xtkao"}: a list with \code{stats},
#'   \code{p} (one-sided left-tail p-values \code{pnorm(stat)}), \code{rho},
#'   \code{t_rho}, \code{t_adf}, \code{adf_lags}, \code{kernel_lags},
#'   \code{sigma2_nu}, \code{omega2_nu}, \code{N}, \code{T} and \code{call}.
#'
#' @details
#' Following Stata, the sample first loses one period to the AR(1)
#' regressions of each variable on its own lag (pooled, without constant)
#' whose residuals \eqn{w_{it} = (u_{it}, \epsilon_{it}')'} enter the
#' long-run covariance matrices. The matrices are computed per panel on the
#' within-panel demeaned residuals with a Bartlett kernel, averaged over
#' panels, and the conditional variances
#' \eqn{\hat\sigma_\nu^2 = \hat\sigma_u^2 - \hat\Sigma_{u\epsilon}
#' \hat\Sigma_{\epsilon\epsilon}^{-1} \hat\Sigma_{\epsilon u}} and
#' \eqn{\hat\omega_\nu^2} (same with the long-run matrix) are taken from the
#' averaged matrices. The fixed-effects regression of \eqn{y} on \eqn{x} is
#' estimated on the remaining \eqn{T_0 - 1} periods, and the pooled DF
#' regression \eqn{\hat e_{it} = \rho \hat e_{i,t-1} + \nu_{it}} (no
#' constant) on \eqn{T = T_0 - 2} periods per panel gives \eqn{\hat\rho} and
#' \eqn{t_\rho = (\hat\rho - 1)/se}, with the residual variance
#' \eqn{s_e^2 = (NT)^{-1}\sum\sum \hat\nu_{it}^2} (Kao 1999, p. 7). The
#' ADF regression (Kao 1999, eq. 9) adds \eqn{p} lagged differences and is
#' estimated on the \eqn{T - p} periods per panel for which the lags exist;
#' its t statistic is \eqn{t_{ADF} = (\hat\rho - 1)[\sum_i \hat e_i' Q_i
#' \hat e_i]^{1/2} / s_\nu} with \eqn{s_\nu^2 = (NT)^{-1}\sum_i\sum_t
#' \hat\nu_{itp}^2} (Kao 1999, p. 9-10), where \eqn{T} is the number of
#' periods of the DF step reported as \code{T} (not the number of ADF
#' observations \eqn{N(T - p)}). Versions up to 1.0.7 divided by the
#' number of ADF observations.
#'
#' All five statistics are standard normal under the null and diverge to
#' minus infinity under the alternative (Kao 1999, Theorem 2, Corollary 2
#' and Remarks 6 and 8), so the p-values are \eqn{\Phi(z)}, \eqn{z} being the standardized statistic.
#' Versions up to 1.0.7 returned \eqn{1 - \Phi(|z|)}, which
#' rejects a true null at twice the nominal rate.
#'
#' Kao (1999) tabulates finite-sample sizes for one regressor only (Table
#' 4); the function warns when there are more than seven regressors (no
#' finite-sample evidence at all), when \eqn{N < 2} (the standardization is
#' asymptotic in \eqn{N}) or when the panel has fewer than 10 periods.
#'
#' \strong{Agreement with Stata.} The unadjusted statistics reproduce
#' \code{xtcointtest kao} to at least 1e-6; the modified DF t and DF t do
#' so when the same kernel bandwidth is used for all panels (verified on a
#' six-country panel, four specifications, bandwidths 0, 1 and 2). Two
#' conventions of that command are kept for this agreement and are the
#' package's choice: the innovations entering the long-run covariance
#' matrices are the residuals of pooled AR(1) regressions, demeaned within
#' panel, rather than the first differences of Kao's eqs. (10)-(12), and
#' the fixed-effects regression uses the \eqn{T_0 - 1} periods after that
#' step. The ADF t follows the paper's \eqn{s_\nu^2 = RSS/(NT)}; Stata's
#' ADF t differs from it by up to 0.19 (standardized units) on the
#' reference data, and its internal convention could not be identified from
#' the documentation (known gap). Stata's default bandwidth is chosen per
#' panel by the Newey and West (1994) algorithm; it is not reproduced here.
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
                  kernel_lags = NULL) {
  lag_select <- match.arg(lag_select)
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
  ## 1. pooled AR(1) residuals (no constant) of each variable, t >= 2
  ar_res <- function(v) { out <- rep(NA_real_, nrow(data)); yy <- c(); xx <- c(); rows <- c()
    for (i in seq_len(N)) { r <- rows_i[[i]]; z <- data[[v]][r]; yy <- c(yy, z[-1]); xx <- c(xx, z[-length(z)]); rows <- c(rows, r[-1]) }
    out[rows] <- yy - sum(xx * yy) / sum(xx^2) * xx; out }
  U <- cbind(ar_res(y), sapply(x, ar_res)); keep <- stats::complete.cases(U)
  TT1 <- T0 - 1L
  if (is.null(kernel_lags)) kernel_lags <- floor(4 * (TT1 / 100)^(2 / 9))
  kernel_lags <- as.integer(kernel_lags)
  bart_mat <- function(W, m) { n <- nrow(W); S <- crossprod(W) / n; O <- S
    if (m > 0) for (j in seq_len(m)) { G <- crossprod(W[(j + 1):n, , drop = FALSE], W[1:(n - j), , drop = FALSE]) / n; O <- O + (1 - j / (m + 1)) * (G + t(G)) }
    list(S = S, O = O) }
  S <- 0; O <- 0
  for (i in seq_len(N)) { Wi <- U[keep & data[[pid_col]] == ids[i], , drop = FALSE]; Wi <- scale(Wi, scale = FALSE)
    cv <- bart_mat(Wi, kernel_lags); S <- S + cv$S / N; O <- O + cv$O / N }
  cond <- function(M) as.numeric(M[1, 1] - M[1, -1, drop = FALSE] %*% solve(M[-1, -1, drop = FALSE]) %*% M[-1, 1, drop = FALSE])
  sigma2 <- cond(S); omega2 <- cond(O)
  ## 2. fixed-effects regression on t >= 2
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
  ## 4. pooled ADF regression, common lag
  adf_fit <- function(p, pmax) { ys <- c(); Xs <- NULL
    for (i in seq_len(N)) { z <- E[[i]]; Tz <- length(z); dz <- diff(z); idx <- (pmax + 2):Tz
      Xi <- cbind(z[idx - 1]); for (j in seq_len(p)) Xi <- cbind(Xi, dz[idx - 1 - j])
      ys <- c(ys, z[idx]); Xs <- rbind(Xs, Xi) }
    fit <- stats::lm.fit(Xs, ys); kk <- ncol(Xs); nn <- length(ys); rss <- sum(fit$residuals^2)
    ll <- -nn / 2 * (log(2 * pi) + log(rss / nn) + 1)
    ## s_v^2 = RSS / (N T) with T the number of DF periods (Kao 1999, p. 9-10),
    ## not RSS / (number of ADF observations) as in versions up to 1.0.7
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
  st <- c("Modified DF t" = (sqrt(N) * TT * (rho - 1) + 3 * sqrt(N) * sigma2 / omega2) / sqrt(3 + 36 * sigma2^2 / (5 * omega2^2)),
          "DF t" = (t_rho + add) / den,
          "ADF t" = (t_adf + add) / den,
          "Unadjusted modified DF t" = (sqrt(N) * TT * (rho - 1) + 3 * sqrt(N)) / sqrt(51 / 5),
          "Unadjusted DF t" = sqrt(5 / 4) * t_rho + sqrt(15 * N / 8))
  ## left tail: all five statistics diverge to minus infinity under the
  ## alternative (Kao 1999, Theorem 2, Corollary 2, Remarks 6 and 8)
  out <- list(stats = st, p = stats::pnorm(st), rho = rho, t_rho = t_rho, t_adf = t_adf,
              adf_lags = p, kernel_lags = kernel_lags, sigma2_nu = sigma2, omega2_nu = omega2,
              N = N, T = TT, call = match.call())
  class(out) <- "xtkao"
  out
}

#' @export
#' @method print xtkao
print.xtkao <- function(x, ...) {
  cat("\nKao (1999) panel cointegration tests (common cointegrating vector, fixed effects)\n")
  cat(sprintf("N = %d panels, T = %g periods | Bartlett kernel lags: %d | ADF lags: %d\n", x$N, x$T, x$kernel_lags, x$adf_lags))
  cat("H0: no cointegration\n\n")
  tab <- data.frame(Statistic = sprintf("%9.4f", x$stats), `p-value` = sprintf("%7.4f", x$p), row.names = names(x$stats), check.names = FALSE)
  print(tab)
  cat("\np-values are Phi(stat): all five statistics reject for large negative values (Kao 1999).\n")
  invisible(x)
}
