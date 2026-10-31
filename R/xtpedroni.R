#' Pedroni (1999, 2004) Panel Cointegration Tests
#'
#' Residual-based panel cointegration tests of Pedroni (1999, 2004) with
#' panel-specific cointegrating vectors: four "panel" (same-AR) statistics
#' (modified variance ratio, modified Phillips-Perron t, Phillips-Perron t,
#' ADF t) and three "group" (panel-specific-AR) statistics (modified
#' Phillips-Perron t, Phillips-Perron t, ADF t). The statistics are those of
#' Pedroni (1999, Table 1); the sample conventions are those of the Stata
#' command \code{xtcointtest pedroni}, which the statistics reproduce
#' numerically (see Details). The p-values follow the tails stated by
#' Pedroni (1999, p. 668).
#'
#' @param data A balanced panel \code{data.frame} in long format.
#' @param formula \code{y ~ x1 + x2 + ...}; the cointegrating regression is
#'   estimated panel by panel with a constant (and a trend if
#'   \code{trend = TRUE}). Terms are evaluated with \code{model.frame} and
#'   \code{model.matrix}, so transformations such as \code{log(x)} or
#'   \code{I(x^2)} are allowed; each column of the model matrix counts as
#'   one regressor.
#' @param index Character vector \code{c("panel_id", "time_id")}.
#' @param trend Logical. Include panel-specific linear trends (default
#'   \code{FALSE}).
#' @param lags Integer. Number of lagged differences in the ADF regressions
#'   (default 1). With \code{lag_select = "aic"} or \code{"bic"} it is the
#'   maximum, and the common lag (shared by all panels) is chosen on the
#'   pooled ADF regression over 1, ..., \code{lags}, as Stata does.
#' @param lag_select \code{"none"} (default), \code{"aic"} or \code{"bic"}.
#' @param kernel_lags Integer bandwidth of the Bartlett kernel used for the
#'   long-run variances, or \code{NULL} (default) for the rule
#'   \code{floor(4 (T/100)^(2/9))} where T is the number of periods after
#'   differencing. Stata reports the bandwidth it used ("Lags:" in the
#'   header); pass that value to reproduce Stata output exactly.
#'
#' @return An object of class \code{"xtpedroni"}: a list with \code{stats}
#'   (standardized statistics), \code{p} (one-sided p-values in the tail
#'   stated by Pedroni (1999, p. 668): \code{1 - pnorm(z)} for the panel
#'   variance ratio, which diverges to plus infinity under the alternative,
#'   and \code{pnorm(z)} for the other six statistics, which diverge to
#'   minus infinity), \code{tail} (the tail used for each statistic),
#'   \code{raw} (the unstandardized statistics), \code{moments} (means and
#'   variances used), \code{N}, \code{T}, \code{k}, \code{adf_lags},
#'   \code{kernel_lags} and \code{call}.
#'
#' @details
#' Let \eqn{\hat e_{it}} be the residuals of the panel-specific level
#' regression on all \eqn{T_0} periods. All sums run over the \eqn{T = T_0 -
#' 1} periods for which first differences exist. For each panel the DF
#' regression \eqn{\Delta\hat e_{it} = (\rho_i - 1) \hat e_{i,t-1} + u_{it}}
#' (no constant) is estimated; its residuals are demeaned and give
#' \eqn{\hat s_i^2 = \sum u_{it}^2 / T}, the Bartlett long-run variance
#' \eqn{\hat\sigma_i^2} of \eqn{u_{it}} (autocovariances divided by
#' \eqn{T}) and \eqn{\hat\lambda_i = (\hat\sigma_i^2 - \hat s_i^2)/2}.
#' \eqn{\hat L_{11i}^2} is the conditional long-run variance of
#' \eqn{\Delta y} given \eqn{\Delta x}, computed from the Bartlett long-run
#' covariance matrix of the demeaned first differences. The ADF regression
#' (panel-specific \eqn{\rho_i}, no constant) on the \eqn{T - p} periods
#' with \code{lags} \eqn{= p} lagged differences available gives, from its
#' demeaned residuals, \eqn{\hat s_i^{*2}}; the ADF statistics use the sums
#' \eqn{\sum \hat e_{i,t-1}\Delta\hat e_{it}} and \eqn{\sum \hat e_{i,t-1}^2}
#' over that reduced sample.
#'
#' Standardization uses \eqn{(\chi - \mu\sqrt N)/\sqrt\nu}. For \eqn{k = 1}
#' regressor the moments are those of Pedroni (2004, Corollary 1): case 2
#' (estimated fixed effects) when \code{trend = FALSE} and case 3 (estimated
#' fixed effects and trends) when \code{trend = TRUE}. For \eqn{k = 2,
#' \dots, 7} they are those of Pedroni (1999, Table 2) with \eqn{m = k},
#' rows "heterogeneous intercepts included" and "heterogeneous deterministic
#' trends and intercepts included". More than seven regressors are refused
#' because no moments are tabulated. The moments are large-\eqn{T}
#' approximations; the function warns when \eqn{N < 2} (the standardization
#' is asymptotic in \eqn{N}) or when the panel has fewer than 10 periods.
#'
#' \strong{Direction of rejection and p-values.} The raw statistics are
#' Pedroni's: under the alternative the panel variance ratio diverges to
#' plus infinity and the other six statistics diverge to minus infinity
#' (Pedroni 1999, p. 668). The p-values are therefore
#' \eqn{1 - \Phi(z)} for the panel variance ratio and \eqn{\Phi(z)} for
#' the other six. Versions up to 1.0.7 returned \eqn{1 - \Phi(|z|)} for all
#' seven statistics, which rejects a true null at twice the nominal rate
#' and gives a small p-value to a statistic of the wrong sign.
#'
#' \strong{Agreement with Stata.} With the same kernel bandwidth and ADF
#' lag, all seven standardized statistics reproduce \code{xtcointtest
#' pedroni} to at least 1e-6 (verified on a six-country panel, four
#' specifications with four to six regressors, bandwidths 0, 1 and 2).
#' Two conventions of that command are kept because they are needed for
#' this agreement and are stated here as the package's choice: the
#' residuals of the DF and ADF regressions are demeaned before the
#' variances are computed, and the panel ADF t uses
#' \eqn{\tilde s^{*2} = N^{-1}\sum_i \hat L_{11i}^{-2}\hat s_i^{*2}}
#' (weighted), whereas Pedroni (1999, Table 1) writes the unweighted
#' average. The p-values are not part of the comparison; they follow the
#' paper as described above. Stata's default bandwidth is chosen
#' per panel by the Newey and West (1994) algorithm and the header reports
#' the average; that algorithm is not reproduced here, so the default rule
#' of this function (or an integer \code{kernel_lags}) will in general differ
#' from a Stata run with \code{kernel(bartlett nwest)} unless all panels
#' received the same bandwidth (as when the reported average is an integer
#' such as 0.00).
#'
#' @references
#' Pedroni, P. (1999). Critical values for cointegration tests in
#' heterogeneous panels with multiple regressors. \emph{Oxford Bulletin of
#' Economics and Statistics}, 61(S1), 653-670.
#' \doi{10.1111/1468-0084.61.s1.14}
#'
#' Pedroni, P. (2004). Panel cointegration: asymptotic and finite sample
#' properties of pooled time series tests with an application to the PPP
#' hypothesis. \emph{Econometric Theory}, 20(3), 597-625.
#' \doi{10.1017/S0266466604203073}
#'
#' @examples
#' set.seed(1)
#' N <- 6; TT <- 20
#' d <- do.call(rbind, lapply(1:N, function(i) {
#'   x <- cumsum(rnorm(TT)); y <- 1 + 0.5 * x + rnorm(TT)
#'   data.frame(id = i, t = 1:TT, y = y, x = x) }))
#' r <- xtpedroni(d, y ~ x, index = c("id", "t"), kernel_lags = 2)
#' r$stats
#' @export
xtpedroni <- function(data, formula, index, trend = FALSE, lags = 1L,
                      lag_select = c("none", "aic", "bic"), kernel_lags = NULL) {
  lag_select <- match.arg(lag_select)
  if (!inherits(data, "data.frame")) stop("'data' must be a data.frame.")
  if (length(index) != 2L) stop("'index' must be a character vector of length 2.")
  pid_col <- index[1L]; tid_col <- index[2L]
  if (!all(c(pid_col, tid_col) %in% names(data))) stop("index columns not found in 'data'.")
  data <- .coint_frame(formula, data, pid_col, tid_col, "xtpedroni")
  y <- "y"; x <- setdiff(names(data), c(pid_col, tid_col, "y")); k <- length(x)
  if (k > 7L) stop("moments are available for at most 7 regressors (Pedroni 1999, Table 2, m = 2, ..., 7; Pedroni 2004 for m = 1)")
  ids <- unique(data[[pid_col]]); N <- length(ids)
  Ts <- table(data[[pid_col]])
  if (length(unique(Ts)) != 1L) stop("xtpedroni requires a balanced panel")
  T0 <- as.integer(Ts[1]); TT <- T0 - 1L
  if (N < 2L) warning("N = 1 panel: the standardization (chi - mu sqrt(N)) / sqrt(nu) is asymptotic in N and the p-values are not meaningful for a single panel")
  if (T0 < 10L) warning(sprintf("only %d periods per panel: the moments are large-T approximations and the standardized statistics are unreliable for short panels", T0))
  if (is.null(kernel_lags)) kernel_lags <- floor(4 * (TT / 100)^(2 / 9))
  kernel_lags <- as.integer(kernel_lags)
  bart <- function(u, m) { n <- length(u); s <- sum(u^2) / n
    if (m > 0) for (j in seq_len(m)) s <- s + 2 * (1 - j / (m + 1)) * sum(u[(j + 1):n] * u[1:(n - j)]) / n
    s }
  bart_mat <- function(W, m) { n <- nrow(W); S <- crossprod(W) / n
    if (m > 0) for (j in seq_len(m)) { G <- crossprod(W[(j + 1):n, , drop = FALSE], W[1:(n - j), , drop = FALSE]) / n; S <- S + (1 - j / (m + 1)) * (G + t(G)) }
    S }
  E <- vector("list", N); L11 <- s2 <- sig2 <- lam <- numeric(N)
  for (i in seq_len(N)) {
    dd <- data[data[[pid_col]] == ids[i], ]
    X <- as.matrix(dd[, x, drop = FALSE])
    Zdet <- if (trend) cbind(1, seq_len(T0)) else matrix(1, T0, 1)
    e <- stats::lm.fit(cbind(Zdet, X), dd[[y]])$residuals
    E[[i]] <- e
    el <- e[-T0]; de <- diff(e)
    ## DF regression without constant; residuals are demeaned before the
    ## variance and long-run variance are computed (convention of Stata)
    u <- de - el * sum(el * de) / sum(el^2); u <- u - mean(u)
    s2[i] <- sum(u^2) / length(u); sig2[i] <- bart(u, kernel_lags); lam[i] <- (sig2[i] - s2[i]) / 2
    W <- cbind(diff(dd[[y]]), apply(X, 2, diff)); W <- scale(W, scale = FALSE)
    Om <- bart_mat(W, kernel_lags)
    L11[i] <- as.numeric(Om[1, 1] - Om[1, -1, drop = FALSE] %*% solve(Om[-1, -1, drop = FALSE]) %*% Om[-1, 1, drop = FALSE])
  }
  ## ADF lag: common, chosen on the pooled ADF regression (Stata), lags 1..lags
  adf_pool <- function(p, pmax) { ys <- c(); Xs <- NULL
    for (i in seq_len(N)) { e <- E[[i]]; de <- diff(e); idx <- (pmax + 2):T0
      Xi <- cbind(e[idx - 1]); for (j in seq_len(p)) Xi <- cbind(Xi, de[idx - 1 - j])
      ys <- c(ys, e[idx]); Xs <- rbind(Xs, Xi) }
    fit <- stats::lm.fit(Xs, ys); n <- length(ys); kk <- ncol(Xs); rss <- sum(fit$residuals^2)
    ll <- -n / 2 * (log(2 * pi) + log(rss / n) + 1)
    c(aic = -2 * ll + 2 * kk, bic = -2 * ll + log(n) * kk) }
  if (lag_select == "none") p <- as.integer(lags) else {
    crit <- sapply(seq_len(lags), function(p) adf_pool(p, lags)[lag_select]); p <- which.min(crit) }
  ## panel-specific ADF regressions (with constant) on the reduced sample
  anum <- aden <- sstar2 <- numeric(N)
  for (i in seq_len(N)) { e <- E[[i]]; de <- diff(e); idx <- (p + 2):T0
    el2 <- e[idx - 1]; de2 <- de[idx - 1]
    Za <- cbind(el2); for (j in seq_len(p)) Za <- cbind(Za, de[idx - 1 - j])
    ua <- stats::lm.fit(Za, de2)$residuals; ua <- ua - mean(ua)
    sstar2[i] <- sum(ua^2) / length(ua); anum[i] <- sum(el2 * de2); aden[i] <- sum(el2^2) }
  num <- sapply(E, function(e) sum(e[-T0] * diff(e))) - TT * lam
  den <- sapply(E, function(e) sum(e[-T0]^2))
  w <- 1 / L11
  raw <- c(panel_v = TT^2 * N^1.5 / sum(w * den),
           panel_rho = TT * sqrt(N) * sum(w * num) / sum(w * den),
           panel_t = sum(w * num) / sqrt(mean(sig2 * w) * sum(w * den)),
           panel_adf = sum(w * anum) / sqrt(mean(sstar2 * w) * sum(w * aden)),
           group_rho = TT / sqrt(N) * sum(num / den),
           group_t = sum(num / sqrt(sig2 * den)) / sqrt(N),
           group_adf = sum(anum / sqrt(sstar2 * aden)) / sqrt(N))
  mom <- .pedroni_moments(k, if (trend) "trend" else "intercept")
  z <- (raw - mom$mu * sqrt(N)) / sqrt(mom$v)
  stats <- c("Modified variance ratio (panel)" = z[["panel_v"]], "Modified PP t (panel)" = z[["panel_rho"]],
             "PP t (panel)" = z[["panel_t"]], "ADF t (panel)" = z[["panel_adf"]],
             "Modified PP t (group)" = z[["group_rho"]], "PP t (group)" = z[["group_t"]], "ADF t (group)" = z[["group_adf"]])
  ## Tails (Pedroni 1999, p. 668): the panel variance statistic diverges to
  ## plus infinity under the alternative (right tail); the other six diverge
  ## to minus infinity (left tail).
  tail <- c("right", rep("left", 6L)); names(tail) <- names(stats)
  pval <- ifelse(tail == "right", stats::pnorm(stats, lower.tail = FALSE), stats::pnorm(stats))
  names(pval) <- names(stats)
  out <- list(stats = stats, p = pval, tail = tail, raw = raw, z = z, moments = mom,
              N = N, T = TT, k = k, adf_lags = p, kernel_lags = kernel_lags, trend = trend,
              L11 = L11, lambda = lam, sigma2 = sig2, s2 = s2, sstar2 = sstar2, call = match.call())
  class(out) <- "xtpedroni"
  out
}

#' @keywords internal
.pedroni_moments <- function(k, case = c("intercept", "trend", "none")) {
  case <- match.arg(case)
  stat_names <- c("panel_v", "panel_rho", "panel_t", "panel_adf", "group_rho", "group_t", "group_adf")
  if (k == 1L) {
    ## Pedroni (2004), Corollary 1, p. 606: case 1 (standard case, no
    ## deterministics), case 2 (estimated fixed effects), case 3 (estimated
    ## fixed effects and trends). The five printed moments are for
    ## (panel v, panel rho, panel t, group rho, group t); the ADF t
    ## statistics share the moments of the corresponding t statistics, as
    ## in Pedroni (1999, Table 2). Versions up to 1.0.7 held the Pedroni
    ## (1999) m = 2 rows for the "none" and "trend" cases.
    mu <- switch(case,
      intercept = c(8.62, -6.02, -1.73, -1.73, -9.05, -2.03, -2.03),
      none      = c(4.00, -2.77, -1.01, -1.01, -6.84, -1.39, -1.39),
      trend     = c(17.86, -10.54, -2.29, -2.29, -13.65, -2.53, -2.53))
    v <- switch(case,
      intercept = c(60.75, 31.27, 0.93, 0.93, 35.98, 0.66, 0.66),
      none      = c(27.81, 24.91, 1.50, 1.50, 26.78, 0.78, 0.78),
      trend     = c(101.68, 39.52, 0.66, 0.66, 50.91, 0.56, 0.56))
  } else {
    ## Pedroni (1999), Table 2, m = k regressors (m = 2, ..., 7), columns:
    ## standard case, heterogeneous intercepts, heterogeneous trends and
    ## intercepts. Checked cell by cell against the printed table (the
    ## m = 3 intercept panel t mean is -2.576; versions up to 1.0.7 held -2.567).
    j <- k - 1L
    MU <- list(
      cbind(c(6.982,-6.388,-1.662,-1.662,-9.889,-1.992,-1.992), c(11.754,-9.495,-2.177,-2.177,-12.938,-2.453,-2.453), c(21.162,-14.011,-2.648,-2.648,-17.359,-2.872,-2.872)),
      cbind(c(10.402,-10.191,-2.156,-2.156,-13.865,-2.44,-2.44), c(15.197,-13.256,-2.576,-2.576,-16.888,-2.827,-2.827), c(24.556,-17.6,-2.967,-2.967,-21.116,-3.179,-3.179)),
      cbind(c(14.254,-14.136,-2.571,-2.571,-17.834,-2.819,-2.819), c(18.91,-17.163,-2.93,-2.93,-20.841,-3.157,-3.157), c(28.046,-21.287,-3.262,-3.262,-24.93,-3.464,-3.464)),
      cbind(c(18.198,-18.042,-2.926,-2.926,-21.805,-3.151,-3.151), c(22.715,-21.013,-3.241,-3.241,-24.775,-3.452,-3.452), c(31.738,-25.13,-3.545,-3.545,-28.849,-3.737,-3.737)),
      cbind(c(22.169,-21.985,-3.244,-3.244,-25.75,-3.45,-3.45), c(26.603,-24.944,-3.531,-3.531,-28.72,-3.726,-3.726), c(35.537,-28.981,-3.806,-3.806,-32.716,-3.986,-3.986)),
      cbind(c(26.12,-25.889,-3.533,-3.533,-29.627,-3.723,-3.723), c(30.457,-28.795,-3.795,-3.795,-32.538,-3.976,-3.976), c(39.231,-32.756,-4.047,-4.047,-36.494,-4.217,-4.217)))
    VV <- list(
      cbind(c(81.145,64.288,1.559,1.559,41.943,0.649,0.649), c(104.546,57.61,0.964,0.964,51.49,0.618,0.618), c(160.249,64.219,0.69,0.69,66.387,0.555,0.555)),
      cbind(c(140.804,89.962,1.286,1.286,57.801,0.6,0.6), c(151.094,81.772,0.923,0.923,67.123,0.585,0.585), c(198.167,83.815,0.686,0.686,81.832,0.548,0.548)),
      cbind(c(182.45,103.176,1.028,1.028,72.097,0.567,0.567), c(190.661,99.331,0.843,0.843,81.835,0.56,0.56), c(239.425,103.905,0.688,0.688,97.362,0.543,0.543)),
      cbind(c(217.784,120.787,0.928,0.928,88.611,0.559,0.559), c(231.864,119.546,0.8,0.8,98.278,0.553,0.553), c(276.997,124.613,0.686,0.686,113.145,0.538,0.538)),
      cbind(c(256.53,132.499,0.82,0.82,103.371,0.544,0.544), c(270.451,134.341,0.75,0.75,113.131,0.542,0.542), c(310.982,138.227,0.654,0.654,127.989,0.53,0.53)),
      cbind(c(277.429,143.561,0.75,0.75,117.059,0.53,0.53), c(293.431,144.615,0.685,0.685,126.059,0.525,0.525), c(348.217,154.378,0.638,0.638,140.756,0.518,0.518)))
    col <- switch(case, none = 1L, intercept = 2L, trend = 3L)
    mu <- MU[[j]][, col]; v <- VV[[j]][, col]
  }
  names(mu) <- names(v) <- stat_names
  list(mu = mu, v = v, k = k, case = case)
}

#' @export
#' @method print xtpedroni
print.xtpedroni <- function(x, ...) {
  cat("\nPedroni panel cointegration tests (panel-specific cointegrating vectors)\n")
  cat(sprintf("N = %d panels, T = %d periods (after differencing), k = %d regressor(s)\n", x$N, x$T, x$k))
  cat(sprintf("Deterministics: %s | Bartlett kernel lags: %d | ADF lags: %d\n", if (x$trend) "constant and trend" else "constant", x$kernel_lags, x$adf_lags))
  cat("H0: no cointegration\n\n")
  tab <- data.frame(Statistic = sprintf("%9.4f", x$stats), `p-value` = sprintf("%7.4f", x$p), Tail = unname(x$tail), row.names = names(x$stats), check.names = FALSE)
  print(tab)
  cat("\nStandardized statistics (chi - mu sqrt(N)) / sqrt(nu).\n")
  cat("p-values: 1 - Phi(z) for the variance ratio (rejects for large positive z),\n")
  cat("Phi(z) for the other six statistics (reject for large negative z); Pedroni (1999, p. 668).\n")
  invisible(x)
}

#' @keywords internal
.coint_frame <- function(formula, data, pid_col, tid_col, fname) {
  ## Evaluate the formula with model.frame / model.matrix (so log(), I() and
  ## interactions work), reject missing values with their positions, and
  ## return a data.frame ordered by panel and time with columns pid, tid,
  ## y and the model-matrix columns.
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  yv <- stats::model.response(mf)
  if (is.null(yv)) stop("the formula needs a response")
  Xm <- stats::model.matrix(stats::terms(mf), mf)
  Xm <- Xm[, colnames(Xm) != "(Intercept)", drop = FALSE]
  if (ncol(Xm) < 1L) stop("at least one regressor must be specified")
  na_rows <- which(!stats::complete.cases(cbind(yv, Xm)))
  if (length(na_rows))
    stop(sprintf("%s: %d row(s) have missing values in the variables (rows %s); remove or impute them first",
                 fname, length(na_rows), .fmt_positions(na_rows)))
  if (anyNA(data[[pid_col]]) || anyNA(data[[tid_col]])) stop(fname, ": missing values in the index columns")
  out <- data.frame(data[[pid_col]], data[[tid_col]], as.numeric(yv), stringsAsFactors = FALSE)
  names(out) <- c(pid_col, tid_col, "y")
  for (j in seq_len(ncol(Xm))) out[[colnames(Xm)[j]]] <- as.numeric(Xm[, j])
  out[order(out[[pid_col]], out[[tid_col]]), , drop = FALSE]
}

#' @keywords internal
.fmt_positions <- function(idx, max_show = 20L) {
  if (length(idx) <= max_show) paste(idx, collapse = ", ")
  else paste0(paste(idx[seq_len(max_show)], collapse = ", "), ", ... (", length(idx) - max_show, " more)")
}
