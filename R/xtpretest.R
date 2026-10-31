#' Comprehensive Panel Data Pre-Testing Suite
#'
#' Performs a full battery of panel data pre-tests: Hsiao (2014) homogeneity
#' F-tests, robust (HC1) versions, Swamy (1970) parameter heterogeneity test,
#' cross-sectional dependence (Pesaran 2004), and panel summary statistics.
#'
#' @param data A \code{data.frame} in long format containing all variables.
#' @param formula A two-sided formula of the form \code{y ~ x1 + x2 + ...}
#'   specifying the dependent and independent variables. Terms are evaluated
#'   with \code{model.frame} and \code{model.matrix}, so \code{log(x)},
#'   \code{I(x^2)} and interactions are allowed; every column of the model
#'   matrix (without intercept) counts as one regressor in \eqn{K}.
#' @param index Character vector of length 2: \code{c("panel_id", "time_id")}.
#' @param tests Character vector. Which modules to run. Possible values:
#'   \code{"summary"}, \code{"hsiao"}, \code{"robust"}, \code{"heterogeneity"},
#'   \code{"csd"}.  Default \code{"ALL"} runs everything.
#' @param level Numeric. Significance level for decisions (default 0.05).
#'
#' @return A list (invisibly) with components:
#'   \describe{
#'     \item{\code{summary}}{Panel summary statistics.}
#'     \item{\code{hsiao}}{Hsiao homogeneity F-test results.}
#'     \item{\code{robust}}{Robust HC1 Wald versions of the Hsiao F2 and F3
#'       tests (see Details), with the names of any collinear columns dropped.}
#'     \item{\code{swamy}}{Swamy heterogeneity test results.}
#'     \item{\code{csd}}{Cross-sectional dependence test results.}
#'     \item{\code{recommendation}}{Character. Suggested estimator.}
#'     \item{\code{failed}}{Named list of error messages for modules that
#'       failed and were skipped (absent when all modules ran).}
#'     \item{\code{dropped}}{Character vector of regressor columns removed
#'       because they were exactly collinear with the others (absent when
#'       none).}
#'   }
#'
#' @details
#' \strong{Input checks.} Missing values in the response or the regressors
#' are refused with an error that reports their number and row positions;
#' the tests have no defined treatment of missing cells. Fewer than two
#' panels, or any panel with fewer than \eqn{K + 2} observations (so that
#' the individual regression of Hsiao's \eqn{S_1} has no residual degrees of
#' freedom), are refused with an error naming the panels. Exactly collinear
#' regressor columns are detected from the rank of the QR decomposition of
#' the pooled design and removed with a warning naming them; \eqn{K} is the
#' rank, so the degrees of freedom refer to the regressors that remain. A
#' regressor that is constant within one panel is not collinear in the
#' pooled design; the HSIAO and HETEROGENEITY modules then stop with a
#' message naming the panel (recorded in \code{failed}) because the
#' individual regressions are not estimable, while the ROBUST module drops
#' the aliased interaction columns and reports them.
#'
#' \strong{Hsiao's F tests and their labels.} The module HSIAO reports
#' three F statistics labelled \code{F1}, \code{F2} and \code{F3} in the
#' order overall homogeneity, slope homogeneity and intercept homogeneity.
#' In Hsiao (2014, section 2.2.2) these are, respectively, eq. (2.2.14)
#' (there called \eqn{F_3}), eq. (2.2.15) (there called \eqn{F_1}) and eq.
#' (2.2.16) (there called \eqn{F_4}). The formulas and degrees of freedom
#' are those of the book; only the labels differ, and they are kept for
#' backward compatibility of the returned list.
#'
#' \strong{Swamy's test.} The module HETEROGENEITY computes Swamy (1970,
#' eq. 6.2): \eqn{H = \sum_i (b_i - \hat\beta)' X_i' X_i (b_i - \hat\beta) /
#' s_{ii}}, where \eqn{b_i} are the slopes of the individual OLS
#' regressions, \eqn{X_i} the within-panel demeaned regressors, \eqn{s_{ii}
#' = RSS_i / (T_i - K - 1)}, and \eqn{\hat\beta = [\sum_i X_i'X_i /
#' s_{ii}]^{-1} \sum_i (X_i'X_i / s_{ii}) b_i} is the weighted pooled
#' estimator of the paper. Versions up to 1.0.7 centred the statistic on the
#' unweighted fixed-effects estimator, which over-rejects under cross-panel
#' heteroskedasticity. \eqn{H} is referred to a chi-square with \eqn{K(N -
#' 1)} degrees of freedom.
#'
#' \strong{Cross-sectional dependence.} The module CSD computes Pesaran's
#' CD statistic \eqn{\sqrt{2/(N(N-1))} \sum_{i<j} \sqrt{T_{ij}}
#' \hat\rho_{ij}} from the residuals of the fixed-effects (within)
#' regression with common slopes, using the pairwise-complete \eqn{T_{ij}}
#' overlapping periods for each pair, with a two-sided normal p-value. The
#' residuals of N individual OLS regressions would give a different value;
#' the fixed-effects residuals are the package's choice.
#'
#' \strong{ROBUST module and its size.} The \code{"robust"} module is an extension of this package and is not
#' part of Hsiao (2014), who only notes that the F-tests assume i.i.d.
#' errors. It fits the fully interacted regression
#' \eqn{y = a + Xb + \sum_{i>1} D_i g_i + \sum_{i>1} (D_i X) h_i + u},
#' with the first panel as reference, and computes HC1 Wald statistics
#' (divided by the number of restrictions and referred to an F distribution)
#' for slope homogeneity (\eqn{h_i = 0}, \eqn{(N-1)K} restrictions) and for
#' intercept homogeneity given heterogeneous slopes (\eqn{g_i = 0},
#' \eqn{N-1} restrictions). Columns aliased in the interacted regression
#' (for example a regressor that is constant within a panel) are dropped and
#' reported. A failure in any module no longer aborts the other modules.
#'
#' The HC1 Wald statistic has \eqn{(N - 1)K} restrictions and \eqn{(K + 1)N}
#' parameters, and it over-rejects a true null at small \eqn{T}. In the
#' simulations of the 1.0.8 validation (i.i.d. normal errors, heterogeneous
#' intercepts, homogeneous slopes, \eqn{K = 2}, nominal 5 percent) the
#' rejection rates of the robust slope test were 33.4 percent (\eqn{N = 10},
#' \eqn{T = 20}; 8000 replications, Monte Carlo standard error 0.5
#' percentage points), 64.8 percent (\eqn{N = 25}, \eqn{T = 20}), 15.4
#' percent (\eqn{N = 10}, \eqn{T = 50}) and 11.8 percent (\eqn{N = 25},
#' \eqn{T = 100}); the last three come from 500 replications (standard error
#' about 2 percentage points). The classical F test in the same cells rejects
#' at 3.4 to 6.4 percent.
#' The printed output repeats this note. Use the classical F test (with
#' caution under heteroskedasticity, where it reached 8.4 percent) or a
#' bootstrap of the Wald statistic; the estimator itself is unchanged.
#'
#' @references
#' Hsiao, C. (2014). \emph{Analysis of Panel Data} (3rd ed.).
#' Cambridge University Press.
#' \doi{10.1017/CBO9781139839327}
#'
#' Swamy, P. A. V. B. (1970). Efficient inference in a random coefficient
#' regression model. \emph{Econometrica}, 38(2), 311-323.
#' \doi{10.2307/1913012}
#'
#' Pesaran, M. H. (2004). General diagnostic tests for cross section
#' dependence in panels. \emph{Cambridge Working Paper in Economics}, No. 0435.
#' \doi{10.2139/ssrn.572504}
#'
#' @examples
#' set.seed(10)
#' n <- 5; t <- 10
#' df <- data.frame(
#'   id   = rep(1:n, each = t),
#'   time = rep(1:t, times = n),
#'   y    = rnorm(n * t),
#'   x1   = rnorm(n * t)
#' )
#' res <- xtpretest(df, y ~ x1, index = c("id", "time"),
#'                  tests = c("hsiao", "csd"))
#'
#' @export
xtpretest <- function(data,
                      formula,
                      index,
                      tests = "ALL",
                      level = 0.05) {

  ## ---- Validate --------------------------------------------------------
  if (!inherits(data, "data.frame")) stop("'data' must be a data.frame.")
  if (length(index) != 2L)          stop("'index' must be a character vector of length 2.")
  pid_col <- index[1L]; tid_col <- index[2L]
  if (!all(c(pid_col, tid_col) %in% names(data)))
    stop("index columns not found in 'data'.")

  tests <- toupper(tests)
  if (identical(tests, "ALL")) {
    tests <- c("SUMMARY", "HSIAO", "ROBUST", "HETEROGENEITY", "CSD")
  }

  ## ---- Extract variables -----------------------------------------------
  ## model.frame / model.matrix so that log(), I() and interactions work
  mf   <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  lhs  <- names(mf)[1L]
  y    <- stats::model.response(mf)
  if (is.null(y)) stop("the formula needs a response")
  y    <- as.numeric(y)
  Xm   <- stats::model.matrix(stats::terms(mf), mf)
  Xm   <- Xm[, colnames(Xm) != "(Intercept)", drop = FALSE]
  if (ncol(Xm) < 1L) stop("at least one regressor must be specified")
  pid  <- data[[pid_col]]
  tid  <- data[[tid_col]]

  ## Missing values: the tests define no treatment of missing cells, so
  ## they are refused with their count and positions.
  na_rows <- which(!stats::complete.cases(cbind(y, Xm)))
  if (length(na_rows))
    stop(sprintf("%d row(s) have missing values in the response or regressors (rows %s); remove or impute them before testing",
                 length(na_rows), .xtp_positions(na_rows)))
  if (anyNA(pid) || anyNA(tid)) stop("missing values in the index columns")

  ## Exactly collinear regressors: detected from the QR rank of the pooled
  ## design; the aliased columns are removed and K is the rank.
  qx <- qr(cbind(1, Xm))
  dropped <- character(0)
  if (qx$rank < ncol(Xm) + 1L) {
    bad <- qx$pivot[-seq_len(qx$rank)]
    if (1L %in% bad) stop("the intercept is aliased in the pooled design: the regressors do not vary")
    dropped <- colnames(Xm)[bad - 1L]
    warning(sprintf("%d regressor column(s) exactly collinear with the others were dropped: %s; K is now %d",
                    length(dropped), paste(dropped, collapse = ", "), qx$rank - 1L))
    Xm <- Xm[, setdiff(seq_len(ncol(Xm)), bad - 1L), drop = FALSE]
  }
  rhs  <- colnames(Xm)
  X    <- as.data.frame(Xm, stringsAsFactors = FALSE, optional = TRUE)
  names(X) <- rhs

  panels <- sort(unique(pid))
  N      <- length(panels)
  T_avg  <- nrow(data) / N
  K      <- length(rhs)

  if (N < 2L) stop("at least two panels are required (N = 1): the homogeneity tests compare panels")
  Ti <- table(pid)
  short <- names(Ti)[Ti < K + 2L]
  if (length(short))
    stop(sprintf("%d panel(s) have fewer than K + 2 = %d observations, so the individual regressions have no residual degrees of freedom: %s",
                 length(short), K + 2L, paste(short, collapse = ", ")))

  .xtp_header(lhs, rhs, pid_col, tid_col, N, T_avg, nrow(data))

  out <- list()
  if (length(dropped)) out$dropped <- dropped

  ## A failure inside one module must not abort the others: the module is
  ## skipped with a message, its result is NULL and the failure is recorded
  ## in out$failed.
  run_module <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      message("  MODULE ", name, " failed and was skipped: ", conditionMessage(e))
      out$failed[[name]] <<- conditionMessage(e)
      NULL
    })
  }

  ## ---- MODULE: Summary ------------------------------------------------
  if ("SUMMARY" %in% tests) {
    out$summary <- run_module("SUMMARY", .xtp_summary(cbind(data.frame(y), X), lhs, pid, level))
  }

  ## ---- MODULE: Hsiao --------------------------------------------------
  if ("HSIAO" %in% tests) {
    out$hsiao <- run_module("HSIAO", .xtp_hsiao(y, X, pid, panels, K, level))
  }

  ## ---- MODULE: Robust -------------------------------------------------
  if ("ROBUST" %in% tests) {
    out$robust <- run_module("ROBUST", .xtp_robust(y, X, pid, data, pid_col, rhs, level))
  }

  ## ---- MODULE: Heterogeneity ------------------------------------------
  if ("HETEROGENEITY" %in% tests) {
    out$swamy <- run_module("HETEROGENEITY", .xtp_swamy(y, X, pid, panels, K))
  }

  ## ---- MODULE: CSD ----------------------------------------------------
  if ("CSD" %in% tests) {
    out$csd <- run_module("CSD", .xtp_csd_pesaran(y, X, pid, tid, panels, level))
  }

  ## ---- Recommendation -------------------------------------------------
  out$recommendation <- .xtp_recommend(out, level)
  .xtp_print_summary(out, level)

  invisible(out)
}


## ==========================================================================
## Internal helpers
## ==========================================================================

#' @keywords internal
.xtp_header <- function(lhs, rhs, pid, tid, N, T_avg, NT) {
  line <- strrep("-", 72)
  message(line)
  message("  xtpretest -- Comprehensive Panel Data Pre-Testing Suite v",
          as.character(utils::packageVersion("paneltests")))
  message(line)
  message("  Dependent  : ", lhs)
  message("  Regressors : ", paste(rhs, collapse = ", "))
  message("  Panel id   : ", pid, "  |  Time id: ", tid)
  message("  N panels   : ", N, "  |  Avg T: ", round(T_avg, 1),
          "  |  Obs: ", NT)
  message(line)
}

#' @keywords internal
.xtp_summary <- function(data, lhs, pid, level) {
  ## data: the evaluated response (first column) and model-matrix columns
  names(data)[1L] <- lhs
  vars <- names(data)
  out  <- list()

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Panel Summary Statistics")
  message(line)
  message(sprintf("  %-16s %-10s %-10s %-10s %-10s %-10s",
                  "Variable", "Overall SD", "Between SD", "Within SD", "Min", "Max"))
  message(line)

  for (v in vars) {
    x       <- data[[v]]
    o_sd    <- stats::sd(x, na.rm = TRUE)
    pmeans  <- tapply(x, pid, mean, na.rm = TRUE)
    b_sd    <- stats::sd(pmeans, na.rm = TRUE)
    grand   <- mean(x, na.rm = TRUE)
    within  <- x - pmeans[match(pid, names(pmeans))] + grand
    w_sd    <- stats::sd(within, na.rm = TRUE)
    message(sprintf("  %-16s %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f",
                    utils::head(v, 1L),
                    o_sd, b_sd, w_sd,
                    min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
    out[[v]] <- list(overall_sd = o_sd, between_sd = b_sd, within_sd = w_sd)
  }
  message(line)
  out
}

#' @keywords internal
.xtp_hsiao <- function(y, X, pid, panels, K, level) {
  N <- length(panels)
  ## S1: unrestricted (individual OLS)
  S1 <- 0; df1 <- 0
  for (p in panels) {
    idx  <- pid == p
    yi   <- y[idx]; Xi <- X[idx, , drop = FALSE]
    mat  <- cbind(1, as.matrix(Xi))
    fit  <- stats::lm.fit(mat, yi)
    if (fit$rank < ncol(mat))
      stop(sprintf("the individual regression of panel %s is rank deficient (a regressor is constant or collinear within the panel); S1 and its degrees of freedom are undefined", as.character(p)))
    S1  <- S1  + sum(fit$residuals^2)
    df1 <- df1 + (sum(idx) - K - 1L)
  }

  ## S2: FE (common slopes, different intercepts)
  fe_mat <- cbind(stats::model.matrix(~ factor(pid) - 1), as.matrix(X))
  fit2   <- stats::lm.fit(fe_mat, y)
  if (fit2$rank < ncol(fe_mat)) stop("the fixed-effects design is rank deficient (a regressor is constant within every panel)")
  S2     <- sum(fit2$residuals^2)
  df2    <- length(y) - ncol(fe_mat)

  ## S3: pooled OLS
  mat3 <- cbind(1, as.matrix(X))
  fit3 <- stats::lm.fit(mat3, y)
  S3   <- sum(fit3$residuals^2)
  df3  <- length(y) - ncol(mat3)

  ## F-statistics
  ## Labels: F1, F2, F3 here are Hsiao's (2.2.14) F3, (2.2.15) F1 and
  ## (2.2.16) F4 respectively (see the Details of the man page).
  F1_nd <- (N - 1L) * (K + 1L); F1_dd <- df1
  F1  <- ((S3 - S1) / F1_nd) / (S1 / F1_dd)
  F1p <- stats::pf(F1, F1_nd, F1_dd, lower.tail = FALSE)

  F2_nd <- (N - 1L) * K; F2_dd <- df1
  F2  <- ((S2 - S1) / F2_nd) / (S1 / F2_dd)
  F2p <- stats::pf(F2, F2_nd, F2_dd, lower.tail = FALSE)

  F3_nd <- N - 1L; F3_dd <- df2
  F3  <- ((S3 - S2) / F3_nd) / (S2 / F3_dd)
  F3p <- stats::pf(F3, F3_nd, F3_dd, lower.tail = FALSE)

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Hsiao (2014) Homogeneity Tests")
  message(line)
  message("  RSS Unrestricted (S1): ", formatC(S1, digits = 4, format = "f"))
  message("  RSS FE (S2):           ", formatC(S2, digits = 4, format = "f"))
  message("  RSS Pooled (S3):       ", formatC(S3, digits = 4, format = "f"))
  message(line)
  .xtp_print_ftest("F1: Overall homogeneity",   F1, F1_nd, F1_dd, F1p, level)
  .xtp_print_ftest("F2: Slope homogeneity",     F2, F2_nd, F2_dd, F2p, level)
  .xtp_print_ftest("F3: Intercept homogeneity", F3, F3_nd, F3_dd, F3p, level)
  message("  (F1, F2, F3 are Hsiao's eqs. 2.2.14, 2.2.15, 2.2.16, there named F3, F1, F4)")
  message(line)

  list(F1 = F1, F1p = F1p, F1_nd = F1_nd, F1_dd = F1_dd,
       F2 = F2, F2p = F2p, F2_nd = F2_nd, F2_dd = F2_dd,
       F3 = F3, F3p = F3p, F3_nd = F3_nd, F3_dd = F3_dd,
       S1 = S1, S2 = S2, S3 = S3)
}

#' @keywords internal
.xtp_print_ftest <- function(label, F, df1, df2, pval, level) {
  if (is.na(F) || is.na(pval)) {
    message(sprintf("  %-38s not computed (NA)", label))
    return(invisible(NULL))
  }
  sig <- if (pval < level) " [REJECT]" else " [accept]"
  message(sprintf("  %-38s F=%8.4f  df=(%d,%d)  p=%6.4f%s",
                  label, F, df1, df2, pval, sig))
}

#' @keywords internal
.xtp_robust <- function(y, X, pid, data, pid_col, rhs, level) {
  ## Heteroskedasticity-robust (HC1) Wald versions of the Hsiao F2 and F3
  ## tests. These are an extension of the package, not part of Hsiao (2014),
  ## who notes only that the F-tests assume i.i.d. errors (3rd ed., p. 24).
  ##
  ## Parameterisation (deviations from a reference panel, the first one):
  ##   y = a + X b + sum_{i>1} D_i g_i + sum_{i>1} (D_i * X) h_i + u
  ## Slope homogeneity    H0: h_i = 0 for all i > 1, (N - 1) K restrictions
  ## Intercept homogeneity H0: g_i = 0 for all i > 1, (N - 1) restrictions,
  ##   conditional on heterogeneous slopes.
  ## Versions <= 1.0.6 included the interactions of ALL N panels together
  ## with X itself; the columns were then exactly collinear (the N
  ## interactions of a regressor sum to the regressor) and solve() failed
  ## with a singular matrix on every data set.
  Xm   <- as.matrix(X)
  K    <- ncol(Xm)
  pf   <- factor(pid)
  N    <- nlevels(pf)
  Dall <- stats::model.matrix(~ pf - 1)
  D    <- Dall[, -1L, drop = FALSE]                       # N - 1 dummies
  colnames(D) <- paste0("D_", levels(pf)[-1L])
  int_terms <- do.call(cbind, lapply(seq_len(K), function(j) {
    M <- D * Xm[, j]
    colnames(M) <- paste0(colnames(D), ":", colnames(Xm)[j])
    M
  }))
  full_mat <- cbind("(Intercept)" = 1, Xm, D, int_terms)
  slope_idx <- ncol(Xm) + 1L + ncol(D) + seq_len(ncol(int_terms))
  int_idx   <- ncol(Xm) + 1L + seq_len(ncol(D))

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Robust (HC1) Homogeneity Tests")
  message(line)

  fit_full <- tryCatch(stats::lm.fit(full_mat, y), error = function(e) NULL)
  if (is.null(fit_full)) {
    message("  Robust tests skipped: the interacted regression could not be fitted.")
    message(line)
    return(NULL)
  }
  ## Columns aliased by lm.fit (NA coefficients) are dropped from the
  ## restriction sets and reported, instead of producing a singular matrix.
  aliased <- which(is.na(fit_full$coefficients))
  if (length(aliased)) {
    message("  Note: ", length(aliased), " collinear column(s) dropped: ",
            paste(colnames(full_mat)[aliased], collapse = ", "))
  }
  keep <- setdiff(seq_len(ncol(full_mat)), aliased)
  hc1  <- .xtp_hc1(full_mat[, keep, drop = FALSE], fit_full$residuals)
  if (is.null(hc1)) {
    message("  Robust tests skipped: the HC1 covariance matrix is singular.")
    message(line)
    return(NULL)
  }
  pos  <- match(seq_len(ncol(full_mat)), keep)      # new position of each column
  b_all <- fit_full$coefficients[keep]
  df2   <- length(y) - length(keep)

  rF2 <- .xtp_wald(b_all, hc1, stats::na.omit(pos[slope_idx]), df2)
  rF3 <- .xtp_wald(b_all, hc1, stats::na.omit(pos[int_idx]),   df2)

  if (!is.null(rF2)) .xtp_print_ftest("Robust F2: Slope homogeneity", rF2$F, rF2$df1, rF2$df2, rF2$p, level)
  if (!is.null(rF3)) .xtp_print_ftest("Robust F3: Intercept homogeneity", rF3$F, rF3$df1, rF3$df2, rF3$p, level)
  message("  Note: the HC1 Wald test over-rejects at small T (simulated size at nominal 5%:")
  message("  33.4% for N=10, T=20 (8000 replications); 64.8% for N=25, T=20; 15.4% for N=10, T=50;")
  message("  11.8% for N=25, T=100 (500 replications each)).")
  message("  Prefer the classical F test (with caution under heteroskedasticity) or a bootstrap.")
  message(line)

  list(rF2 = rF2, rF3 = rF3, aliased = colnames(full_mat)[aliased])
}

#' @keywords internal
.xtp_hc1 <- function(X, resid) {
  ## HC1 sandwich with a QR-based inverse; returns NULL when X'X is singular
  n <- nrow(X); k <- ncol(X)
  qx <- qr(X)
  if (qx$rank < k) return(NULL)
  XtX_inv <- chol2inv(qr.R(qx))
  piv <- qx$pivot
  XtX_inv[piv, piv] <- XtX_inv
  e2 <- resid^2 * (n / (n - k))
  Xe <- X * e2
  XtX_inv %*% crossprod(Xe, X) %*% XtX_inv
}

#' @keywords internal
.xtp_wald <- function(b_all, vcov, idx, df2) {
  if (!length(idx)) return(NULL)
  b   <- b_all[idx]
  V   <- vcov[idx, idx, drop = FALSE]
  df1 <- length(idx)
  Wst <- tryCatch(as.numeric(t(b) %*% solve(V) %*% b) / df1,
                  error = function(e) NA_real_)
  pval <- if (!is.na(Wst)) stats::pf(Wst, df1, df2, lower.tail = FALSE) else NA_real_
  list(F = Wst, df1 = df1, df2 = df2, p = pval)
}

#' @keywords internal
.xtp_swamy <- function(y, X, pid, panels, K) {
  N <- length(panels)
  ## Swamy (1970), eq. (6.2):
  ##   H = sum_i (b_i - beta_hat)' (X_i'X_i / s_ii) (b_i - beta_hat),
  ##   beta_hat = [sum_i X_i'X_i / s_ii]^-1 sum_i (X_i'X_i / s_ii) b_i,
  ## with X_i the within-panel demeaned regressors (Frisch-Waugh: the
  ## inverse of the slope block of s_ii (X_i'X_i)^-1 for the regression with
  ## intercept equals the demeaned X_i'X_i / s_ii), s_ii = RSS_i / (T_i - K - 1).
  ## Versions up to 1.0.7 centred on the unweighted fixed-effects estimator.
  B <- matrix(NA_real_, N, K); W <- vector("list", N)
  for (j in seq_along(panels)) {
    idx <- pid == panels[j]
    yi  <- y[idx]; Xi <- as.matrix(X[idx, , drop = FALSE])
    mat <- cbind(1, Xi)
    fit <- stats::lm.fit(mat, yi)
    if (fit$rank < ncol(mat))
      stop(sprintf("the individual regression of panel %s is rank deficient (a regressor is constant or collinear within the panel)", as.character(panels[j])))
    s_ii <- sum(fit$residuals^2) / (sum(idx) - K - 1L)
    Xd   <- sweep(Xi, 2L, colMeans(Xi))
    B[j, ] <- fit$coefficients[-1L]
    W[[j]] <- crossprod(Xd) / s_ii
  }
  Wsum <- Reduce(`+`, W)
  beta_hat <- solve(Wsum, Reduce(`+`, lapply(seq_len(N), function(j) W[[j]] %*% B[j, ])))
  chi2 <- 0
  for (j in seq_len(N)) {
    d <- B[j, ] - as.numeric(beta_hat)
    chi2 <- chi2 + as.numeric(t(d) %*% W[[j]] %*% d)
  }
  df   <- K * (N - 1L)
  pval <- stats::pchisq(chi2, df = df, lower.tail = FALSE)

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Swamy (1970) Parameter Heterogeneity")
  message(line)
  message(sprintf("  chi2(%d) = %.4f   p-value = %.4f  %s",
                  df, chi2, pval,
                  if (pval < 0.05) "[REJECT: heterogeneity]" else "[accept: homogeneity]"))
  message(line)

  list(chi2 = chi2, df = df, pval = pval, beta_hat = as.numeric(beta_hat))
}

#' @keywords internal
.xtp_csd_pesaran <- function(y, X, pid, tid, panels, level) {
  ## Pesaran (2004) CD statistic computed from the residuals of the
  ## fixed-effects (within) regression with common slopes (package choice,
  ## stated in the man page); pairwise-complete overlapping periods T_ij.
  fe_mat <- cbind(stats::model.matrix(~ factor(pid) - 1), as.matrix(X))
  fit_fe <- tryCatch(stats::lm.fit(fe_mat, y), error = function(e) NULL)
  if (is.null(fit_fe)) {
    message("  MODULE: CSD skipped.")
    return(NULL)
  }
  resid <- fit_fe$residuals
  N     <- length(panels)
  tvals <- sort(unique(tid))
  T_len <- length(tvals)

  ## Build residual matrix T x N
  R <- matrix(NA_real_, nrow = T_len, ncol = N)
  for (j in seq_along(panels)) {
    idx <- pid == panels[j]
    ti  <- match(tid[idx], tvals)
    R[ti, j] <- resid[idx]
  }

  ## Pairwise correlations
  cd_sum  <- 0; rho_sum <- 0; npairs <- 0
  for (i in seq_len(N - 1L)) {
    for (j in (i + 1L):N) {
      valid <- !is.na(R[, i]) & !is.na(R[, j])
      Tij   <- sum(valid)
      if (Tij < 3L) next
      rij   <- stats::cor(R[valid, i], R[valid, j])
      cd_sum  <- cd_sum  + sqrt(Tij) * rij
      rho_sum <- rho_sum + abs(rij)
      npairs  <- npairs + 1L
    }
  }
  cd_stat <- sqrt(2 / (N * (N - 1L))) * cd_sum
  pval    <- 2 * stats::pnorm(abs(cd_stat), lower.tail = FALSE)
  avg_rho <- rho_sum / max(npairs, 1L)

  line <- strrep("-", 72)
  message(line)
  message("  MODULE: Cross-Sectional Dependence (Pesaran 2004)")
  message(line)
  message(sprintf("  CD statistic   = %.4f", cd_stat))
  message(sprintf("  p-value        = %.4f", pval))
  message(sprintf("  Avg |rho_ij|   = %.4f", avg_rho))
  if (pval < level) {
    message("  Cross-sectional dependence detected. Consider CCE/CCEMG.")
  } else {
    message("  No significant cross-sectional dependence.")
  }
  message(line)

  list(cd_stat = cd_stat, pval = pval, avg_rho = avg_rho)
}

#' @keywords internal
.xtp_recommend <- function(out, level) {
  rec <- "Pooled OLS"

  if (!is.null(out$hsiao)) {
    if (!is.na(out$hsiao$F2p) && out$hsiao$F2p < level) rec <- "Mean Group (MG) / CCEMG"
    else if (!is.na(out$hsiao$F3p) && out$hsiao$F3p < level) rec <- "Fixed Effects (FE)"
  }
  if (!is.null(out$swamy)) {
    if (!is.na(out$swamy$pval) && out$swamy$pval < level && rec == "Pooled OLS")
      rec <- "Fixed Effects (FE)"
  }
  if (!is.null(out$csd)) {
    if (!is.na(out$csd$pval) && out$csd$pval < level) {
      if (rec == "Mean Group (MG) / CCEMG") rec <- "CCEMG (cross-sectional dependence)"
      else if (rec == "Fixed Effects (FE)") rec <- "FE with Driscoll-Kraay SE"
    }
  }
  rec
}

#' @keywords internal
.xtp_print_summary <- function(out, level) {
  line <- strrep("-", 72)
  message(line)
  message("  OVERALL PRE-TEST SUMMARY")
  message(line)
  if (!is.null(out$hsiao)) {
    message("  Hsiao F1 (overall)   : p=", round(out$hsiao$F1p, 3L))
    message("  Hsiao F2 (slopes)    : p=", round(out$hsiao$F2p, 3L))
    message("  Hsiao F3 (intercepts): p=", round(out$hsiao$F3p, 3L))
  }
  if (!is.null(out$swamy)) {
    message("  Swamy test           : p=", round(out$swamy$pval, 3L))
  }
  if (!is.null(out$csd)) {
    message("  CSD (Pesaran)        : p=", round(out$csd$pval, 3L))
  }
  message("  >>> Recommended estimator: ", out$recommendation)
  message(line)
}

#' @keywords internal
.xtp_positions <- function(idx, max_show = 20L) {
  if (length(idx) <= max_show) paste(idx, collapse = ", ")
  else paste0(paste(idx[seq_len(max_show)], collapse = ", "), ", ... (", length(idx) - max_show, " more)")
}
