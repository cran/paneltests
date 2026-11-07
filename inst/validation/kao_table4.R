## Replication of Kao (1999, Table 4) with paneltests::xtkao()
##
## Design of Kao (1999, section 5.2, pp. 15-16): y_it - alpha_i - beta x_it = z_it,
## a1 y_it - a2 x_it = w_it, z_it = rho z_i,t-1 + e_zit, w_it = w_i,t-1 + e_wit,
## with a1 = 1, a2 = -1, beta = 2, alpha_i ~ U[0, 10], rho = 1 (no cointegration),
## sigma = 1, pi = 0, theta = 0; T + 1000 observations generated and the first
## 1000 discarded. Two ADF lags, Bartlett window of lag length 5, rejection
## when the statistic is below the asymptotic critical value -1.645.
##
## Usage: Rscript kao_table4.R [replications] [convention]
## Not run during R CMD check (about 20 minutes for 5000 replications).

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 5000L
convention <- if (length(args) >= 2) args[2] else "kao"

library(paneltests)

gen <- function(N, TT) {
  do.call(rbind, lapply(seq_len(N), function(i) {
    z <- utils::tail(cumsum(stats::rnorm(TT + 1000)), TT)
    w <- utils::tail(cumsum(stats::rnorm(TT + 1000)), TT)
    a <- stats::runif(1, 0, 10)
    x <- (w - a - z) / 3
    data.frame(id = i, t = seq_len(TT), y = w - x, x = x)
  }))
}

## Kao (1999, Table 4), columns DF_rho, DF_t, DF_rho*, DF_t*, ADF
kao <- read.table(header = TRUE, text = "
T  N   DFrho DFt   DFrho_s DFt_s ADF
10 2   0.257 0.233 0.028   0.079 0.138
10 5   0.182 0.191 0.021   0.084 0.167
10 10  0.149 0.179 0.018   0.096 0.210
10 25  0.149 0.212 0.024   0.121 0.367
10 50  0.178 0.270 0.041   0.192 0.575
25 2   0.307 0.180 0.130   0.092 0.125
25 5   0.195 0.139 0.088   0.089 0.115
25 10  0.142 0.114 0.068   0.079 0.113
25 25  0.115 0.110 0.052   0.078 0.129
25 50  0.119 0.122 0.053   0.089 0.161
50 2   0.322 0.159 0.249   0.115 0.135
50 5   0.202 0.118 0.153   0.093 0.107
50 10  0.149 0.101 0.106   0.079 0.092
50 25  0.113 0.091 0.077   0.074 0.094
50 50  0.105 0.090 0.072   0.075 0.098")

stat_names <- c("Unadjusted modified DF t", "Unadjusted DF t", "Modified DF t", "DF t", "ADF t")
sim <- kao
sim[, 3:7] <- NA_real_
for (j in seq_len(nrow(kao))) {
  set.seed(1000 * kao$T[j] + kao$N[j])
  rej <- 0
  for (r in seq_len(R)) {
    d <- gen(kao$N[j], kao$T[j])
    k <- suppressWarnings(xtkao(d, y ~ x, index = c("id", "t"), lags = 2,
                                kernel_lags = 5, convention = convention))
    rej <- rej + (k$stats[stat_names] < -1.645)
  }
  sim[j, 3:7] <- rej / R
}

cat(sprintf("\nconvention = \"%s\", %d replications per cell\n\n", convention, R))
out <- cbind(kao[, 1:2], kao = kao[, 3:7], sim = round(sim[, 3:7], 4))
print(out, row.names = FALSE)
cat("\nMean absolute difference from Kao (1999, Table 4):\n")
print(round(colMeans(abs(sim[, 3:7] - kao[, 3:7])), 3))
