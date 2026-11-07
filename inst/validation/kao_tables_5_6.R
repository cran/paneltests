## Replication of Kao (1999, Tables 5 and 6) with paneltests::xtkao()
##
## Design of Kao (1999, section 5.2, pp. 15-16): y_it - alpha_i - beta x_it = z_it,
## a1 y_it - a2 x_it = w_it, z_it = rho z_i,t-1 + e_zit, w_it = w_i,t-1 + e_wit,
## e_wit = phi_it + pi phi_i,t-1, (e_zit, phi_it) normal with variances 1 and
## sigma^2 and covariance theta sigma; a1 = 1, a2 = -1, beta = 2,
## alpha_i ~ U[0, 10]; T + 1000 observations generated and the first 1000
## discarded. Two ADF lags, Bartlett window 5, critical value -1.645.
## Table 5: power, rho = 0.9, sigma = 1, pi = 0, theta = 0.
## Table 6: size, rho = 1, N = T = 50, sigma in (0.25, 1, 4), pi in
## (-0.8, 0, 0.8), theta in (-0.5, 0, 0.5).
## In Kao (1999, Table 5) the row T = 50, N = 2 repeats the row T = 25, N = 10
## digit for digit; it is kept below as printed.
##
## Usage: Rscript kao_tables_5_6.R [replications] [convention]
## Not run during R CMD check (several hours for 1500 replications).

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1500L
convention <- if (length(args) >= 2) args[2] else "kao"

library(paneltests)

gen <- function(N, TT, rho = 1, sigma = 1, pi = 0, theta = 0) {
  do.call(rbind, lapply(seq_len(N), function(i) {
    n <- TT + 1001
    ez <- stats::rnorm(n)
    phi <- theta * sigma * ez + sigma * sqrt(1 - theta^2) * stats::rnorm(n)
    ew <- phi[-1] + pi * phi[-n]; ez <- ez[-1]
    z <- utils::tail(as.numeric(stats::filter(ez, rho, method = "recursive")), TT)
    w <- utils::tail(cumsum(ew), TT)
    a <- stats::runif(1, 0, 10)
    x <- (w - a - z) / 3
    data.frame(id = i, t = seq_len(TT), y = w - x, x = x)
  }))
}

stat_names <- c("Unadjusted modified DF t", "Unadjusted DF t", "Modified DF t", "DF t", "ADF t")
rate <- function(N, TT, ...) {
  rej <- 0
  for (r in seq_len(R)) {
    k <- suppressWarnings(xtkao(gen(N, TT, ...), y ~ x, index = c("id", "t"), lags = 2,
                                kernel_lags = 5, convention = convention))
    rej <- rej + (k$stats[stat_names] < -1.645)
  }
  rej / R
}

## Kao (1999, Table 5), columns DF_rho, DF_t, DF_rho*, DF_t*, ADF
t5 <- read.table(header = TRUE, text = "
T  N  DFrho DFt   DFrho_s DFt_s ADF
10 2  0.296 0.252 0.022   0.073 0.136
10 10 0.286 0.277 0.014   0.090 0.248
10 50 0.626 0.581 0.027   0.195 0.665
25 2  0.451 0.234 0.154   0.091 0.145
25 10 0.569 0.341 0.255   0.163 0.225
25 50 0.966 0.799 0.774   0.507 0.566
50 2  0.569 0.341 0.255   0.163 0.225
50 10 0.936 0.698 0.853   0.561 0.496
50 50 1.000 1.000 1.000   0.996 0.977")
s5 <- t5
for (j in seq_len(nrow(t5))) { set.seed(j); s5[j, 3:7] <- rate(t5$N[j], t5$T[j], rho = 0.9) }

## Kao (1999, Table 6), N = T = 50; one row per (pi, statistic), columns
## sigma = 0.25, 1, 4, each with theta = -0.5, 0, 0.5
v6 <- c(rep(1, 36), .814, .038, .062, 1, .985, .998, .993, .672, .466,
        .100, .102, .102, .104, .105, .104, .104, .104, .101,
        .089, .089, .089, .087, .090, .088, .089, .088, .088,
        .069, .069, .069, .073, .072, .073, .073, .070, .075,
        .072, .072, .072, .077, .090, .088, .074, .076, .073,
        .097, .096, .096, .099, .098, .100, .099, .099, .100,
        .005, .000, .000, .145, .020, .067, .287, .086, .306,
        .014, .000, .000, .113, .028, .064, .206, .078, .215,
        .287, .095, .041, .208, .163, .273, .101, .076, .089,
        .161, .070, .048, .149, .116, .174, .098, .078, .088,
        .519, .615, .416, .157, .317, .269, .049, .113, .044)
t6 <- expand.grid(theta = c(-0.5, 0, 0.5), sigma = c(0.25, 1, 4),
                  stat = c("DFrho", "DFt", "DFrho_s", "DFt_s", "ADF"), pi = c(-0.8, 0, 0.8))
t6$kao <- v6
t6 <- stats::reshape(t6, idvar = c("theta", "sigma", "pi"), timevar = "stat", direction = "wide")
names(t6) <- sub("kao.", "", names(t6), fixed = TRUE)
s6 <- t6
for (j in seq_len(nrow(t6))) {
  set.seed(100 + j)
  s6[j, 4:8] <- rate(50, 50, rho = 1, sigma = t6$sigma[j], pi = t6$pi[j], theta = t6$theta[j])
}

cat(sprintf("\nconvention = \"%s\", %d replications per cell\n", convention, R))
cat("\nTable 5, mean absolute difference from Kao (1999):\n")
print(round(colMeans(abs(s5[, 3:7] - t5[, 3:7])), 3))
cat("\nTable 6, mean absolute difference from Kao (1999):\n")
print(round(colMeans(abs(s6[, 4:8] - t6[, 4:8])), 3))
