# Validation Checks

This vignette gives short, runnable checks so users can validate that
`blockjack` behaves as expected for OLS, WLS, and clustered data.

## 1. OLS: coefficients match `lm`, SE are close

``` r
set.seed(101)
n <- 20000
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 0.5 + 0.7 * x1 - 0.2 * x2 + rnorm(n)
d <- data.frame(y = y, x1 = x1, x2 = x2)

fit_bj <- bjlm(y ~ x1 + x2, data = d, n_blocks = 200)
fit_lm <- lm(y ~ x1 + x2, data = d)

coef_diff <- max(abs(coef(fit_bj) - coef(fit_lm)))
se_ols <- coef(summary(fit_lm))[, 2]
se_abs_diff <- max(abs(fit_bj$se - se_ols))
se_rel_diff <- max(abs(fit_bj$se - se_ols) / pmax(se_ols, 1e-12))

c(coef_max_diff = coef_diff,
  se_max_abs_diff = se_abs_diff,
  se_max_rel_diff = se_rel_diff)
#>   coef_max_diff se_max_abs_diff se_max_rel_diff 
#>    4.440892e-16    5.927894e-04    8.434721e-02
```

Expected: coefficient difference near machine precision; SE differences
small.

## 2. WLS: weighted coefficients match `lm(..., weights=)`

``` r
set.seed(202)
n <- 20000
x1 <- rnorm(n)
x2 <- rnorm(n)
v <- exp(0.7 * x1)
w <- 1 / v
y <- 0.4 + 0.8 * x1 - 0.3 * x2 + rnorm(n, sd = sqrt(v))
d <- data.frame(y = y, x1 = x1, x2 = x2, w = w)

fit_bj_w <- bjlm(y ~ x1 + x2, data = d, weights = w, n_blocks = 200)
fit_wls <- lm(y ~ x1 + x2, data = d, weights = w)

coef_table <- cbind(blockjack = coef(fit_bj_w),
                    wls = coef(fit_wls),
                    diff = coef(fit_bj_w) - coef(fit_wls))

se_wls <- coef(summary(fit_wls))[, 2]
se_table <- cbind(blockjack_jk_se = fit_bj_w$se,
                  wls_se = se_wls,
                  abs_diff = abs(fit_bj_w$se - se_wls),
                  rel_diff = abs(fit_bj_w$se - se_wls) / pmax(se_wls, 1e-12))

round(coef_table, 10)
#>              blockjack        wls diff
#> (Intercept)  0.3972959  0.3972959    0
#> x1           0.8034756  0.8034756    0
#> x2          -0.3017081 -0.3017081    0
round(se_table, 10)
#>             blockjack_jk_se      wls_se     abs_diff    rel_diff
#> (Intercept)     0.007628293 0.007616547 0.0000117460 0.001542166
#> x1              0.006215225 0.006197613 0.0000176120 0.002841740
#> x2              0.005839611 0.006169317 0.0003297062 0.053442891
```

Expected: coefficients match `lm` very closely; JK SE and WLS model SE
are close, not identical.

## 3. Cluster-aware JK: compare against `sandwich::vcovCL`

When `cluster=` is provided, `bjlm()` keeps clusters intact when forming
jackknife blocks.

``` r
if (!requireNamespace("sandwich", quietly = TRUE)) {
  stop("Please install 'sandwich' to run cluster validation.")
}

set.seed(303)
n <- 20000
cluster_size <- 10
n_clusters <- n / cluster_size
n_runs <- 100

ratio <- matrix(NA_real_, n_runs, 3)
colnames(ratio) <- c("(Intercept)", "x1", "x2")

for (r in seq_len(n_runs)) {
  cl <- rep(seq_len(n_clusters), each = cluster_size)
  g1 <- rep(rnorm(n_clusters), each = cluster_size)
  g2 <- rep(rnorm(n_clusters), each = cluster_size)
  x1 <- sqrt(0.8) * g1 + sqrt(0.2) * rnorm(n)
  x2 <- sqrt(0.8) * g2 + sqrt(0.2) * rnorm(n)
  u <- rep(rnorm(n_clusters, sd = 2), each = cluster_size)
  y <- 0.5 + 0.7 * x1 - 0.2 * x2 + u + rnorm(n)
  d <- data.frame(y = y, x1 = x1, x2 = x2, cl = cl)

  fit_bj <- bjlm(y ~ x1 + x2, data = d, cluster = "cl", n_blocks = 200)
  fit_lm <- lm(y ~ x1 + x2, data = d)
  se_cl <- sqrt(diag(sandwich::vcovCL(fit_lm, cluster = d$cl)))

  ratio[r, ] <- fit_bj$se / se_cl
}

round(apply(ratio, 2, quantile, probs = c(0.05, 0.50, 0.95)), 4)
#>     (Intercept)     x1     x2
#> 5%       0.9142 0.9330 0.9200
#> 50%      1.0071 1.0103 1.0021
#> 95%      1.0811 1.0839 1.0874
```

Expected: ratios near 1 (with sampling variability), showing
cluster-aware JK SE tracks `vcovCL` in this setting.
