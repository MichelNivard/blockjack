ols_se_classical <- function(X, y) {
  n <- nrow(X)
  p <- ncol(X)
  fit <- .lm.fit(X, y)
  rss <- sum(fit$residuals^2)
  sigma2 <- rss / (n - p)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  sqrt(diag(sigma2 * XtX_inv))
}

test_that("coefficients are close to truth in simulation", {
  set.seed(101)
  n <- 20000L
  p <- 3L
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  beta_true <- c(0.5, -0.2, 0.1)
  y <- drop(X %*% beta_true + rnorm(n))

  fit <- block_jackknife_fit(X, y, n_blocks = 200L)
  expect_lt(max(abs(fit$coefficients - beta_true)), 0.05)
})

test_that("jackknife SE is close to classical OLS SE", {
  set.seed(202)
  n <- 20000L
  p <- 3L
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  beta_true <- c(0.5, -0.2, 0.1)
  y <- drop(X %*% beta_true + rnorm(n))

  fit <- block_jackknife_fit(X, y, n_blocks = 200L)
  ols_se <- ols_se_classical(X, y)

  abs_diff <- abs(fit$se - ols_se)
  rel_diff <- abs_diff / pmax(ols_se, 1e-12)

  expect_lte(max(abs_diff), 0.002)
  expect_lte(max(rel_diff), 0.20)
})

test_that("formula interface and summary work", {
  set.seed(303)
  n <- 15000L
  d <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  d$y <- 0.3 + 0.5 * d$x1 - 0.2 * d$x2 + d$y

  fit <- bjlm(y ~ x1 + x2, data = d, n_blocks = 100L)
  s <- summary(fit)

  expect_s3_class(fit, "bjlm")
  expect_s3_class(s, "summary.bjlm")
  expect_equal(colnames(s$coefficients), c("Estimate", "Jackknife SE", "z value", "Pr(>|z|)"))
  expect_equal(fit$backend, "Rcpp")
  expect_equal(s$backend, "Rcpp")

  lm_fit <- lm(y ~ x1 + x2, data = d)
  expect_lt(max(abs(coef(fit) - coef(lm_fit))), 1e-8)
})

test_that("weighted bjlm matches WLS coefficients and has close SE", {
  set.seed(404)
  n <- 20000L
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  v <- exp(0.7 * x1)
  w <- 1 / v
  e <- rnorm(n, sd = sqrt(v))
  y <- 0.4 + 0.8 * x1 - 0.3 * x2 + e
  d <- data.frame(y = y, x1 = x1, x2 = x2, w = w)

  fit_bj <- bjlm(y ~ x1 + x2, data = d, weights = w, n_blocks = 200L)
  fit_wls <- lm(y ~ x1 + x2, data = d, weights = w)

  expect_lt(max(abs(coef(fit_bj) - coef(fit_wls))), 1e-8)

  se_wls <- coef(summary(fit_wls))[, 2]
  abs_diff <- abs(fit_bj$se - se_wls)
  rel_diff <- abs_diff / pmax(se_wls, 1e-12)

  expect_lte(max(abs_diff), 0.002)
  expect_lte(max(rel_diff), 0.20)
})

test_that("cluster option keeps clusters intact and triggers warnings", {
  set.seed(505)
  n <- 100L
  cluster <- c(rep("A", 40), rep("B", 30), rep("C", 30))
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- 0.2 + 0.3 * x1 - 0.1 * x2 + rnorm(n)
  d <- data.frame(y = y, x1 = x1, x2 = x2, cl = cluster)

  expect_warning(
    fit <- bjlm(y ~ x1 + x2, data = d, cluster = "cl", n_blocks = 10L),
    "Number of clusters is less than n_blocks"
  )
  expect_warning(
    bjlm(y ~ x1 + x2, data = d, cluster = "cl", n_blocks = 10L),
    "Largest cluster is larger than n / n_blocks"
  )

  expect_equal(length(unique(fit$block_id)), 3L)
  expect_equal(length(unique(tapply(fit$block_id, d$cl, function(z) length(unique(z))))), 1L)
})

test_that("cluster-jackknife SE tracks vcovCL under strong cluster correlation", {
  skip_if_not_installed("sandwich")

  set.seed(606)
  n <- 20000L
  cluster_size <- 10L
  n_clusters <- n / cluster_size
  n_runs <- 100L
  p <- 3L

  ratio <- matrix(NA_real_, n_runs, p)
  colnames(ratio) <- c("(Intercept)", "x1", "x2")

  for (r in seq_len(n_runs)) {
    cl <- rep(seq_len(n_clusters), each = cluster_size)
    g1 <- rep(rnorm(n_clusters), each = cluster_size)
    g2 <- rep(rnorm(n_clusters), each = cluster_size)
    x1 <- sqrt(0.8) * g1 + sqrt(0.2) * rnorm(n)
    x2 <- sqrt(0.8) * g2 + sqrt(0.2) * rnorm(n)
    u <- rep(rnorm(n_clusters, sd = 2), each = cluster_size)
    e <- rnorm(n, sd = 1)
    y <- 0.5 + 0.7 * x1 - 0.2 * x2 + u + e
    d <- data.frame(y = y, x1 = x1, x2 = x2, cl = cl)

    fit_bj <- bjlm(y ~ x1 + x2, data = d, cluster = "cl", n_blocks = 200L)
    fit_lm <- lm(y ~ x1 + x2, data = d)
    se_cl <- sqrt(diag(sandwich::vcovCL(fit_lm, cluster = d$cl)))
    se_ols <- coef(summary(fit_lm))[, 2]

    ratio[r, ] <- fit_bj$se / se_cl
    expect_true(median(se_cl / se_ols) > 1.5)
  }

  q <- apply(ratio, 2, quantile, probs = c(0.05, 0.50, 0.95))
  expect_true(all(is.finite(q)))
  expect_true(all(q > 0.5))
  expect_true(all(q < 1.8))
})

test_that("OLS JK/model-SE ratio quantiles are stable across 200 runs", {
  set.seed(911)
  n <- 20000L
  n_runs <- 200L
  ratio <- matrix(NA_real_, n_runs, 3L)
  colnames(ratio) <- c("(Intercept)", "x1", "x2")

  for (r in seq_len(n_runs)) {
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    d <- data.frame(
      y = 0.5 - 0.2 * x1 + 0.1 * x2 + rnorm(n),
      x1 = x1,
      x2 = x2
    )
    fit_bj <- bjlm(y ~ x1 + x2, data = d, n_blocks = 200L)
    fit_lm <- lm(y ~ x1 + x2, data = d)
    se_lm <- coef(summary(fit_lm))[, 2]
    ratio[r, ] <- fit_bj$se / se_lm
  }

  q <- apply(ratio, 2, quantile, probs = c(0.05, 0.50, 0.95))
  expect_true(all(is.finite(q)))
  expect_true(all(q[1, ] > 0.85))
  expect_true(all(q[2, ] > 0.95 & q[2, ] < 1.05))
  expect_true(all(q[3, ] < 1.15))
})

test_that("WLS JK/model-SE ratio quantiles are stable across 200 runs", {
  set.seed(912)
  n <- 20000L
  n_runs <- 200L
  ratio <- matrix(NA_real_, n_runs, 3L)
  colnames(ratio) <- c("(Intercept)", "x1", "x2")

  for (r in seq_len(n_runs)) {
    x1 <- rnorm(n)
    x2 <- rnorm(n)
    v <- exp(0.7 * x1)
    w <- 1 / v
    y <- 0.4 + 0.8 * x1 - 0.3 * x2 + rnorm(n, sd = sqrt(v))
    d <- data.frame(y = y, x1 = x1, x2 = x2, w = w)
    fit_bj <- bjlm(y ~ x1 + x2, data = d, weights = w, n_blocks = 200L)
    fit_wls <- lm(y ~ x1 + x2, data = d, weights = w)
    se_wls <- coef(summary(fit_wls))[, 2]
    ratio[r, ] <- fit_bj$se / se_wls
  }

  q <- apply(ratio, 2, quantile, probs = c(0.05, 0.50, 0.95))
  expect_true(all(is.finite(q)))
  expect_true(all(q[1, ] > 0.85))
  expect_true(all(q[2, ] > 0.95 & q[2, ] < 1.05))
  expect_true(all(q[3, ] < 1.15))
})

test_that("Rcpp and R backends are numerically equivalent (OLS/WLS/cluster)", {
  set.seed(707)
  n <- 12000L

  # OLS
  X1 <- cbind(1, matrix(rnorm(n * 2L), n, 2L))
  y1 <- drop(X1 %*% c(0.5, -0.2, 0.1) + rnorm(n))
  fR1 <- block_jackknife_fit(X1, y1, n_blocks = 100L, backend = "R")
  fC1 <- block_jackknife_fit(X1, y1, n_blocks = 100L, backend = "Rcpp")
  expect_lte(max(abs(fR1$coefficients - fC1$coefficients)), 1e-10)
  expect_lte(max(abs(fR1$se - fC1$se)), 1e-10)
  expect_lte(max(abs(fR1$cov - fC1$cov)), 1e-10)

  # WLS
  x1 <- rnorm(n); x2 <- rnorm(n)
  X2 <- cbind(1, x1, x2)
  v <- exp(0.7 * x1)
  w <- 1 / v
  y2 <- drop(X2 %*% c(0.4, 0.8, -0.3) + rnorm(n, sd = sqrt(v)))
  fR2 <- block_jackknife_fit(X2, y2, n_blocks = 100L, weights = w, backend = "R")
  fC2 <- block_jackknife_fit(X2, y2, n_blocks = 100L, weights = w, backend = "Rcpp")
  expect_lte(max(abs(fR2$coefficients - fC2$coefficients)), 1e-10)
  expect_lte(max(abs(fR2$se - fC2$se)), 1e-10)
  expect_lte(max(abs(fR2$cov - fC2$cov)), 1e-10)

  # Cluster/block-id
  cluster_size <- 10L
  n_clusters <- n / cluster_size
  cl <- rep(seq_len(n_clusters), each = cluster_size)
  block_id <- ((cl - 1L) %% 100L) + 1L
  g1 <- rep(rnorm(n_clusters), each = cluster_size)
  g2 <- rep(rnorm(n_clusters), each = cluster_size)
  X3 <- cbind(1, sqrt(0.8) * g1 + sqrt(0.2) * rnorm(n), sqrt(0.8) * g2 + sqrt(0.2) * rnorm(n))
  u <- rep(rnorm(n_clusters, sd = 2), each = cluster_size)
  y3 <- drop(X3 %*% c(0.5, 0.7, -0.2) + u + rnorm(n))
  fR3 <- block_jackknife_fit(X3, y3, block_id = block_id, backend = "R")
  fC3 <- block_jackknife_fit(X3, y3, block_id = block_id, backend = "Rcpp")
  expect_lte(max(abs(fR3$coefficients - fC3$coefficients)), 1e-10)
  expect_lte(max(abs(fR3$se - fC3$se)), 1e-10)
  expect_lte(max(abs(fR3$cov - fC3$cov)), 1e-10)
})

test_that("Rcpp backend is not slower than R backend on representative workloads", {
  skip_on_cran()

  bench <- function(expr_fun, reps = 4L, iterations = 15L) {
    t <- numeric(reps)
    for (i in seq_len(reps)) {
      gc(FALSE)
      tt <- system.time({
        for (j in seq_len(iterations)) expr_fun()
      })
      t[i] <- unname(tt[["elapsed"]]) / iterations
    }
    median(t)
  }

  set.seed(808)
  n <- 12000L
  x1 <- rnorm(n); x2 <- rnorm(n)
  X <- cbind(1, x1, x2)
  y <- drop(X %*% c(0.5, -0.2, 0.1) + rnorm(n))
  v <- exp(0.5 * x1); w <- 1 / v
  cl <- rep(seq_len(n / 10L), each = 10L)
  block_id <- ((cl - 1L) %% 100L) + 1L

  med_R_ols <- bench(function() block_jackknife_fit(X, y, n_blocks = 100L, backend = "R"))
  med_C_ols <- bench(function() block_jackknife_fit(X, y, n_blocks = 100L, backend = "Rcpp"))
  expect_lte(med_C_ols, med_R_ols * 1.10)

  med_R_wls <- bench(function() block_jackknife_fit(X, y, n_blocks = 100L, weights = w, backend = "R"))
  med_C_wls <- bench(function() block_jackknife_fit(X, y, n_blocks = 100L, weights = w, backend = "Rcpp"))
  expect_lte(med_C_wls, med_R_wls * 1.10)

  med_R_cl <- bench(function() block_jackknife_fit(X, y, block_id = block_id, backend = "R"))
  med_C_cl <- bench(function() block_jackknife_fit(X, y, block_id = block_id, backend = "Rcpp"))
  expect_lte(med_C_cl, med_R_cl * 1.10)
})
