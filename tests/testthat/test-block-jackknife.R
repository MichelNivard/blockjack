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
