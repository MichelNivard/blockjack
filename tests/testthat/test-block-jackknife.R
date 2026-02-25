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
