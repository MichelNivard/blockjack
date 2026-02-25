# Baseline from prompt
block_jackknife_lm_baseline <- function(X, y, n_blocks = 200L, separators = NULL) {
  X <- as.matrix(X)
  y <- matrix(as.numeric(y), ncol = 1)
  n <- nrow(X); p <- ncol(X)
  stopifnot(nrow(y) == n, ncol(y) == 1, p <= n)

  if (is.null(separators)) {
    separators <- floor(seq(0, n, length.out = n_blocks + 1L))
    separators[length(separators)] <- n
  }
  n_blocks <- length(separators) - 1L
  stopifnot(n_blocks >= 2L, separators[1] == 0L, separators[length(separators)] == n)

  XtY_blk <- matrix(0, n_blocks, p)
  XtX_blk <- array(0, dim = c(n_blocks, p, p))
  for (b in seq_len(n_blocks)) {
    i0 <- separators[b] + 1L
    i1 <- separators[b + 1L]
    Xb <- X[i0:i1, , drop = FALSE]
    yb <- y[i0:i1, , drop = FALSE]
    XtY_blk[b, ] <- drop(crossprod(Xb, yb))
    XtX_blk[b, , ] <- crossprod(Xb)
  }

  XtY_tot <- colSums(XtY_blk)
  XtX_tot <- apply(XtX_blk, c(2, 3), sum)
  est <- solve(XtX_tot, XtY_tot)

  delete_vals <- matrix(0, n_blocks, p)
  for (b in seq_len(n_blocks)) {
    XtY_del <- XtY_tot - XtY_blk[b, ]
    XtX_del <- XtX_tot - XtX_blk[b, , ]
    delete_vals[b, ] <- solve(XtX_del, XtY_del)
  }

  pseudo <- n_blocks * matrix(est, n_blocks, p, byrow = TRUE) - (n_blocks - 1L) * delete_vals
  cov_jk <- stats::cov(pseudo) / n_blocks
  se_jk <- sqrt(diag(cov_jk))

  list(
    est = est,
    se = se_jk,
    cov = cov_jk,
    delete_values = delete_vals,
    separators = separators
  )
}

# Optimized v1: avoid y matrix conversion, avoid apply over 3d array, avoid repeated matrix(est,...)
block_jackknife_lm_opt1 <- function(X, y, n_blocks = 200L, separators = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(y) == n, p <= n)

  if (is.null(separators)) {
    separators <- floor(seq(0L, n, length.out = n_blocks + 1L))
    separators[length(separators)] <- n
  }
  n_blocks <- length(separators) - 1L
  stopifnot(n_blocks >= 2L, separators[1] == 0L, separators[length(separators)] == n)

  XtY_blk <- matrix(0, n_blocks, p)
  XtX_blk <- array(0, dim = c(n_blocks, p, p))
  XtY_tot <- numeric(p)
  XtX_tot <- matrix(0, p, p)

  for (b in seq_len(n_blocks)) {
    i0 <- separators[b] + 1L
    i1 <- separators[b + 1L]
    idx <- i0:i1
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    xty <- drop(crossprod(Xb, yb))
    xtx <- crossprod(Xb)
    XtY_blk[b, ] <- xty
    XtX_blk[b, , ] <- xtx
    XtY_tot <- XtY_tot + xty
    XtX_tot <- XtX_tot + xtx
  }

  est <- solve(XtX_tot, XtY_tot)

  delete_vals <- matrix(0, n_blocks, p)
  for (b in seq_len(n_blocks)) {
    delete_vals[b, ] <- solve(XtX_tot - XtX_blk[b, , ], XtY_tot - XtY_blk[b, ])
  }

  pseudo <- n_blocks * rep(est, each = n_blocks) - (n_blocks - 1L) * c(delete_vals)
  dim(pseudo) <- c(n_blocks, p)

  # unbiased covariance with denominator (n_blocks - 1), then divide by n_blocks
  pseudo_centered <- sweep(pseudo, 2L, colMeans(pseudo), "-")
  cov_jk <- crossprod(pseudo_centered) / ((n_blocks - 1) * n_blocks)
  se_jk <- sqrt(diag(cov_jk))

  list(
    est = est,
    se = se_jk,
    cov = cov_jk,
    delete_values = delete_vals,
    separators = separators
  )
}

# Optimized v2: no per-block XtX/XtY storage unless requested
block_jackknife_lm_opt2 <- function(X, y, n_blocks = 200L, separators = NULL, keep_delete = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(y) == n, p <= n)

  if (is.null(separators)) {
    separators <- floor(seq(0L, n, length.out = n_blocks + 1L))
    separators[length(separators)] <- n
  }
  n_blocks <- length(separators) - 1L
  stopifnot(n_blocks >= 2L, separators[1] == 0L, separators[length(separators)] == n)

  XtY_blk <- matrix(0, n_blocks, p)
  XtX_blk <- array(0, dim = c(n_blocks, p, p))

  XtY_tot <- numeric(p)
  XtX_tot <- matrix(0, p, p)
  for (b in seq_len(n_blocks)) {
    i0 <- separators[b] + 1L
    i1 <- separators[b + 1L]
    idx <- i0:i1
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    xty <- drop(crossprod(Xb, yb))
    xtx <- crossprod(Xb)

    XtY_blk[b, ] <- xty
    XtX_blk[b, , ] <- xtx
    XtY_tot <- XtY_tot + xty
    XtX_tot <- XtX_tot + xtx
  }

  est <- solve(XtX_tot, XtY_tot)

  if (keep_delete) {
    delete_vals <- matrix(0, n_blocks, p)
  } else {
    delete_vals <- NULL
  }
  pseudo <- matrix(0, n_blocks, p)

  for (b in seq_len(n_blocks)) {
    beta_b <- solve(XtX_tot - XtX_blk[b, , ], XtY_tot - XtY_blk[b, ])
    if (keep_delete) delete_vals[b, ] <- beta_b
    pseudo[b, ] <- n_blocks * est - (n_blocks - 1L) * beta_b
  }

  pseudo_centered <- sweep(pseudo, 2L, colMeans(pseudo), "-")
  cov_jk <- crossprod(pseudo_centered) / ((n_blocks - 1) * n_blocks)
  se_jk <- sqrt(diag(cov_jk))

  list(
    est = est,
    se = se_jk,
    cov = cov_jk,
    delete_values = delete_vals,
    separators = separators
  )
}

# Optimized v3: use Cholesky solves for SPD XtX matrices
block_jackknife_lm_opt3 <- function(X, y, n_blocks = 200L, separators = NULL, keep_delete = TRUE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  stopifnot(length(y) == n, p <= n)

  if (is.null(separators)) {
    separators <- floor(seq(0L, n, length.out = n_blocks + 1L))
    separators[length(separators)] <- n
  }
  n_blocks <- length(separators) - 1L
  stopifnot(n_blocks >= 2L, separators[1] == 0L, separators[length(separators)] == n)

  XtY_blk <- matrix(0, n_blocks, p)
  XtX_blk <- array(0, dim = c(n_blocks, p, p))
  XtY_tot <- numeric(p)
  XtX_tot <- matrix(0, p, p)

  for (b in seq_len(n_blocks)) {
    i0 <- separators[b] + 1L
    i1 <- separators[b + 1L]
    idx <- i0:i1
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]

    xty <- drop(crossprod(Xb, yb))
    xtx <- crossprod(Xb)

    XtY_blk[b, ] <- xty
    XtX_blk[b, , ] <- xtx
    XtY_tot <- XtY_tot + xty
    XtX_tot <- XtX_tot + xtx
  }

  R_tot <- chol(XtX_tot)
  est <- backsolve(R_tot, forwardsolve(t(R_tot), XtY_tot))

  if (keep_delete) {
    delete_vals <- matrix(0, n_blocks, p)
  } else {
    delete_vals <- NULL
  }
  pseudo <- matrix(0, n_blocks, p)

  for (b in seq_len(n_blocks)) {
    XtX_del <- XtX_tot - XtX_blk[b, , ]
    XtY_del <- XtY_tot - XtY_blk[b, ]
    R_del <- chol(XtX_del)
    beta_b <- backsolve(R_del, forwardsolve(t(R_del), XtY_del))
    if (keep_delete) delete_vals[b, ] <- beta_b
    pseudo[b, ] <- n_blocks * est - (n_blocks - 1L) * beta_b
  }

  pseudo_centered <- sweep(pseudo, 2L, colMeans(pseudo), "-")
  cov_jk <- crossprod(pseudo_centered) / ((n_blocks - 1) * n_blocks)
  se_jk <- sqrt(diag(cov_jk))

  list(
    est = est,
    se = se_jk,
    cov = cov_jk,
    delete_values = delete_vals,
    separators = separators
  )
}

bench_once <- function(fun, X, y, n_blocks, reps = 3L, iterations = 1L, ...) {
  times <- numeric(reps)
  out <- NULL
  for (i in seq_len(reps)) {
    gc(FALSE)
    t <- system.time({
      for (j in seq_len(iterations)) {
        out <- fun(X, y, n_blocks = n_blocks, ...)
      }
    })
    times[i] <- unname(t[["elapsed"]]) / iterations
  }
  list(time = times, out = out)
}

run_benchmark <- function(n = 20000L, p = 3L, n_blocks = 200L, reps = 5L, iterations = 20L, seed = 1L) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
  beta_true <- c(0.5, seq(-0.2, by = 0.03, length.out = p - 1L))
  y <- drop(X %*% beta_true + rnorm(n))

  b0 <- bench_once(block_jackknife_lm_baseline, X, y, n_blocks, reps = reps, iterations = iterations)
  b1 <- bench_once(block_jackknife_lm_opt1, X, y, n_blocks, reps = reps, iterations = iterations)
  b2 <- bench_once(block_jackknife_lm_opt2, X, y, n_blocks, reps = reps, iterations = iterations)
  b3 <- bench_once(block_jackknife_lm_opt3, X, y, n_blocks, reps = reps, iterations = iterations)

  tol <- 1e-7
  stopifnot(max(abs(b0$out$est - b1$out$est)) < tol)
  stopifnot(max(abs(b0$out$se - b1$out$se)) < tol)
  stopifnot(max(abs(b0$out$cov - b1$out$cov)) < 1e-7)

  stopifnot(max(abs(b0$out$est - b2$out$est)) < tol)
  stopifnot(max(abs(b0$out$se - b2$out$se)) < tol)
  stopifnot(max(abs(b0$out$cov - b2$out$cov)) < 1e-7)

  stopifnot(max(abs(b0$out$est - b3$out$est)) < tol)
  stopifnot(max(abs(b0$out$se - b3$out$se)) < tol)
  stopifnot(max(abs(b0$out$cov - b3$out$cov)) < 1e-7)

  stats <- function(x) c(min = min(x), median = median(x), mean = mean(x))
  res <- rbind(
    baseline = stats(b0$time),
    opt1 = stats(b1$time),
    opt2 = stats(b2$time),
    opt3_chol = stats(b3$time)
  )

  speedup_vs_baseline <- res["baseline", "median"] / res[, "median"]

  list(
    n = n,
    p = p,
    n_blocks = n_blocks,
    reps = reps,
    iterations = iterations,
    timing = res,
    speedup = speedup_vs_baseline,
    example_est_se = cbind(beta_hat = b3$out$est, jk_se = b3$out$se)
  )
}

print_bench <- function(obj) {
  cat(sprintf("n=%d, p=%d, n_blocks=%d, reps=%d, iterations=%d\n\n", obj$n, obj$p, obj$n_blocks, obj$reps, obj$iterations))
  print(round(obj$timing, 4))
  cat("\nSpeedup vs baseline (median elapsed):\n")
  print(round(obj$speedup, 2))
  cat("\nExample estimates and JK SE (opt3_chol):\n")
  print(round(obj$example_est_se, 6))
}

# Classical OLS SE under homoskedastic model, computed from X and y directly.
ols_se_classical <- function(X, y) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  fit <- .lm.fit(X, y)
  rss <- sum(fit$residuals^2)
  sigma2 <- rss / (n - p)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  sqrt(diag(sigma2 * XtX_inv))
}

compare_jk_to_ols <- function(X, y, n_blocks = 200L, tol_abs = 0.005, tol_rel = 0.2) {
  jk <- block_jackknife_lm_opt1(X, y, n_blocks = n_blocks)
  ols_se <- ols_se_classical(X, y)
  abs_diff <- abs(jk$se - ols_se)
  rel_diff <- abs_diff / pmax(ols_se, 1e-12)
  list(
    pass = (max(abs_diff) <= tol_abs) && (max(rel_diff) <= tol_rel),
    max_abs_diff = max(abs_diff),
    max_rel_diff = max(rel_diff),
    jk_se = jk$se,
    ols_se = ols_se
  )
}

derive_jk_ols_tolerance <- function(R = 200L, n = 20000L, p = 3L, n_blocks = 200L, seed = 123L) {
  set.seed(seed)
  max_abs <- numeric(R)
  max_rel <- numeric(R)
  for (i in seq_len(R)) {
    X <- cbind(1, matrix(rnorm(n * (p - 1L)), n, p - 1L))
    beta <- c(0.5, seq(-0.2, by = 0.03, length.out = p - 1L))
    y <- drop(X %*% beta + rnorm(n))
    cmp <- compare_jk_to_ols(X, y, n_blocks = n_blocks, tol_abs = Inf, tol_rel = Inf)
    max_abs[i] <- cmp$max_abs_diff
    max_rel[i] <- cmp$max_rel_diff
  }
  list(
    abs_quantiles = quantile(max_abs, c(0.5, 0.9, 0.95, 0.99, 1)),
    rel_quantiles = quantile(max_rel, c(0.5, 0.9, 0.95, 0.99, 1))
  )
}

if (sys.nframe() == 0L) {
  cat("Benchmark 1: prompt-sized demo\n")
  r1 <- run_benchmark(n = 20000L, p = 3L, n_blocks = 200L, reps = 5L, iterations = 50L, seed = 1L)
  print_bench(r1)

  cat("\nBenchmark 2: larger stress test\n")
  r2 <- run_benchmark(n = 100000L, p = 10L, n_blocks = 200L, reps = 3L, iterations = 10L, seed = 2L)
  print_bench(r2)

  cat("\nJK vs OLS tolerance research (simulation)\n")
  t1 <- derive_jk_ols_tolerance(R = 100L, n = 20000L, p = 3L, n_blocks = 200L, seed = 11L)
  cat("n=20000, p=3 max-abs quantiles:\n"); print(round(t1$abs_quantiles, 6))
  cat("n=20000, p=3 max-rel quantiles:\n"); print(round(t1$rel_quantiles, 4))

  t2 <- derive_jk_ols_tolerance(R = 60L, n = 100000L, p = 10L, n_blocks = 200L, seed = 22L)
  cat("n=100000, p=10 max-abs quantiles:\n"); print(round(t2$abs_quantiles, 6))
  cat("n=100000, p=10 max-rel quantiles:\n"); print(round(t2$rel_quantiles, 4))
}
