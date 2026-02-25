# Standalone fast block jackknife function for users.
#
# This script can be sourced directly:
# source("block_jackknife_lm_fast.R")
# fit <- bjlm(X, y, n_blocks = 200)

bjlm <- function(X, y, n_blocks = 200L, separators = NULL) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
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
  pseudo_centered <- sweep(pseudo, 2L, colMeans(pseudo), "-")
  cov_jk <- crossprod(pseudo_centered) / ((n_blocks - 1L) * n_blocks)
  se_jk <- sqrt(diag(cov_jk))

  list(est = est, se = se_jk, cov = cov_jk, delete_values = delete_vals, separators = separators)
}
