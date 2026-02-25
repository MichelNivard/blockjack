# Fast block jackknife for linear regression.
#
# block_jackknife_fit() is the matrix interface.
# bjlm() is a formula/data wrapper with an lm-like summary.

.block_jackknife_prepare <- function(X, y, n_blocks, separators, weights, block_id) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  n <- nrow(X)
  p <- ncol(X)
  stopifnot(length(y) == n, p <= n)

  if (is.null(weights)) {
    weights <- rep(1, n)
  } else {
    weights <- as.numeric(weights)
    stopifnot(length(weights) == n, all(is.finite(weights)), all(weights > 0))
  }

  if (is.null(block_id)) {
    if (is.null(separators)) {
      separators <- floor(seq(0L, n, length.out = n_blocks + 1L))
      separators[length(separators)] <- n
    }
    n_blocks <- length(separators) - 1L
    stopifnot(n_blocks >= 2L, separators[1] == 0L, separators[length(separators)] == n)
    block_id <- integer(n)
    for (b in seq_len(n_blocks)) {
      block_id[(separators[b] + 1L):separators[b + 1L]] <- b
    }
  } else {
    block_id <- as.integer(as.factor(block_id))
    stopifnot(length(block_id) == n)
    n_blocks <- length(unique(block_id))
    stopifnot(n_blocks >= 2L)
    separators <- NULL
  }

  list(
    X = X, y = y, weights = weights, block_id = block_id,
    separators = separators, n = n, p = p, n_blocks = n_blocks
  )
}

.block_jackknife_fit_R <- function(X, y, weights, block_id, keep_delete) {
  n_blocks <- length(unique(block_id))
  p <- ncol(X)
  block_index <- split(seq_len(nrow(X)), block_id)
  XtY_blk <- matrix(0, n_blocks, p)
  XtX_blk <- array(0, dim = c(n_blocks, p, p))
  XtY_tot <- numeric(p)
  XtX_tot <- matrix(0, p, p)

  for (b in seq_len(n_blocks)) {
    idx <- block_index[[b]]
    Xb <- X[idx, , drop = FALSE]
    yb <- y[idx]
    wb <- sqrt(weights[idx])
    Xbw <- Xb * wb
    ybw <- yb * wb

    xty <- drop(crossprod(Xbw, ybw))
    xtx <- crossprod(Xbw)
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
    if (keep_delete) {
      delete_vals[b, ] <- beta_b
    }
    pseudo[b, ] <- n_blocks * est - (n_blocks - 1L) * beta_b
  }

  pseudo_centered <- sweep(pseudo, 2L, colMeans(pseudo), "-")
  cov_jk <- crossprod(pseudo_centered) / ((n_blocks - 1L) * n_blocks)
  se_jk <- sqrt(diag(cov_jk))

  list(
    coefficients = est,
    se = se_jk,
    cov = cov_jk,
    delete_values = delete_vals,
    n_blocks = n_blocks
  )
}

block_jackknife_fit <- function(X, y, n_blocks = 200L, separators = NULL, keep_delete = TRUE,
                                weights = NULL, block_id = NULL, backend = c("Rcpp", "R")) {
  backend <- match.arg(backend)
  prep <- .block_jackknife_prepare(
    X = X, y = y, n_blocks = n_blocks, separators = separators, weights = weights, block_id = block_id
  )

  use_rcpp <- identical(backend, "Rcpp")
  if (use_rcpp && !exists("bj_core_cpp", mode = "function")) {
    warning("Rcpp backend requested but compiled backend is unavailable; falling back to backend='R'.")
    use_rcpp <- FALSE
  }

  fit_core <- if (use_rcpp) {
    bj_core_cpp(prep$X, prep$y, prep$weights, prep$block_id, keep_delete)
  } else {
    .block_jackknife_fit_R(prep$X, prep$y, prep$weights, prep$block_id, keep_delete)
  }

  fit_core$separators <- prep$separators
  fit_core$block_id <- prep$block_id
  fit_core$weights <- prep$weights
  fit_core$n <- prep$n
  fit_core$p <- prep$p
  fit_core$n_blocks <- prep$n_blocks
  fit_core$backend <- if (use_rcpp) "Rcpp" else "R"
  fit_core
}

bjlm <- function(formula, data = NULL, n_blocks = 200L, separators = NULL,
                 na.action = na.omit, keep_delete = FALSE, weights = NULL, cluster = NULL,
                 backend = c("Rcpp", "R")) {
  cl <- match.call()
  backend <- match.arg(backend)

  if (inherits(formula, "formula")) {
    mf <- stats::model.frame(formula = formula, data = data, na.action = na.action)
    mt <- stats::terms(mf)
    y <- stats::model.response(mf)
    X <- stats::model.matrix(mt, mf)
    if (is.null(weights)) {
      w <- rep(1, nrow(X))
    } else {
      if (is.character(weights) && length(weights) == 1L && !is.null(data)) {
        w_all <- as.numeric(data[[weights]])
      } else {
        w_all <- as.numeric(weights)
      }
      if (!is.null(data)) {
        data_rows <- rownames(data)
        if (is.null(data_rows)) {
          data_rows <- as.character(seq_len(nrow(data)))
        }
        idx <- match(rownames(mf), data_rows)
      } else {
        idx <- as.integer(rownames(mf))
      }
      if (anyNA(idx)) {
        stop("Could not align weights with model frame rows.")
      }
      w <- w_all[idx]
    }
    if (is.null(cluster)) {
      cluster_vec <- NULL
    } else {
      if (is.character(cluster) && length(cluster) == 1L && !is.null(data)) {
        cluster_all <- data[[cluster]]
      } else {
        cluster_all <- cluster
      }
      if (!is.null(data)) {
        data_rows <- rownames(data)
        if (is.null(data_rows)) {
          data_rows <- as.character(seq_len(nrow(data)))
        }
        idx <- match(rownames(mf), data_rows)
      } else {
        idx <- as.integer(rownames(mf))
      }
      if (anyNA(idx)) {
        stop("Could not align cluster variable with model frame rows.")
      }
      cluster_vec <- cluster_all[idx]
    }
    coef_names <- colnames(X)
  } else {
    # Matrix interface through this wrapper: bjlm(X, y, ...)
    X <- as.matrix(formula)
    y <- as.numeric(data)
    if (is.null(weights)) {
      w <- rep(1, nrow(X))
    } else {
      w <- as.numeric(weights)
    }
    if (is.null(cluster)) {
      cluster_vec <- NULL
    } else {
      cluster_vec <- cluster
    }
    if (is.null(colnames(X))) {
      coef_names <- paste0("V", seq_len(ncol(X)))
      colnames(X) <- coef_names
    } else {
      coef_names <- colnames(X)
    }
    mt <- NULL
    mf <- NULL
  }

  if (!is.null(cluster_vec)) {
    cluster_vec <- as.vector(cluster_vec)
    stopifnot(length(cluster_vec) == nrow(X))
    cluster_fac <- as.factor(cluster_vec)
    cluster_sizes <- as.integer(table(cluster_fac))
    n_clusters <- length(cluster_sizes)
    n_blocks_requested <- n_blocks

    if (n_clusters < n_blocks) {
      warning("Number of clusters is less than n_blocks; reducing number of blocks to n_clusters.")
      n_blocks <- n_clusters
    }
    if (max(cluster_sizes) > (nrow(X) / n_blocks_requested)) {
      warning("Largest cluster is larger than n / n_blocks.")
    }

    # Assign whole clusters to jackknife blocks without splitting clusters.
    cluster_levels <- levels(cluster_fac)
    target <- nrow(X) / n_blocks
    block_for_cluster <- integer(length(cluster_levels))
    current_block <- 1L
    current_size <- 0L
    for (i in seq_along(cluster_levels)) {
      remaining_clusters <- length(cluster_levels) - i + 1L
      remaining_blocks <- n_blocks - current_block + 1L
      if (remaining_clusters == remaining_blocks) {
        block_for_cluster[i:length(cluster_levels)] <- current_block:n_blocks
        break
      }
      sz <- cluster_sizes[i]
      if (current_block < n_blocks && current_size > 0L && (current_size + sz) > target) {
        current_block <- current_block + 1L
        current_size <- 0L
      }
      block_for_cluster[i] <- current_block
      current_size <- current_size + sz
    }
    names(block_for_cluster) <- cluster_levels
    block_id <- block_for_cluster[as.character(cluster_fac)]
  } else {
    block_id <- NULL
  }

  fit <- block_jackknife_fit(
    X = X,
    y = y,
    n_blocks = n_blocks,
    separators = separators,
    keep_delete = keep_delete,
    weights = w,
    block_id = block_id,
    backend = backend
  )

  names(fit$coefficients) <- coef_names
  names(fit$se) <- coef_names
  dimnames(fit$cov) <- list(coef_names, coef_names)

  fit$call <- cl
  fit$terms <- mt
  fit$model <- mf
  fit$weights <- w
  fit$cluster <- cluster_vec
  fit$fitted.values <- drop(X %*% fit$coefficients)
  fit$residuals <- y - fit$fitted.values
  fit$df.residual <- nrow(X) - ncol(X)
  fit$rank <- ncol(X)
  fit$coeff_names <- coef_names
  class(fit) <- "bjlm"
  fit
}

summary.bjlm <- function(object, ...) {
  z <- object$coefficients / object$se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  tab <- cbind(
    Estimate = object$coefficients,
    `Jackknife SE` = object$se,
    `z value` = z,
    `Pr(>|z|)` = p
  )

  out <- list(
    call = object$call,
    coefficients = tab,
    sigma = sqrt(sum(object$weights * object$residuals^2) / max(1, object$df.residual)),
    df = c(object$rank, object$df.residual),
    n_blocks = object$n_blocks,
    backend = object$backend
  )
  class(out) <- "summary.bjlm"
  out
}

print.bjlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(cbind(Estimate = x$coefficients, `Jackknife SE` = x$se), digits = digits)
  invisible(x)
}

print.summary.bjlm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = getOption("show.signif.stars"))
  cat(
    sprintf(
      "\nResidual standard error: %.4f on %d degrees of freedom\n",
      x$sigma,
      x$df[2]
    )
  )
  cat(sprintf("Jackknife blocks: %d\n", x$n_blocks))
  cat(sprintf("Backend: %s\n", x$backend))
  invisible(x)
}

# Backward-compatible alias.
block_jackknife_lm <- bjlm

summary.block_jackknife_lm <- summary.bjlm
print.block_jackknife_lm <- print.bjlm
print.summary.block_jackknife_lm <- print.summary.bjlm
