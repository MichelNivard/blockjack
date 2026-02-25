# blockjack

Fast block jackknife standard errors for linear regression.

## What It Does

`blockjack` computes linear-model coefficients and block jackknife standard
errors efficiently by:

1. Splitting rows into blocks.
2. Precomputing block-wise `X'X` and `X'y`.
3. Reusing those crossproducts for all leave-one-block-out fits.

This avoids refitting the model from scratch for each block and is much faster
than naive jackknife loops when `n_blocks` is large.

## Method Source

The computational trick follows the LD Score regression jackknife approach:

Bulik-Sullivan, B., Loh, PR., Finucane, H. et al. LD Score regression
distinguishes confounding from polygenicity in genome-wide association studies.
Nat Genet 47, 291-295 (2015). https://doi.org/10.1038/ng.3211

## Install locally

```r
# install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("MichelNivard/blockjack")
```

```r
# or install from a local checkout (from package root)
install.packages(".", repos = NULL, type = "source")
```

## Usage

```r
library(blockjack)

set.seed(1)
n <- 20000
X <- cbind(1, matrix(rnorm(n * 2), n, 2))
beta <- c(0.5, -0.2, 0.1)
y <- drop(X %*% beta + rnorm(n))

fit <- block_jackknife_fit(X, y, n_blocks = 200)
fit$coefficients
fit$se
```

Formula interface:

```r
d <- data.frame(y = y, x1 = X[, 2], x2 = X[, 3])
fit2 <- bjlm(y ~ x1 + x2, data = d, n_blocks = 200)
summary(fit2)
```
