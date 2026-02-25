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
x1 <- rnorm(n)
x2 <- rnorm(n)
y <- 0.5 + 0.7 * x1 - 0.2 * x2 + rnorm(n)
d <- data.frame(y = y, x1 = x1, x2 = x2)

fit <- bjlm(y ~ x1 + x2, data = d, n_blocks = 200)
summary(fit)
```

Example output:

```text
Call:
bjlm(formula = y ~ x1 + x2, data = d, n_blocks = 200)

Coefficients:
             Estimate Jackknife SE z value Pr(>|z|)    
(Intercept)  0.493175     0.006641   74.26   <2e-16 ***
x1           0.686438     0.007263   94.51   <2e-16 ***
x2          -0.196584     0.007320  -26.86   <2e-16 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 1.0071 on 19997 degrees of freedom
Jackknife blocks: 200
```
