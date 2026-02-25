# blockjack

Fast block jackknife standard errors for linear regression.

## Install locally

```r
# from the package root
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
