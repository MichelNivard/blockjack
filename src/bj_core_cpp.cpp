#include <Rcpp.h>
using namespace Rcpp;

static std::vector<double> solve_linear(std::vector<double> A, std::vector<double> b, int p) {
  for (int i = 0; i < p; ++i) {
    int piv = i;
    double best = std::abs(A[i * p + i]);
    for (int r = i + 1; r < p; ++r) {
      double v = std::abs(A[r * p + i]);
      if (v > best) {
        best = v;
        piv = r;
      }
    }
    if (best < 1e-12) stop("Singular linear system in jackknife solve.");

    if (piv != i) {
      for (int c = 0; c < p; ++c) {
        std::swap(A[i * p + c], A[piv * p + c]);
      }
      std::swap(b[i], b[piv]);
    }

    double diag = A[i * p + i];
    for (int c = i; c < p; ++c) A[i * p + c] /= diag;
    b[i] /= diag;

    for (int r = 0; r < p; ++r) {
      if (r == i) continue;
      double f = A[r * p + i];
      if (f == 0.0) continue;
      for (int c = i; c < p; ++c) A[r * p + c] -= f * A[i * p + c];
      b[r] -= f * b[i];
    }
  }
  return b;
}

// [[Rcpp::export]]
Rcpp::List bj_core_cpp(const NumericMatrix& X,
                       const NumericVector& y,
                       const NumericVector& weights,
                       const IntegerVector& block_id,
                       const bool keep_delete) {
  const int n = X.nrow();
  const int p = X.ncol();

  int B = 0;
  for (int i = 0; i < n; ++i) {
    if (block_id[i] > B) B = block_id[i];
  }
  if (B < 2) stop("Need at least 2 blocks.");

  std::vector<double> XtY_blk(B * p, 0.0);
  std::vector<double> XtX_blk(B * p * p, 0.0);

  for (int i = 0; i < n; ++i) {
    int b = block_id[i] - 1;
    double wi = weights[i];
    double yi = y[i];

    for (int j = 0; j < p; ++j) {
      double xij = X(i, j);
      double wxij = wi * xij;
      XtY_blk[b * p + j] += wxij * yi;
      for (int k = 0; k < p; ++k) {
        XtX_blk[(b * p + j) * p + k] += wxij * X(i, k);
      }
    }
  }

  std::vector<double> XtY_tot(p, 0.0);
  std::vector<double> XtX_tot(p * p, 0.0);
  for (int b = 0; b < B; ++b) {
    for (int j = 0; j < p; ++j) {
      XtY_tot[j] += XtY_blk[b * p + j];
      for (int k = 0; k < p; ++k) {
        XtX_tot[j * p + k] += XtX_blk[(b * p + j) * p + k];
      }
    }
  }

  std::vector<double> est = solve_linear(XtX_tot, XtY_tot, p);

  NumericMatrix delete_vals;
  if (keep_delete) delete_vals = NumericMatrix(B, p);
  NumericMatrix pseudo(B, p);

  for (int b = 0; b < B; ++b) {
    std::vector<double> XtY_del = XtY_tot;
    std::vector<double> XtX_del = XtX_tot;

    for (int j = 0; j < p; ++j) {
      XtY_del[j] -= XtY_blk[b * p + j];
      for (int k = 0; k < p; ++k) {
        XtX_del[j * p + k] -= XtX_blk[(b * p + j) * p + k];
      }
    }

    std::vector<double> beta_b = solve_linear(XtX_del, XtY_del, p);
    for (int j = 0; j < p; ++j) {
      if (keep_delete) delete_vals(b, j) = beta_b[j];
      pseudo(b, j) = (double)B * est[j] - ((double)B - 1.0) * beta_b[j];
    }
  }

  std::vector<double> mu(p, 0.0);
  for (int j = 0; j < p; ++j) {
    for (int b = 0; b < B; ++b) mu[j] += pseudo(b, j);
    mu[j] /= (double)B;
  }

  NumericMatrix cov_jk(p, p);
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      double s = 0.0;
      for (int b = 0; b < B; ++b) {
        s += (pseudo(b, j) - mu[j]) * (pseudo(b, k) - mu[k]);
      }
      cov_jk(j, k) = s / (((double)B - 1.0) * (double)B);
    }
  }

  NumericVector est_out(p), se_out(p);
  for (int j = 0; j < p; ++j) {
    est_out[j] = est[j];
    se_out[j] = std::sqrt(std::max(0.0, cov_jk(j, j)));
  }

  List out = List::create(
    Named("coefficients") = est_out,
    Named("se") = se_out,
    Named("cov") = cov_jk
  );
  if (keep_delete) out["delete_values"] = delete_vals;
  else out["delete_values"] = R_NilValue;
  return out;
}
