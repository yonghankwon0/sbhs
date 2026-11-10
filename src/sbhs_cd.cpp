// -----------------------------------------------------------------------------
// One sweep of coordinate descent for the SBHS objective
//
//   (1/2n) ||y - X b||^2  +  sum_j MCP(|b_j|; lambda1, gamma)  +  (lambda2/2) b' H b
//
// The columns of X are assumed to be scaled so that ||x_j||^2 = n. With that
// scaling the one-dimensional problem in b_j has a closed-form solution (the
// "MCP thresholding operator", Section 4.1 of the paper), which is what
// mcp_update() computes.
//
// This file is compiled from R with Rcpp::sourceCpp(); see R/sbhs.R.
// -----------------------------------------------------------------------------
#include <Rcpp.h>
using namespace Rcpp;

// Soft thresholding: S(z, lambda) = sign(z) * max(|z| - lambda, 0)
// [[Rcpp::export]]
double soft_thresh(double z, double lambda) {
  if (z > lambda) return z - lambda;
  else if (z < -lambda) return z + lambda;
  else return 0.0;
}

// Minimiser of  0.5 * w * b^2 - z * b + MCP(|b|; lambda1, gamma)  over b.
//   z : the "partial residual" for coordinate j (data fit plus graph term)
//   w : the curvature 1 + lambda2 * H_jj  (always >= 1 > 1/gamma, so the
//       one-dimensional problem is strictly convex)
// [[Rcpp::export]]
double mcp_update(double z, double w, double lambda1, double gamma) {
  double thresh = gamma * lambda1 * w;
  if (std::abs(z) <= thresh) {
    // inside the concave part of the MCP: rescaled soft thresholding
    double denom = w - 1.0 / gamma;
    if (denom > 0) return soft_thresh(z, lambda1) / denom;
    else return 0.0;
  } else {
    // beyond gamma * lambda1 the MCP is flat: plain ridge-type update
    return z / w;
  }
}

// One full pass over the p coordinates.
//   X_scaled : n x p design with ||x_j||^2 = n
//   H        : p x p penalty matrix (Bethe-Hessian, Laplacian, ...)
//   beta     : current coefficients (on the scaled-X scale)
//   residual : current residual y - X_scaled %*% beta
//   H_diag   : diagonal of H
// Returns the updated coefficients, the updated residual and the largest
// absolute change of a coefficient in this pass (used as stopping rule).
// [[Rcpp::export]]
List cd_sweep(NumericMatrix X_scaled, NumericMatrix H, NumericVector beta,
              NumericVector residual, NumericVector H_diag,
              double lambda1, double lambda2, double gamma, int n, int p) {

  NumericVector beta_new = clone(beta);
  NumericVector resid_new = clone(residual);
  double max_change = 0.0;

  for (int j = 0; j < p; j++) {
    // data-fit part: x_j' r / n + b_j
    double dot_xj_r = 0.0;
    for (int i = 0; i < n; i++) dot_xj_r += X_scaled(i, j) * resid_new[i];
    double z_j = dot_xj_r / n + beta_new[j];

    // graph part: - lambda2 * sum_{k != j} H_jk b_k
    // (for the Bethe-Hessian H_jk = -r A_jk, so neighbours pull b_j towards
    //  the sign-coherent value)
    double h_contrib = 0.0;
    for (int k = 0; k < p; k++) if (k != j) h_contrib += H(j, k) * beta_new[k];
    z_j -= lambda2 * h_contrib;

    double w_j = 1.0 + lambda2 * H_diag[j];
    double beta_new_j = mcp_update(z_j, w_j, lambda1, gamma);

    if (beta_new_j != beta_new[j]) {
      double diff = beta_new_j - beta_new[j];
      for (int i = 0; i < n; i++) resid_new[i] -= diff * X_scaled(i, j);
      double change = std::abs(diff);
      if (change > max_change) max_change = change;
      beta_new[j] = beta_new_j;
    }
  }

  return List::create(Named("beta") = beta_new,
                      Named("residual") = resid_new,
                      Named("max_change") = max_change);
}
