// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// Ordinary least squares via Gauss-Jordan elimination (with partial
// pivoting) on the normal equations (X'X) beta = X'y. X is n x p
// (p = number of predictors including the intercept column). Returns
// false if the system is (numerically) singular — i.e. collinear
// predictors — mirroring R's own lm() failing/producing NA coefficients
// in that situation.
static bool ols_fit(const std::vector<std::vector<double> >& X, const std::vector<double>& y,
                     int n, int p, std::vector<double>& beta, double& r2) {
  std::vector<std::vector<double> > A(p, std::vector<double>(p, 0.0));
  std::vector<double> b(p, 0.0);
  for (int i = 0; i < p; i++) {
    for (int j = i; j < p; j++) {
      double s = 0.0;
      for (int r = 0; r < n; r++) s += X[r][i] * X[r][j];
      A[i][j] = s; A[j][i] = s;
    }
    double sb = 0.0;
    for (int r = 0; r < n; r++) sb += X[r][i] * y[r];
    b[i] = sb;
  }
  for (int col = 0; col < p; col++) {
    int piv = col; double best = std::fabs(A[col][col]);
    for (int r = col + 1; r < p; r++) if (std::fabs(A[r][col]) > best) { best = std::fabs(A[r][col]); piv = r; }
    if (best < 1e-10) return false;
    if (piv != col) { std::swap(A[piv], A[col]); std::swap(b[piv], b[col]); }
    double d = A[col][col];
    for (int j = 0; j < p; j++) A[col][j] /= d;
    b[col] /= d;
    for (int r = 0; r < p; r++) {
      if (r == col) continue;
      double f = A[r][col];
      if (f == 0.0) continue;
      for (int j = 0; j < p; j++) A[r][j] -= f * A[col][j];
      b[r] -= f * b[col];
    }
  }
  beta = b;
  double ybar = 0.0; for (int r = 0; r < n; r++) ybar += y[r]; ybar /= n;
  double sst = 0.0, sse = 0.0;
  for (int r = 0; r < n; r++) {
    double yhat = 0.0;
    for (int j = 0; j < p; j++) yhat += X[r][j] * beta[j];
    double resid = y[r] - yhat;
    sse += resid * resid;
    double dd = y[r] - ybar;
    sst += dd * dd;
  }
  r2 = (sst > 1e-12) ? (1.0 - sse / sst) : NA_REAL;
  return true;
}

static void standardize_vec(std::vector<double>& v) {
  int n = (int) v.size();
  double m = 0.0; for (double x : v) m += x; m /= n;
  double s = 0.0; for (double x : v) s += (x - m) * (x - m);
  s = std::sqrt(s / (n - 1));
  if (s < 1e-12) { for (auto& x : v) x -= m; return; }
  for (auto& x : v) x = (x - m) / s;
}

// ═════════════════════════════════════════════════════════════════════════
//  MRM (multiple regression on distance matrices; Legendre, Lapointe &
//  Casgrain 1994) — native C++ engine, same statistic and permutation
//  scheme as this app's R implementation: joint row/column relabelling of
//  the response (Y) matrix, refitting Y ~ X1 + ... + Xk (OLS, with an
//  intercept) on each replicate. Xmats is a list of k square matrices
//  (predictors); standardize z-scores Y and every predictor over the
//  current complete-case subset before fitting, exactly mirroring the R
//  version's .pm_std() behaviour.
// ═════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List mrm_cpp(NumericMatrix Ymat, List Xmats, int nperm, bool standardize) {
  RNGScope scope;

  int n = Ymat.nrow();
  int k = Xmats.size();
  if (k < 1 || k > 10) stop("Between 1 and 10 predictor matrices are supported.");
  std::vector<NumericMatrix> Xm(k);
  for (int m = 0; m < k; m++) {
    NumericMatrix Xi = as<NumericMatrix>(Xmats[m]);
    if (Xi.nrow() != n || Xi.ncol() != n) stop("All matrices must be square and the same size as Y.");
    Xm[m] = Xi;
  }

  std::vector<int> pi_, pj_;
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      double y = Ymat(i, j);
      bool ok = R_finite(y);
      for (int m = 0; m < k && ok; m++) if (!R_finite(Xm[m](i, j))) ok = false;
      if (ok) { pi_.push_back(i); pj_.push_back(j); }
    }
  }
  int npairs = (int) pi_.size();
  if (npairs < k + 3) stop("Not enough complete dyads across Y and all predictors for this many predictors.");

  std::vector<std::vector<double> > Xraw(npairs, std::vector<double>(k));
  for (int r = 0; r < npairs; r++)
    for (int m = 0; m < k; m++) Xraw[r][m] = Xm[m](pi_[r], pj_[r]);
  std::vector<double> Yraw(npairs);
  for (int r = 0; r < npairs; r++) Yraw[r] = Ymat(pi_[r], pj_[r]);

  auto fitOnce = [&](const std::vector<double>& yv, std::vector<double>& coefOut, double& r2Out) -> bool {
    std::vector<double> yStd = yv;
    std::vector<std::vector<double> > Xcols(k, std::vector<double>(npairs));
    for (int m = 0; m < k; m++) for (int r = 0; r < npairs; r++) Xcols[m][r] = Xraw[r][m];
    if (standardize) {
      standardize_vec(yStd);
      for (int m = 0; m < k; m++) standardize_vec(Xcols[m]);
    }
    std::vector<std::vector<double> > Xdesign(npairs, std::vector<double>(k + 1));
    for (int r = 0; r < npairs; r++) {
      Xdesign[r][0] = 1.0;
      for (int m = 0; m < k; m++) Xdesign[r][m + 1] = Xcols[m][r];
    }
    return ols_fit(Xdesign, yStd, npairs, k + 1, coefOut, r2Out);
  };

  std::vector<double> obsCoef; double obsR2 = NA_REAL;
  bool ok0 = fitOnce(Yraw, obsCoef, obsR2);
  if (!ok0) stop("Could not fit the model (check for collinear predictors).");

  std::vector<int> p(n);
  for (int i = 0; i < n; i++) p[i] = i;

  NumericMatrix permCoef(nperm, k);
  NumericVector permR2(nperm);
  std::vector<double> yPerm(npairs);

  for (int rep = 0; rep < nperm; rep++) {
    for (int t = n - 1; t > 0; t--) {
      int kk = (int) std::floor(unif_rand() * (t + 1));
      std::swap(p[kk], p[t]);
    }
    for (int r = 0; r < npairs; r++) {
      int a = p[pi_[r]], b = p[pj_[r]];
      int lo = std::min(a, b), hi = std::max(a, b);
      yPerm[r] = Ymat(hi, lo);
    }
    std::vector<double> coefv; double r2v = NA_REAL;
    bool okk = fitOnce(yPerm, coefv, r2v);
    if (okk) {
      for (int m = 0; m < k; m++) permCoef(rep, m) = coefv[m + 1];
      permR2[rep] = r2v;
    } else {
      for (int m = 0; m < k; m++) permCoef(rep, m) = NA_REAL;
      permR2[rep] = NA_REAL;
    }
  }

  NumericVector obsCoefOut(k);
  for (int m = 0; m < k; m++) obsCoefOut[m] = obsCoef[m + 1];

  return List::create(
    _["obs_coef"]  = obsCoefOut,
    _["obs_r2"]    = obsR2,
    _["perm_coef"] = permCoef,
    _["perm_r2"]   = permR2,
    _["n_pairs"]   = npairs,
    _["n_pops"]    = n
  );
}
