// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

// Pearson correlation of two equal-length vectors (no missing-value
// handling here — callers only pass already-filtered complete cases).
static double corr3(const std::vector<double>& a, const std::vector<double>& b, int n) {
  double ma = 0.0, mb = 0.0;
  for (int i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  double sab = 0.0, saa = 0.0, sbb = 0.0;
  for (int i = 0; i < n; i++) {
    double da = a[i] - ma, db = b[i] - mb;
    sab += da * db; saa += da * da; sbb += db * db;
  }
  double denom = std::sqrt(saa * sbb);
  return (denom > 1e-15) ? (sab / denom) : NA_REAL;
}

// ═════════════════════════════════════════════════════════════════════════
//  Classic partial Mantel test (Yule 1907 partial-correlation formula),
//  native C++ engine — same statistic and permutation scheme as this app's
//  R implementation (which itself matches vegan::mantel.partial() and
//  ade4::mantel.rtest()'s own conventions): permute the row/column labels
//  of X (the variable of interest) only, keep Y (response) and Z (control)
//  fixed, recompute part.cor(rxy, rxz, ryz) each time. Spearman-style
//  testing is obtained by rank-transforming X, Y, Z in R before calling
//  this function — it always just computes Pearson correlations on
//  whatever values it receives.
// ═════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
List classic_partial_mantel_cpp(NumericMatrix Xmat, NumericMatrix Ymat, NumericMatrix Zmat, int nperm) {
  RNGScope scope;  // synchronise with R's set.seed()

  int n = Xmat.nrow();
  if (Xmat.ncol() != n || Ymat.nrow() != n || Ymat.ncol() != n || Zmat.nrow() != n || Zmat.ncol() != n)
    stop("Xmat, Ymat and Zmat must be square matrices of the same size.");

  std::vector<int> pi_, pj_;
  std::vector<double> xv, yv, zv;
  for (int i = 1; i < n; i++) {
    for (int j = 0; j < i; j++) {
      double x = Xmat(i, j), y = Ymat(i, j), z = Zmat(i, j);
      if (R_finite(x) && R_finite(y) && R_finite(z)) {
        pi_.push_back(i); pj_.push_back(j);
        xv.push_back(x);  yv.push_back(y);  zv.push_back(z);
      }
    }
  }
  int npairs = (int) xv.size();
  if (npairs < 4) stop("Not enough complete (X, Y, Z) triples to run the test.");

  double rxy = corr3(xv, yv, npairs);
  double rxz = corr3(xv, zv, npairs);
  double ryz = corr3(yv, zv, npairs);
  double statistic = (rxy - rxz * ryz) / std::sqrt((1.0 - rxz * rxz) * (1.0 - ryz * ryz));
  if (!R_finite(statistic))
    stop("Could not compute the partial correlation (check for collinearity between X and Z).");

  std::vector<int> p(n);
  for (int i = 0; i < n; i++) p[i] = i;

  NumericVector permStats(nperm);
  std::vector<double> xp(npairs);
  const double EPS = 1e-9;
  int bpos = 0, bneg = 0, mValid = 0;

  for (int rep = 0; rep < nperm; rep++) {
    for (int t = n - 1; t > 0; t--) {
      int k = (int) std::floor(unif_rand() * (t + 1));
      std::swap(p[k], p[t]);
    }
    for (int k = 0; k < npairs; k++) {
      int a = p[pi_[k]], b = p[pj_[k]];
      int lo = std::min(a, b), hi = std::max(a, b);
      xp[k] = Xmat(hi, lo);
    }
    double rxy_p = corr3(xp, yv, npairs);
    double rxz_p = corr3(xp, zv, npairs);
    double s = (rxy_p - rxz_p * ryz) / std::sqrt((1.0 - rxz_p * rxz_p) * (1.0 - ryz * ryz));
    permStats[rep] = s;
    if (R_finite(s)) {
      mValid++;
      if (s >= statistic - EPS) bpos++;
      if (s <= statistic + EPS) bneg++;
    }
  }

  double pPos = (mValid > 0) ? (double) (bpos + 1) / (double) (mValid + 1) : NA_REAL;
  double pNeg = (mValid > 0) ? (double) (bneg + 1) / (double) (mValid + 1) : NA_REAL;

  return List::create(
    _["statistic"]  = statistic,
    _["rxy"]        = rxy,
    _["rxz"]        = rxz,
    _["ryz"]        = ryz,
    _["p_pos"]      = pPos,
    _["p_neg"]      = pNeg,
    _["perm_stats"] = permStats,
    _["n_pairs"]    = npairs,
    _["n_pops"]     = n
  );
}
