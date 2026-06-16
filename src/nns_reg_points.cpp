// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "central_tendencies.h"

using namespace Rcpp;

// ===========================================================================
// NNS.reg multivariate regression-point builder (native fast path)
//
// Reproduces, bit-for-bit, the regression.points pipeline in R/Regression.R
// (the block that NNS.reg returns for multivariate.call = TRUE) for the gated
// configuration: type = NULL, noise.reduction = "off", dependence < 1,
// dep.reduced.order != "max", smooth = FALSE.  All other configurations fall
// back to the pure-R path.
//
// This exists because NNS.ARMA / NNS.stack / NNS.boost call NNS.reg hundreds of
// times on small inputs, where per-call data.table overhead dominates.
// ===========================================================================

// --- mean via long double accumulation.  Matches both base R mean() and
//     Rcpp sugar mean() (verified identical), the latter used by fast_lm.
static double r_mean_v(const std::vector<double>& v) {
  long double s = 0.0L; for (double d : v) s += d;
  return (double) (s / (long double) v.size());
}

// --- base R mean(c(a, b)) : two values via long double, matching R exactly
static double r_mean2(double a, double b) {
  long double s = (long double) a + (long double) b;
  return (double) (s / 2.0L);
}

// --- gravity() : filter to finite then call shared core (continuous)
static double grav(const std::vector<double>& v) {
  std::vector<double> f; f.reserve(v.size());
  for (double d : v) if (R_finite(d)) f.push_back(d);
  return gravity_value(f, false);
}

// --- OLS fit replicating src/fast_lm.cpp exactly (intercept a, slope b)
static void ols_fit(const std::vector<double>& xs, const std::vector<double>& ys,
                    double& a, double& b) {
  const double mx = r_mean_v(xs);
  const double my = r_mean_v(ys);
  double vx = 0.0, cv = 0.0;
  for (size_t i = 0; i < xs.size(); ++i) {
    const double dx = xs[i] - mx, dy = ys[i] - my;
    vx += dx * dx; cv += dx * dy;
  }
  if (vx == 0.0) { a = my; b = 0.0; }
  else { b = cv / vx; a = my - b * mx; }
}

// --- number of distinct values (exact equality, matching base R unique)
static int count_unique(const std::vector<double>& v) {
  std::vector<double> t = v;
  std::sort(t.begin(), t.end());
  t.erase(std::unique(t.begin(), t.end()), t.end());
  return (int) t.size();
}

// --- unique() preserving first-occurrence order (matches base R unique)
static std::vector<double> unique_preserve(const std::vector<double>& v) {
  std::vector<double> out;
  for (double d : v) {
    bool seen = false;
    for (double e : out) if (e == d) { seen = true; break; }
    if (!seen) out.push_back(d);
  }
  return out;
}

// --- consolidate: group by exact x (ascending), y = gravity() of group's y.
//     Equivalent to setkey(rp,x); rp[, y := gravity(y), by="x"]; unique(rp).
static void consolidate(const std::vector<double>& xs, const std::vector<double>& ys,
                        std::vector<double>& ox, std::vector<double>& oy) {
  ox.clear(); oy.clear();
  const size_t n = xs.size();
  if (n == 0) return;
  std::vector<size_t> idx(n);
  for (size_t i = 0; i < n; ++i) idx[i] = i;
  std::stable_sort(idx.begin(), idx.end(),
                   [&](size_t a, size_t b) { return xs[a] < xs[b]; });
  size_t i = 0;
  while (i < n) {
    const double cx = xs[idx[i]];
    std::vector<double> grp;
    size_t j = i;
    while (j < n && xs[idx[j]] == cx) {
      if (R_finite(ys[idx[j]])) grp.push_back(ys[idx[j]]);
      ++j;
    }
    ox.push_back(cx);
    oy.push_back(gravity_value(grp, false));
    i = j;
  }
}

// [[Rcpp::export]]
DataFrame NNS_reg_points_cpp(NumericVector x_, NumericVector y_,
                             NumericVector rpx_, NumericVector rpy_,
                             double dependence, double stn) {
  const std::vector<double> x(x_.begin(), x_.end());
  const std::vector<double> y(y_.begin(), y_.end());

  // min/max of original x, y
  double minx = R_PosInf, maxx = R_NegInf, miny = R_PosInf, maxy = R_NegInf;
  for (double v : x) { if (v < minx) minx = v; if (v > maxx) maxx = v; }
  for (double v : y) { if (v < miny) miny = v; if (v > maxy) maxy = v; }

  // ---- Step A: clamp rp x into [minx, maxx], consolidate ----
  std::vector<double> cx(rpx_.size()), cy(rpy_.begin(), rpy_.end());
  for (int i = 0; i < rpx_.size(); ++i)
    cx[i] = std::min(maxx, std::max(rpx_[i], minx));     // pmin(maxx, pmax(rp, minx))
  std::vector<double> rx, ry;
  consolidate(cx, cy, rx, ry);
  const int N = (int) rx.size();

  // ---- Step B: central point (med.rps), type = NULL branch ----
  const double medpos = (N + 1) / 2.0;                   // median(1:N)
  const int m_lo = (int) std::floor(medpos);
  const int m_hi = (int) std::ceil(medpos);
  const double cxlo = rx[m_lo - 1];
  const double cxhi = rx[m_hi - 1];
  const bool two = (m_lo != m_hi);                       // length(unique(central_rows)) > 1

  double central_y;
  if (two) {
    std::vector<double> sub;
    for (size_t i = 0; i < x.size(); ++i)
      if (x[i] >= cxlo && x[i] <= cxhi) sub.push_back(y[i]);
    central_y = grav(sub);
  } else {
    central_y = ry[m_lo - 1];
  }
  std::vector<double> cc = { cxlo, cxhi };               // rp[central_rows,]$x (length 2)
  const double central_x = grav(cc);                     // gravity(central_x)

  // ---- Step C: append med.rps, complete.cases, consolidate ----
  std::vector<double> bx = rx, by = ry;
  bx.push_back(central_x); by.push_back(central_y);
  std::vector<double> fx, fy;
  // complete.cases() keeps Inf/-Inf and drops only NA/NaN (ISNAN covers both).
  for (size_t i = 0; i < bx.size(); ++i)
    if (!ISNAN(bx[i]) && !ISNAN(by[i])) { fx.push_back(bx[i]); fy.push_back(by[i]); }
  std::vector<double> cx2, cy2;
  consolidate(fx, fy, cx2, cy2);

  // ---- Step D: endpoints (dependence < 1, type = NULL) ----
  const double minr = *std::min_element(cx2.begin(), cx2.end());
  const double maxr = *std::max_element(cx2.begin(), cx2.end());
  const double mid_min = r_mean2(minx, minr);            // mean(c(min(x), min(rp$x)))
  const double mid_max = r_mean2(maxx, maxr);

  std::vector<double> y_min, y_midmin, x_midmin;
  std::vector<double> y_max, y_midmax, x_midmax;
  // na.omit() (like complete.cases) keeps Inf/-Inf and drops only NA/NaN.
  for (size_t i = 0; i < x.size(); ++i) {
    const double xi = x[i], yi = y[i];
    if (xi <= minr   && !ISNAN(yi)) y_min.push_back(yi);
    if (xi <= mid_min) { if (!ISNAN(yi)) y_midmin.push_back(yi); if (!ISNAN(xi)) x_midmin.push_back(xi); }
    if (xi >= maxr   && !ISNAN(yi)) y_max.push_back(yi);
    if (xi >= mid_max) { if (!ISNAN(yi)) y_midmax.push_back(yi); if (!ISNAN(xi)) x_midmax.push_back(xi); }
  }
  const int l_y_min  = (int) y_min.size();
  const int l_y_midmin = (int) y_midmin.size();
  const int l_y_max  = (int) y_max.size();
  const int l_y_midmax = (int) y_midmax.size();
  const int l_x_midmin_unique = count_unique(x_midmin);
  const int l_x_midmax_unique = count_unique(x_midmax);

  // y / x values where original x == min(x) / max(x)
  std::vector<double> y_at_minx, y_at_maxx;
  std::vector<double> xs_le_minr, ys_le_minr, xs_le_midmin, ys_le_midmin;
  std::vector<double> xs_ge_maxr, ys_ge_maxr, xs_ge_midmax, ys_ge_midmax;
  for (size_t i = 0; i < x.size(); ++i) {
    const double xi = x[i], yi = y[i];
    if (xi == minx) y_at_minx.push_back(yi);
    if (xi == maxx) y_at_maxx.push_back(yi);
    if (xi <= minr)   { xs_le_minr.push_back(xi);   ys_le_minr.push_back(yi); }
    if (xi <= mid_min){ xs_le_midmin.push_back(xi); ys_le_midmin.push_back(yi); }
    if (xi >= maxr)   { xs_ge_maxr.push_back(xi);   ys_ge_maxr.push_back(yi); }
    if (xi >= mid_max){ xs_ge_midmax.push_back(xi); ys_ge_midmax.push_back(yi); }
  }

  // --- min endpoint x0 ---
  std::vector<double> x0;
  if (l_x_midmin_unique > 1 && l_y_min > 5) {
    if (dependence < stn) {
      if (l_y_min > 1 && l_y_midmin > 1) {
        double a1, b1, a2, b2;
        ols_fit(xs_le_minr,   ys_le_minr,   a1, b1);
        ols_fit(xs_le_midmin, ys_le_midmin, a2, b2);
        const double f1 = a1 + b1 * (*std::min_element(xs_le_minr.begin(),   xs_le_minr.end()));
        const double f2 = a2 + b2 * (*std::min_element(xs_le_midmin.begin(), xs_le_midmin.end()));
        // R: sum(f1*l_y.min, f2*l_y.mid.min) / sum(l_y.min, l_y.mid.min)
        // each sum() returns a double, then the division is in double.
        const double num = (double)((long double)(f1 * l_y_min) + (long double)(f2 * l_y_midmin));
        x0.push_back(num / (double)(l_y_min + l_y_midmin));
      } else {
        x0 = y_min;
      }
    } else {
      x0 = unique_preserve(y_at_minx);
    }
  } else {
    x0.push_back(grav(y_at_minx));
  }
  const double min_rps_y = r_mean_v(x0);                 // mean(x0)

  // --- max endpoint x.max ---
  std::vector<double> xmaxv;
  if (l_x_midmax_unique > 1 && l_y_max > 5) {
    if (dependence < stn) {
      if (l_y_max > 1 && l_y_midmax > 1) {
        double a1, b1, a2, b2;
        ols_fit(xs_ge_maxr,   ys_ge_maxr,   a1, b1);
        ols_fit(xs_ge_midmax, ys_ge_midmax, a2, b2);
        const double f1 = a1 + b1 * (*std::max_element(xs_ge_maxr.begin(),   xs_ge_maxr.end()));
        const double f2 = a2 + b2 * (*std::max_element(xs_ge_midmax.begin(), xs_ge_midmax.end()));
        const double num = (double)((long double)(f1 * l_y_max) + (long double)(f2 * l_y_midmax));
        xmaxv.push_back(num / (double)(l_y_max + l_y_midmax));
      } else {
        xmaxv = y_max;
      }
    } else {
      xmaxv = unique_preserve(y_at_maxx);
    }
  } else {
    xmaxv.push_back(grav(y_at_maxx));
  }
  const double max_rps_y = r_mean_v(xmaxv);              // mean(x.max)

  // ---- Step E: append min/max/med rps, complete.cases, consolidate ----
  std::vector<double> ex = cx2, ey = cy2;
  ex.push_back(minx);      ey.push_back(min_rps_y);      // min.rps
  ex.push_back(maxx);      ey.push_back(max_rps_y);      // max.rps
  ex.push_back(central_x); ey.push_back(central_y);      // med.rps
  std::vector<double> gx, gy;
  // complete.cases() keeps Inf/-Inf and drops only NA/NaN.
  for (size_t i = 0; i < ex.size(); ++i)
    if (!ISNAN(ex[i]) && !ISNAN(ey[i])) { gx.push_back(ex[i]); gy.push_back(ey[i]); }
  std::vector<double> hx, hy;
  consolidate(gx, gy, hx, hy);

  // ---- Step F/G: single-row tripling, then clamp ----
  std::vector<double> ox, oy;
  if ((int) hx.size() == 1) {
    for (int k = 0; k < 3; ++k) { ox.push_back(hx[0]); oy.push_back(hy[0]); }
  } else {
    ox = hx; oy = hy;
  }
  for (size_t i = 0; i < ox.size(); ++i) {
    ox[i] = std::min(std::max(ox[i], minx), maxx);
    oy[i] = std::min(std::max(oy[i], miny), maxy);
  }

  return DataFrame::create(_["x"] = wrap(ox), _["y"] = wrap(oy));
}
