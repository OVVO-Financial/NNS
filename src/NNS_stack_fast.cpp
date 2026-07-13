#include <Rcpp.h>
#include <RcppParallel.h>
#include <vector>
using namespace Rcpp;
using namespace RcppParallel;

struct FillProjectionWorker : public Worker {
  RMatrix<double> design;
  RMatrix<double> out;
  const std::vector<int>& order;
  const std::vector<double>& coef;
  const std::vector<double>& denom;
  const bool all_zero;

  FillProjectionWorker(const NumericMatrix& design_, NumericMatrix& out_,
                       const std::vector<int>& order_,
                       const std::vector<double>& coef_,
                       const std::vector<double>& denom_, bool all_zero_)
      : design(design_), out(out_), order(order_), coef(coef_),
        denom(denom_), all_zero(all_zero_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const int p = out.ncol();
    for (std::size_t i = begin; i < end; ++i) {
      double numerator = 0.0;
      for (int m = 0; m < p; ++m) {
        const int j = order[m];
        const double c = all_zero ? 1.0 : coef[j];
        numerator += design(i, j) * c;
        out(i, m) = numerator / denom[m];
      }
    }
  }
};

// [[Rcpp::export]]
List NNS_xstar_path_cpp(const NumericMatrix& train_design,
                        const NumericMatrix& test_design,
                        const NumericVector& coefficients,
                        const IntegerVector& column_order,
                        int nthreads) {
  const int ntr = train_design.nrow(), nte = test_design.nrow(), p = train_design.ncol();
  if (test_design.ncol() != p || coefficients.size() != p || column_order.size() != p) {
    stop("incompatible dimensions");
  }
  NumericMatrix train(ntr, p), test(nte, p);
  NumericVector denominator(p);
  IntegerVector representative(p);
  std::vector<int> order(p);
  std::vector<double> coef(p), denom(p);
  bool all_zero = true;
  for (int j = 0; j < p; ++j) {
    order[j] = column_order[j] - 1;
    coef[j] = coefficients[j];
    if (coef[j] != 0.0) all_zero = false;
  }
  double den = 0.0;
  for (int m = 0; m < p; ++m) {
    const int j = order[m];
    if (all_zero || coef[j] != 0.0) den += 1.0;
    denom[m] = den;
    denominator[m] = den;
  }
  const int threads = std::max(1, nthreads);
  FillProjectionWorker train_worker(train_design, train, order, coef, denom, all_zero);
  FillProjectionWorker test_worker(test_design, test, order, coef, denom, all_zero);
  if (threads == 1 || ntr < 2) train_worker(0, ntr); else parallelFor(0, ntr, train_worker, 1, threads);
  if (threads == 1 || nte < 2) test_worker(0, nte); else parallelFor(0, nte, test_worker, 1, threads);

  for (int m = 0; m < p; ++m) {
    representative[m] = m + 1;
    for (int r = 0; r < m; ++r) {
      bool same = true;
      for (int i = 0; i < ntr && same; ++i) same = (train(i, m) == train(i, r));
      for (int i = 0; i < nte && same; ++i) same = (test(i, m) == test(i, r));
      if (same) { representative[m] = r + 1; break; }
    }
  }
  return List::create(_["train"] = train, _["test"] = test,
                      _["denominator"] = denominator,
                      _["column_order"] = column_order,
                      _["representative"] = representative);
}
