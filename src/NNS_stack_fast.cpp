#include <Rcpp.h>
#include <unordered_map>
#include <cstring>
using namespace Rcpp;

// [[Rcpp::export]]
List NNS_xstar_path_cpp(const NumericMatrix& train_design,
                        const NumericMatrix& test_design,
                        const NumericVector& coefficients,
                        const IntegerVector& column_order,
                        int nthreads) {
  int ntr=train_design.nrow(), nte=test_design.nrow(), p=train_design.ncol();
  if(test_design.ncol()!=p || coefficients.size()!=p) stop("incompatible dimensions");
  NumericMatrix train(ntr,p), test(nte,p); NumericVector denom(p); IntegerVector rep(p);
  std::vector<double> num_tr(ntr,0.0), num_te(nte,0.0); double den=0.0; bool any_nonzero=false;
  for(int c=0;c<p;++c) if(coefficients[c]!=0.0) any_nonzero=true;
  for(int m=0;m<p;++m){
    int j=column_order[m]-1; double coef=coefficients[j];
    if(!any_nonzero) { coef=1.0; den += 1.0; }
    else if(coef!=0.0) den += 1.0;
    for(int i=0;i<ntr;++i) num_tr[i] += train_design(i,j)*coef;
    for(int i=0;i<nte;++i) num_te[i] += test_design(i,j)*coef;
    denom[m]=den;
    for(int i=0;i<ntr;++i) train(i,m)=num_tr[i]/den;
    for(int i=0;i<nte;++i) test(i,m)=num_te[i]/den;
  }
  for(int m=0;m<p;++m){ rep[m]=m+1; for(int r=0;r<m;++r){ bool same=true; for(int i=0;i<ntr && same;++i) same=(train(i,m)==train(i,r)); for(int i=0;i<nte && same;++i) same=(test(i,m)==test(i,r)); if(same){rep[m]=r+1; break;} } }
  return List::create(_["train"]=train, _["test"]=test, _["denominator"]=denom, _["column_order"]=column_order, _["representative"]=rep);
}
