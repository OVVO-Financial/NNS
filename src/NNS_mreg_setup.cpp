#include <Rcpp.h>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <unordered_map>
using namespace Rcpp;

double gravity_value(std::vector<double> x, bool discrete);

namespace {
double reduce_numeric(std::vector<double> v, int reducer, bool is_class) {
  v.erase(std::remove_if(v.begin(), v.end(), [](double x){ return !R_finite(x); }), v.end());
  if (v.empty()) return NA_REAL;
  std::sort(v.begin(), v.end());
  if (is_class || reducer == 2) {
    double best = v[0]; int bestn = 0;
    for (std::size_t i=0; i<v.size();) { double val=v[i]; int n=0; while(i<v.size() && v[i]==val){++i;++n;} if(n>bestn){bestn=n; best=val;} }
    return best;
  }
  if (reducer == 1) {
    std::size_t n=v.size(); return (n%2) ? v[n/2] : (v[n/2-1]+v[n/2])/2.0;
  }
  if (reducer == 3) return gravity_value(v, false);
  double s=0; for(double x:v) s+=x; return s/static_cast<double>(v.size());
}
int find_interval_one(double x, const NumericVector& b) {
  int n=b.size(); if(n==0) return 0;
  if (x <= b[0]) return 0;
  if (x >= b[n-1]) return n-1;
  return static_cast<int>(std::upper_bound(b.begin(), b.end(), x) - b.begin()) - 1;
}
std::string id_for_row(const NumericMatrix& X, const List& boundaries, int r) {
  std::ostringstream os;
  for (int j=0; j<X.ncol(); ++j) { if(j) os << "."; os << find_interval_one(X(r,j), boundaries[j]); }
  return os.str();
}
}

// [[Rcpp::export]]
IntegerVector NNS_duplicate_column_map_cpp(const NumericMatrix& X) {
  const int n=X.nrow(), p=X.ncol();
  IntegerVector rep(p);
  std::unordered_map<std::string, std::vector<int> > seen;
  for(int j=0;j<p;++j){
    std::string key; key.reserve(n*sizeof(double));
    for(int i=0;i<n;++i){ double v=X(i,j); key.append(reinterpret_cast<const char*>(&v), sizeof(double)); }
    int r=j;
    auto it=seen.find(key);
    if(it!=seen.end()){
      for(int cand: it->second){ bool same=true; for(int i=0;i<n && same;++i) same = (X(i,j)==X(i,cand)); if(same){r=cand; break;} }
      it->second.push_back(j);
    } else seen[key]=std::vector<int>(1,j);
    rep[j]=r+1;
  }
  return rep;
}

// [[Rcpp::export]]
List NNS_mreg_setup_cpp(const NumericMatrix& X, const NumericVector& y,
                        const List& boundaries, int reducer_code, bool is_class) {
  const int nr=X.nrow(), p=X.ncol();
  if(y.size()!=nr) stop("y length must equal nrow(X)");
  CharacterVector ids(nr);
  for(int i=0;i<nr;++i) ids[i]=id_for_row(X,boundaries,i);
  std::vector<std::string> unique;
  unique.reserve(nr);
  for(int i=0;i<nr;++i) unique.push_back(as<std::string>(ids[i]));
  std::sort(unique.begin(), unique.end()); unique.erase(std::unique(unique.begin(), unique.end()), unique.end());
  NumericMatrix rpm(unique.size(), p+1);
  for(std::size_t g=0; g<unique.size(); ++g){
    for(int j=0;j<p;++j){ std::vector<double> vals; for(int i=0;i<nr;++i) if(as<std::string>(ids[i])==unique[g]) vals.push_back(X(i,j)); rpm(g,j)=reduce_numeric(vals,reducer_code,false); }
    std::vector<double> yy; for(int i=0;i<nr;++i) if(as<std::string>(ids[i])==unique[g]) yy.push_back(y[i]); rpm(g,p)=reduce_numeric(yy,reducer_code,is_class);
  }
  CharacterVector row_ids(unique.size()); for(std::size_t i=0;i<unique.size();++i) row_ids[i]=unique[i];
  return List::create(_["RPM"] = rpm, _["ids"] = ids, _["row_ids"] = row_ids);
}
