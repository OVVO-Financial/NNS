// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// forward decls to your central_tendencies (or include the header you made)
SEXP NNS_mode_cpp(SEXP xSEXP, bool discrete, bool multi);
SEXP NNS_gravity_cpp(SEXP xSEXP, bool discrete);

// mean/median helpers
static inline double mean_no_na(const NumericVector& v){
  long double s=0.0L; std::size_t m=0;
  for(double xi: v) if(R_finite(xi)){ s+=xi; ++m; }
  return m? static_cast<double>(s/m): NA_REAL;
}
static inline double median_no_na(NumericVector v){
  std::vector<double> a; a.reserve(v.size());
  for(double xi: v) if(R_finite(xi)) a.push_back(xi);
  if(a.empty()) return NA_REAL;
  std::size_t n=a.size();
  std::nth_element(a.begin(), a.begin()+n/2, a.end());
  double hi=a[n/2];
  if(n&1u) return hi;
  auto lm = std::max_element(a.begin(), a.begin()+n/2);
  return (*lm+hi)*0.5;
}

struct Agg{
  std::string noise;
  inline double mode_disc_single(const NumericVector& v) const{
    return as<double>(NNS_mode_cpp(v, true, false));
  }
  inline double gravity_cont(const NumericVector& v, bool discrete=false) const{
    return as<double>(NNS_gravity_cpp(v, discrete));
  }
  inline double for_x(const NumericVector& v) const{
    if(noise=="mean") return mean_no_na(v);
    if(noise=="median") return median_no_na(v);
    if(noise=="mode") return mode_disc_single(v);
    if(noise=="mode_class") return gravity_cont(v,false);
    return gravity_cont(v,false);
  }
  inline double for_y(const NumericVector& v) const{
    if(noise=="mean") return mean_no_na(v);
    if(noise=="median") return median_no_na(v);
    if(noise=="mode") return mode_disc_single(v);
    if(noise=="mode_class") return mode_disc_single(v);
    return gravity_cont(v,false);
  }
};

struct Pair{ double x; double y; };

// [[Rcpp::export]]
List NNS_part_cpp(NumericVector x,
                  NumericVector y,
                  Nullable<std::string> type,
                  Nullable<int> order_in,
                  int obs_req,
                  bool min_obs_stop,
                  std::string noise_reduction){
  
  const int n=x.size();
  if(y.size()!=n) stop("x and y must have same length");
  
  int default_order = std::max((int)std::ceil(std::log2(std::max(1,n))),1);
  int max_order = order_in.isNotNull()? as<int>(order_in): default_order;
  if(max_order==0) max_order=1;
  bool xonly = type.isNotNull();
  std::transform(noise_reduction.begin(), noise_reduction.end(),
                 noise_reduction.begin(), ::tolower);
  Agg agg{noise_reduction};
  
  std::vector<std::string> quadrant(n,"q"), prior_quadrant(n,"pq");
  int depth=0;
  
  // -------------------- NEW: segment collectors --------------------
  std::vector<double> H_x0, H_x1, H_y;      // horizontal: (x0, x1) at y
  std::vector<double> V_x,  V_y0, V_y1;     // vertical: x at (y0, y1)
  std::vector<double> V_lines;              // XONLY: abline(v = ...)
  // ----------------------------------------------------------------
  
  while(true){
    if(depth>=max_order) break;
    if(depth>=(int)std::floor(std::log2(std::max(1,n)))) break;
    
    // groups
    std::unordered_map<std::string,std::vector<int>> grp; grp.reserve(n*2);
    for(int i=0;i<n;++i) grp[quadrant[i]].push_back(i);
    
    // choose to split
    std::vector<std::string> to_split; to_split.reserve(grp.size());
    for(auto &kv: grp) if((int)kv.second.size()>obs_req) to_split.push_back(kv.first);
    if(to_split.empty()) break;
    
    // centers per split group (and bounds if we draw)
    std::unordered_map<std::string, Pair> centers; centers.reserve(to_split.size());
    for(const auto &q: to_split){
      const auto &idx = grp[q];
      NumericVector xv(idx.size()), yv(idx.size());
      double minx=R_PosInf, maxx=R_NegInf, miny=R_PosInf, maxy=R_NegInf; // NEW
      for(std::size_t k=0;k<idx.size();++k){
        int i=idx[k]; xv[k]=x[i]; yv[k]=y[i];
        if(R_finite(x[i])){ if(x[i]<minx)minx=x[i]; if(x[i]>maxx)maxx=x[i]; } // NEW
        if(R_finite(y[i])){ if(y[i]<miny)miny=y[i]; if(y[i]>maxy)maxy=y[i]; } // NEW
      }
      Pair c{ agg.for_x(xv), agg.for_y(yv) };
      centers[q]=c;
      
      // -------------------- NEW: record segments for 2D case --------------------
      if(!xonly){
        if(R_finite(c.y) && R_finite(minx) && R_finite(maxx)){
          H_x0.push_back(minx); H_x1.push_back(maxx); H_y.push_back(c.y);
        }
        if(R_finite(c.x) && R_finite(miny) && R_finite(maxy)){
          V_x.push_back(c.x); V_y0.push_back(miny); V_y1.push_back(maxy);
        }
      }
      // -------------------------------------------------------------------------
    }
    
    // -------------------- NEW: record ablines for XONLY at this depth ----------
    if(xonly){
      for(auto &kv: grp){
        const auto &idx = kv.second;
        double minx=R_PosInf, maxx=R_NegInf;
        for(int i: idx){
          if(R_finite(x[i])){ if(x[i]<minx)minx=x[i]; if(x[i]>maxx)maxx=x[i]; }
        }
        if(R_finite(minx)) V_lines.push_back(minx);
        if(R_finite(maxx)) V_lines.push_back(maxx);
      }
    }
    // ---------------------------------------------------------------------------
    
    // reassign
    for(const auto &q: to_split){
      const Pair c = centers[q];
      for(int i: grp[q]){
        prior_quadrant[i]=quadrant[i];
        int qn;
        if(!xonly){
          int lox = (R_finite(x[i])&&R_finite(c.x)) ? (x[i]<=c.x) : 0;
          int loy = (R_finite(y[i])&&R_finite(c.y)) ? (y[i]<=c.y) : 0;
          qn = 1 + lox + 2*loy;
        }else{
          int lox = (R_finite(x[i])&&R_finite(c.x)) ? (x[i]>c.x) : 0;
          qn = 1 + lox;
        }
        quadrant[i]+=std::to_string(qn);
      }
    }
    
    ++depth;
    
    if(min_obs_stop){
      std::unordered_map<std::string,int> cnt; cnt.reserve(n*2);
      for(const auto &qstr: quadrant) ++cnt[qstr];
      int minc=n; for(auto &kv: cnt) if(kv.second<minc) minc=kv.second;
      if(minc<=obs_req) break;
    }
  }
  
  // assemble results
  CharacterVector q_cur(n), q_prior(n);
  for(int i=0;i<n;++i){ q_cur[i]=quadrant[i]; q_prior[i]=prior_quadrant[i]; }
  DataFrame part = DataFrame::create(_["x"]=x,_["y"]=y,_["quadrant"]=q_cur,
                                     _["prior.quadrant"]=q_prior,
                                     _["stringsAsFactors"]=false);
  std::unordered_map<std::string,std::vector<int>> by_prior; by_prior.reserve(n*2);
  for(int i=0;i<n;++i) by_prior[prior_quadrant[i]].push_back(i);
  
  std::vector<std::string> rp_q; std::vector<double> rp_x, rp_y;
  for(auto &kv: by_prior){
    const auto &idx=kv.second;
    NumericVector xv(idx.size()), yv(idx.size());
    for(std::size_t k=0;k<idx.size();++k){ int i=idx[k]; xv[k]=x[i]; yv[k]=y[i]; }
    rp_q.push_back(kv.first);
    rp_x.push_back(agg.for_x(xv));
    rp_y.push_back(agg.for_y(yv));
  }
  DataFrame rp = DataFrame::create(_["quadrant"]=wrap(rp_q),
                                   _["x"]=wrap(rp_x),
                                   _["y"]=wrap(rp_y),
                                   _["stringsAsFactors"]=false);
  
  // -------------------- NEW: return segments for plotting ----------------------
  DataFrame seg_h = DataFrame::create(_["x0"]=wrap(H_x0),
                                      _["x1"]=wrap(H_x1),
                                      _["y"]= wrap(H_y),
                                      _["stringsAsFactors"]=false);
  DataFrame seg_v = DataFrame::create(_["x"]= wrap(V_x),
                                      _["y0"]=wrap(V_y0),
                                      _["y1"]=wrap(V_y1),
                                      _["stringsAsFactors"]=false);
  // ---------------------------------------------------------------------------
  
  return List::create(_["order"]=depth,
                      _["dt"]=part,
                      _["regression.points"]=rp,
                      _["segments_h"]=seg_h,      // NEW
                      _["segments_v"]=seg_v,      // NEW
                      _["vlines"]=wrap(V_lines)); // NEW
}
