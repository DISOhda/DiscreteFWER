#include <RcppArmadillo.h>
using namespace Rcpp;

inline void eval_pv(double& eval, double val, const NumericVector& vec, int len, int& pos){
  //if(val < 1){
    while(pos < len && vec[pos] <= 1 && vec[pos] <= val) pos++;
    if(pos) eval = vec[pos - 1];
    else eval = 0;
  //} else eval = 1;
}

// computes the index of the largest element of a vector which is <= a given value
inline int binary_search(const NumericVector &vec, const double value, const int len) {
  int idx_left = 0, idx_right = len - 1, idx_mid = len - 1;
  bool stop = false;
  
  while(!stop) {
    if(vec[idx_mid] > value) {
      if(idx_mid == 0) {
        stop = true;
      } else if(idx_mid - idx_left == 1) {
        stop = true;
        idx_mid = idx_left;
      } else {
        idx_right = idx_mid;
        idx_mid = idx_left + (idx_right - idx_left) / 2;
      }
    } else if(vec[idx_mid] <= value) {
      if(vec[idx_mid] == value || idx_mid == len - 1 || idx_right - idx_mid == 1) {
        stop = true;
      } else {
        idx_left = idx_mid;
        idx_mid = idx_left + (idx_right - idx_left + 1) / 2;
      }
    }
  }
  
  return idx_mid;
}

// sort order
IntegerVector order(const NumericVector &x, bool descending = false) {
  arma::vec y = as<arma::vec>(x);
  IntegerVector ord = as<IntegerVector>(wrap(arma::sort_index(y)));
  
  if(descending)
    return rev(ord);
  else
    return ord;
}
