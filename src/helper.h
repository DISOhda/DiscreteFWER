#include <Rcpp.h>
using namespace Rcpp;

inline void eval_pv(double& eval, double val, const NumericVector& vec, int len, int& pos){
  //if(val < 1){
    while(pos < len && vec[pos] <= 1 && vec[pos] <= val) pos++;
    if(pos) eval = vec[pos - 1];
    else eval = 0;
  //} else eval = 1;
}
