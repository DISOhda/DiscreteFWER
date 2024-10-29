#include "kernel.h"

NumericVector kernel_DFWER_single_fast(const List &pCDFlist, const NumericVector &pvalues, const bool independent, const Nullable<NumericVector> &pCDFcounts) {
  // Number of p-values
  int numValues = pvalues.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution and number of tests
  int numTests;
  // counts of the CDFs
  NumericVector CDFcounts;
  if(numCDF == numValues || pCDFcounts.isNull()) {
    CDFcounts = NumericVector(numCDF, 1.0);
    numTests = numCDF;
  } else {
    CDFcounts = pCDFcounts;
    numTests = sum(CDFcounts);
  }
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // vector to store transformed p-values
  NumericVector pval_transf(numValues);
  // evaluation of current p-value CDF
  NumericVector f_eval(numValues);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    int pos = 0;
    int len = sfuns[i].length();
    for(int j = 0; j < numValues; j++) {
      eval_pv(f_eval[j], pvalues[j], sfuns[i], len, pos);
    }
    
    if(independent)
      pval_transf += CDFcounts[i] * log(1 - f_eval);
    else 
      pval_transf += CDFcounts[i] * f_eval;
  }
  
  if(independent)
    pval_transf = 1 - exp(pval_transf);
  
  // garbage collection
  delete[] sfuns;
  
  return pval_transf;
}

List kernel_DFWER_single_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const double alpha, const bool independent, const Nullable<NumericVector> &pCDFcounts) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // transform support with fast kernel
  NumericVector support_transf = kernel_DFWER_single_fast(pCDFlist, support, independent, CDFcounts);
  
  // number of all attainable p-values in the support
  int numValues = support.length();
  // vector to store critical values indices
  NumericVector crit(numTests);
  
  // get indices of critical values
  int idx_pval = 0;
  for(int i = 0; i < numTests; i++) {
    checkUserInterrupt();
    while(idx_pval < numValues && support_transf[idx_pval] <= alpha) idx_pval++;
    crit[i] = support[idx_pval - 1];
  }
  
  // store transformed sorted pvalues
  NumericVector pval_transf(numTests);
  // search for sorted p-values in 'pv_list' and get their transformations
  idx_pval = 0;
  for(int i = 0; i < numTests; i++) {
    checkUserInterrupt();
    while(idx_pval < numValues && support[idx_pval] != sorted_pv[i]) idx_pval++;
    pval_transf[i] = support_transf[idx_pval];
  }
  
  // garbage collection
  delete[] sfuns;
  
  // return critical values and transformed sorted p-values
  return List::create(Named("crit.consts") = crit, Named("pval.transf") = pval_transf);
}
