#include "kernel.h"

NumericVector kernel_DFWER_single_fast(const List &pCDFlist, const NumericVector &pvalues, const bool independent, const Nullable<NumericVector> &pCDFcounts) {
  // Number of p-values
  int numValues = pvalues.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // counts of the CDFs
  NumericVector CDFcounts;
  if(numCDF == numValues || pCDFcounts.isNull()) 
    CDFcounts = NumericVector(numCDF, 1.0);
  else 
    CDFcounts = pCDFcounts;
  
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
  // number of all attainable p-values in the support
  int numValues = support.length();
  
  // extract p-value CDF vectors
  NumericVector* sfuns = new NumericVector[(unsigned int)numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // get count of each unique p-value distribution
  NumericVector CDFcounts;
  if(pCDFcounts.isNull()) CDFcounts = NumericVector(numCDF, 1.0);
  else CDFcounts = pCDFcounts;
  
  // restrict support to values <= alpha (critical value cannot exceed alpha)
  int index_max = binary_search(support, alpha, numValues);
  //NumericVector pv_list = support[Range(0, index_max)];
  
  // transform support with fast kernel
  NumericVector support_transf = kernel_DFWER_single_fast(pCDFlist, support, independent, CDFcounts);
  
  // get index of critical value
  int idx_pval = binary_search(support_transf, alpha, index_max + 1);
  // vector to store critical value
  NumericVector crit(1, support[idx_pval]);
  
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

NumericVector kernel_DFWER_sd_fast(const List &pCDFlist, const NumericVector &sorted_pv, const bool independent, const Nullable<List> &pCDFindices) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // indices of the CDFs and their counts
  IntegerVector* CDFindices = new IntegerVector[numCDF];
  int* CDFcounts = new int[numCDF];
  if(pCDFindices.isNull()) {
    for(int i = 0; i < numCDF; i++) {
      CDFindices[i] = IntegerVector(1, i + 1);
      CDFcounts[i] = 1;
    }
  } else {
    for(int i = 0; i < numCDF; i++) {
      CDFindices[i] = as<IntegerVector>(as<List>(pCDFindices)[i]);
      CDFcounts[i] = CDFindices[i].length();
    }
  }
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[numCDF];
  for(int i = 0; i < numCDF; i++) sfuns[i] = as<NumericVector>(pCDFlist[i]);
  
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // evaluation of current p-value CDF
  NumericVector f_eval(numTests);
  for(int i = 0; i < numCDF; i++) {
    checkUserInterrupt();
    
    // current position in i-th CDF
    int pos = 0;
    // length of i-th CDF
    int len = sfuns[i].length();
    // current sorted p-value to which i-th CDF belongs
    int k = 0;
     
    for(int j = 0; j < CDFindices[i][CDFcounts[i] - 1]; j++) {
      // evaluate i-th CDF for all RELEVANT p-values
      eval_pv(f_eval[j], sorted_pv[j], sfuns[i], len, pos);
      // multiply with number of relevant p-values
      f_eval[j] *= CDFcounts[i] - k;
      if(CDFindices[i][k] == j + 1) k++;
    }
    for(int j = CDFindices[i][CDFcounts[i] - 1]; j < numTests; j++)
      f_eval[j] = 0;
    
    if(independent)
      // add logs instead of multiplication of non-logs for numerical stability
      pval_transf += log(1 - f_eval);
    else 
      pval_transf += f_eval;
  }
  
  if(independent)
    // revert log
    pval_transf = 1 - exp(pval_transf);
  
  // garbage collection
  delete[] sfuns;
  
  return pval_transf;
}

List kernel_DFWER_sd_crit(const List &pCDFlist, const NumericVector &support, const NumericVector &sorted_pv, const double alpha, const bool independent, const Nullable<List> &pCDFindices) {
  // number of tests
  int numTests = sorted_pv.length();
  // number of unique p-value distributions
  int numCDF = pCDFlist.length();
  // support size
  int numValues = support.length();
  
  // extract p-value CDF vectors and their lengths
  NumericVector* sfuns = new NumericVector[numCDF];
  int* lens = new int[numCDF];
  for(int i = 0; i < numCDF; i++) {
    sfuns[i] = as<NumericVector>(pCDFlist[i]);
    lens[i] = sfuns[i].length();
  }
  
  // indices of the CDFs and their counts
  IntegerVector* CDFindices = new IntegerVector[numCDF];
  int* CDFcounts = new int[numCDF];
  //int* CDFcounts_crit = new int[numCDF];
  if(pCDFindices.isNull()) {
    for(int i = 0; i < numCDF; i++) {
      CDFindices[i] = IntegerVector(1, i + 1);
      CDFcounts[i] = 1;
      //CDFcounts_crit[i] = 1;
    }
  } else {
    for(int i = 0; i < numCDF; i++) {
      CDFindices[i] = as<IntegerVector>(as<List>(pCDFindices)[i]);
      CDFcounts[i] = CDFindices[i].length();
      //CDFcounts_crit[i] = CDFcounts[i];
    }
  }
  
  // critical values indices
  NumericVector crit(numTests);
  // vector to store transformed p-values
  NumericVector pval_transf(numTests);
  // possibly large data size requires to use chunks
  // size of the chunks (i.e. number of elements in a ~512 MiB matrix)
  int size = std::max<int>(1, std::pow(2.0, 26.0) / numCDF);
  // number of chunks
  int chunks = (numValues - 1) / size + 1;
  
  // index of current critical value
  int idx_crit = 0;
  // index of current raw p-value to be transformed
  int idx_transf = 0;
  // last positions in step function evaluations
  int* pos = new int[numCDF]{};
  
  // compute critical values (and transformed raw p-values for step-down)
  for(int i = 0; i < chunks; i++) {
    // the min( , numValues) is here for the last chunk
    NumericVector pv = support[Range(i * size, std::min<int>((i + 1) * size, numValues) - 1)];
    // length of the vector
    int numPV = pv.length();
    // rows:    indices from 1 to numTests
    // columns: p-values
    NumericMatrix mat(numCDF, numPV);
    // compute columns \sum_{j=1}^numTests F_(j)(pv)
    for(int j = 0; j < numCDF; j++) {
      checkUserInterrupt();
      for(int k = 0; k < numPV; k++)
        eval_pv(mat(j, k), pv[k], sfuns[j], lens[j], pos[j]);
    }
    
    // compute transformed p-value support (as in pv_list)
    int j = 0;
    while(j < numPV && (idx_transf < numTests || idx_crit < numTests)) {
      checkUserInterrupt();
      
      // index of current p-value
      int idx_pval = i * size + j; Rcout << idx_pval << " " << idx_crit << " " << pv[j];
      // store column values in a vector
      NumericVector temp = NumericVector(mat(_, j));
      // get order
      //IntegerVector ord = order(temp, true);
      IntegerVector ord = order(temp, false);
      
      // number of remaining needed values
      int rem = numTests - idx_crit - 1;
      // sum
      double s = 0;
      // compute sum
      if(idx_crit < numTests) {
        //int k = 1;
        int k = 0;
        while(k < numCDF && CDFcounts[ord[k]] <= rem && s <= alpha - pv[j]) {
          s += CDFcounts[ord[k]] * temp[ord[k]];
          rem -= CDFcounts[ord[k]];
          k++;
        }
        if(k < numCDF && rem > 0 && s <= alpha - pv[j]) s += rem * temp[ord[k]];
        
        Rcout << ", " << s + pv[j];
        
        // check satisfaction of condition
        if(s <= alpha - pv[j]) {
          // current p-value satisfies condition
          // => go to next p-value in this chunk
          j++;
        } else {
          Rcout << " => Here!";
          // current p-value does not satisfy condition
          // => save previous p-value as critical value
          if(idx_pval) crit[idx_crit] = support[idx_pval - 1];
          else crit[idx_crit] = 0;
          // stay at current p-value but check for next critical value
          idx_crit++;
        }
      } else j++;
      Rcout << "\n";
      
      if(idx_transf < numTests) {
        // current p-value satisfies condition or some transforms not done yet
        // => do transformation for sorted observed p-value that match
        while(idx_transf < numTests && support[idx_pval] > sorted_pv[idx_transf]) idx_transf++;
        while(idx_transf < numTests && support[idx_pval] == sorted_pv[idx_transf]) {
          pval_transf[idx_transf] = 0;
          rem = numTests - idx_transf;
          int k = 0;
          while(k < numCDF && rem > 0) {
            int count_k = CDFcounts[ord[k]];
            //if(idx_transf == 1) Rcout << ord[k] << ": " << count_k << " (" << temp[ord[k]] << "),";
            if(count_k >= 1) {
              int l = 0;
              //if(idx_transf == 1) Rcout << " " << CDFindices[ord[k]][l];
              while(l < CDFcounts[ord[k]] && CDFindices[ord[k]][l++] <= idx_transf)
                count_k--;
            }
            //if(idx_transf == 1) Rcout << " => " << count_k << "\n";
            pval_transf[idx_transf] += std::min<int>(count_k, rem) * temp[ord[k]];
            rem -= count_k;
            k++;
          }
          idx_transf++;
        }
      }
    }
  }
  
  // garbage collection
  delete[] pos;
  delete[] sfuns;
  delete[] lens;
  
  // output results
  return List::create(Named("crit.consts") = crit, Named("pval.transf") = pval_transf);
}