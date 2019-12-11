#include <Rcpp.h>
using namespace Rcpp;

// These functions compute the CDF of the sum of a sequence of
// independent binomial trials with the same index N. 

// function for calculating the CDF of observing y_i (and y_i - 1) individuals 
// at site i by summing over the binomial distribution for all possible 
// combinations of observations across replicate surveys at site i (and it's the 
// CDF, so observing 0 up to y_i and 0 up to y_i - 1)

// [[Rcpp::export]]
NumericVector pbinsumRow(NumericVector y, double N, NumericVector p) {
   // empty vector to hold the three "residual" values
   NumericVector res(3);
   // empty vector to hold the 1 total number of observations at a site
   NumericVector ySum(1);
   // start with -1 observations at the site
   ySum(0) = -1;
   // Handle NAs by setting p to 0 and y  to 1 => dbinom(y,N,p) = 0.
   for (int i=0; i<y.length(); i++) {
      if(NumericVector::is_na(y[i]) || NumericVector::is_na(p[i])) {
         y[i] = 1;
         p[i] = 0.0;
      } else {
         if (ySum[0] == -1) {
            ySum[0] = y[i];
         } else {
            // add up all the observations at the site (across replicates)
            ySum[0] += y[i];
         }
      }
   }
   // if there weren't any observations, return NA
   if (ySum[0] == -1) {
      res[0] = res[1] = res[2] = NA_REAL;
      return res;
   }
   // if was an infinite number of observations, return NAs
   if (!Rcpp::traits::is_finite<REALSXP>(ySum(0))) {
      res[0] = ySum[0];
      res[1] = res[2] = R_NaN;
      return res;
   }
   // create an empty matrix to hold probabilities of observing x observations at each replicate
   NumericMatrix pMat(ySum[0] + 1, p.length());
   // vector for summing probabilities across k values
   NumericVector pT(ySum[0] + 1);
   // another vector for summing probabilities across k values
   NumericVector pF(ySum[0] + 1);
   
   // for each replicate
   for (int col=0; col<pMat.ncol(); col++) {
      // for each possible number of observations at an individual replicate (i.e. 0 up to ySum)
      for (int k = 0; k <= ySum[0]; k++) {
         // calculate the probability of observing that many individuals given true abundance of 
         // N and the detection probability for that site and replicate
         pMat(k, col) = Rf_dbinom(k, N, p[col], false);
      }
   }
   for (int k = 0; k <= ySum[0]; k++) {
      // start pT as the first column of probabilities
      pT[k] = pMat(k, 0);
      // start pF as the first column of probabilities
      pF[k] = pT[k];
   }
   // for each replicate
   for (int col = 1; col<pMat.ncol(); col++) {
      // for each possible number of observations at the replicate
      for (int j = 0; j <= ySum[0]; j++) {
         // start with pF at 0
         pF[j] = 0;
         // for that k value, cycle through all k values up to that k value and 
         // add up the probabilities of the collective observations (across all 
         // replicates up to col)
         for (int i = 0; i<=j; i++) {
            pF[j] += pT[i]*pMat(j-i, col);
         }
      }
      // before starting on a new column (replicate) save the collective probabilities
      // to the pT vector
      for (int k = 0; k <= ySum[0]; k++) {
         pT[k] = pF[k];
      }
   }
   // the last for loop goes over all the replicates and returns a vector (pF)
   // that has the total probability of observing exactly k individuals (that could
   // be distributed across col replicates).
   for (int k = 0; k < ySum[0]; k++)
      res[1] += pF[k]; // add up all the cumulative probability of observing up to (ySum - 1) total individuals
   
   res[2] = res[1] + pF[ySum[0]];
   // cumulative probability of observing ySum
   res[0] = ySum[0];
   return res; 
}  

// [[Rcpp::export]]
NumericMatrix pbinsum(NumericMatrix y, NumericVector N, NumericMatrix p) {
   NumericMatrix cumProb(y.nrow(),3);
   if(y.ncol() != p.ncol() || y.nrow() != p.nrow()) {
      stop("Dimensions of y do not match those of p.");
   }
   if(y.nrow() != N.length()) {
      stop("Length of N does not match the number of rows of y or p.");
   }
   /*
    for (int i = 0; i < y.nrow(); i++) {
    for (int j = 0; j < y.ncol(); j++) {
    if (NumericVector::is_na(y(i,j)) || NumericVector::is_na(p(i,j))) {
    y(i,j) = NA_REAL;
    p(i,j) = NA_REAL;
    }
    }
    }
    */
   // loop over all sites and calculate the cumulative probabilities 
   for (int i = 0; i < y.nrow(); i++) {
      cumProb.row(i) = pbinsumRow(y.row(i), N[i], p.row(i));
   }
   return cumProb;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
 pbinsumRow(c(10,2,2), 10, c(.5,.4,.3))
 pbinsumRow(c(NA, 3, 4), 11, c(.1,.1, NA))
 pbinsum(matrix(c(10,NA,2,3,2,4),nrow = 2), c(10,11), matrix(c(.5,.1,.4,.1,.3,NA), nrow = 2))
 */