#include <RcppArmadillo.h>
using namespace Rcpp;

//' A matrix multiplication function written by RCpp to accelerate computation.
//' @param matrixA A matrix
//' @param matrixB Another matrix
//'
//' @return the matrix multiplication of A and B
//'
//' @export
// [[Rcpp::export]]
arma::mat matmult(const arma::mat& A, const arma::mat& B) {
  return A * B; 
}