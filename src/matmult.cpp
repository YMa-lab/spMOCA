#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat matmult(const arma::mat& A, const arma::mat& B) {
  return A * B; 
}