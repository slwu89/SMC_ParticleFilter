#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <boost/function.hpp>
using namespace Rcpp;

/*
 * sir_rate is function to return event rates based on current system state
 * theta: vector of parameters
 * current_state: vector of state variables
 */
// [[Rcpp::export]]   
arma::vec sir_rate(arma::vec theta, arma::vec current_state){
  
}

