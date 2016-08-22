#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <boost/function.hpp>
using namespace Rcpp;

/*
 * sir_rate is function to return event rates based on current system state for SIR model
 * theta: vector of parameters
 * current_state: vector of state variables
 */
// [[Rcpp::export]]   
arma::vec sir_rate(arma::vec theta, arma::vec current_state){
  
  //define parameters
  double beta = theta(0) / theta(1);
  double gamma = 1/theta(1);
  
  //define states
  int s = current_state(0);
  int i = current_state(1);
  int r = current_state(2);
  int n = s + i + r;
  
  //return rates
  arma::vec rates = arma::zeros(2);
  rates(0) = beta * s * i/n; //susceptible to infectious
  rates(1) = gamma * i; //infectious to recovered
  return(rates);
}

// [[Rcpp::export]]
NumericVector sir_rateRcpp(NumericVector theta, IntegerVector current_state){
  double beta = theta[0] / theta[1];
  double gamma = 1 / theta[1];
  int s = current_state[0];
  int i = current_state[1];
  int r = current_state[2];
  int n = s + i + r;
  NumericVector out(2);
  out[0] = beta * s * i/n;
  out[1] = gamma * i;
  return(out);
}


/*
 * gillespie_frm is a function to simulate a realization of a continuous time Markov process via the first-reaction method
 * theta: vector of parameters
 * init_state: vector of state variables
 * t_end: integer value of simulation length
 * trans: matrix describing transitions
 * rates: function to evaluate rates given current state and theta
 */
// [[Rcpp::export]]
arma::vec test(arma::vec theta1, arma::vec current_state1, boost::function<NumericVector (NumericVector, NumericVector)> rates){
  arma::vec out;
  out = rates(theta1,current_state1);
  return(out);
}