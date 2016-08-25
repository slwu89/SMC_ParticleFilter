#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
#include <boost/function.hpp>
using namespace Rcpp;


/*
 * sir_rate is function to return event rates based on current system state for SIR model with demography
 * theta: vector of parameters
 * current_state: vector of state variables
 */
// [[Rcpp::export]]   
arma::vec sir_demography_rate(arma::vec theta, arma::vec current_state){
  
  //define parameters
  double beta = theta(0) / theta(1);
  double gamma = 1 / theta(1);
  double mu = 1 / theta(2);
  
  //define states
  int s = current_state(0);
  int i = current_state(1);
  int r = current_state(2);
  int n = s + i + r;
  
  //return rates
  arma::vec rates = arma::zeros(6);
  rates(0) = mu * n; //birth into S
  rates(1) = beta * s * i/n; //S to I
  rates(2) = gamma * i; //I to R
  rates(3) = mu * s; //death from S
  rates(4) = mu * i; //death from I
  rates(5) = mu * r; //death from R
  return(rates);
}


// define wrappers to allow any user specified rate function to be used in the SSA algorithm
typedef boost::function<arma::vec(arma::vec,arma::vec)> rate_arguments;

arma::vec rate_wrapper(arma::vec theta, arma::vec current_state, rate_arguments function){
  return(function(theta,current_state));
}

// // [[Rcpp::export]]
// void test(arma::vec theta, arma::vec init_state){
//   Rcout << "testing passing function" << rate_wrapper(theta,init_state,sir_demography_rate) << std::endl;
// }


/*
 * gillespie_first is a function to simulate a realization of a continuous time Markov process via the first-reaction method
 * theta: vector of parameters
 * init_state: vector of state variables
 * trans: matrix describing transitions
 * t_end: integer value of simulation length
 * rates: function to evaluate rates given current state and theta
 */
// [[Rcpp::export]]
arma::mat gillespie_first(arma::vec theta, arma::vec init_state, arma::mat trans, int t_end, bool info = false){

  int n_event = trans.n_rows;
  
  //initialize trace of Monte Carlo simulation
  arma::mat trace = arma::zeros(1,init_state.n_elem);
  trace.row(0) = init_state.t();
  arma::vec current_state;
  
  //main simulation loop
  double time = 0.0;
  int i = 1;
  
  while(time <= t_end){
    
    if(info){ //print simulation details
      Rcout << "time: " << time << ", i: " << i << std::endl;  
    }
    
    current_state = trace.row(i-1).t(); //extract state at beginning of time jump
    arma::vec current_rates = rate_wrapper(theta,current_state,sir_demography_rate); //calculate current rates
    
    arma::vec rand; //create vector of U[0,1]
    rand.randu(n_event);
    
    arma::vec tau(n_event); //calculate vector of times to next event
    for(int j=0; j<n_event; j++){
      tau(j) = (-1 / current_rates(j)) * log(rand(j));
    }
    
    int first_rxn = tau.index_min(); //which even occurs first
    current_state = current_state + trans.row(first_rxn).t(); //update the current state
    trace.insert_rows(i,current_state.t());
    
    time = time + tau(first_rxn); //update time
    i++; //update iterator
  }
  
  return(trace);
}

/***R
#specify stoichiometry/events
trans <- matrix(c(1,0,0, #birth into S
                  -1,1,0, #infection from S into I
                  0,-1,1, #recovery from I into R
                  -1,0,0, #death from S
                  0,-1,0, #death from I
                  0,0,-1), #death from R
                nrow=6,ncol=3,byrow=TRUE)

#specify theta (R0, infectious duration, lifespan)
theta <- c(2.5,5,365*65)

#specifc initial state vector
init_state <- c(1e3,1,0)

sim_firstRxn <- gillespie_first(theta,init_state,trans,20,TRUE)
*/


/*
 * gillespie_direct is a function to simulate a realization of a continuous time Markov process via the direct reaction method
 * theta: vector of parameters
 * init_state: vector of state variables
 * trans: matrix describing transitions
 * t_end: integer value of simulation length
 * rates: function to evaluate rates given current state and theta
 */
arma::mat gillespie_direct(arma::vec theta, arma::vec init_state, arma::mat trans, int t_end, bool info = false){
  
  int n_event = trans.n_rows;
  
  //initialize trace of Monte Carlo simulation
  arma::mat trace = arma::zeros(1,init_state.n_elem);
  trace.row(0) = init_state.t();
  arma::vec current_state;
  
  //main simulation loop
  double time = 0.0;
  int i = 1;
  
  while(time <= t_end){
    
    if(info){ //print simulation details
      Rcout << "time: " << time << ", i: " << i << std::endl;  
    }
    
    current_state = trace.row(i-1).t(); //extract state at beginning of time jump
    arma::vec current_rates = rate_wrapper(theta,current_state,sir_demography_rate); //calculate current rates
    
    double w0 = sum(current_rates); //sum of rate (propensity) functions
    double tau = 1/w0 * log(1/R::runif(0,1)); //calculate time to next reaction
    
    double r_num = R::runif(0,1);
    double w0_rxn = r_num * w0; //instantaneous event probabilities
    
    //decide which event j occured
    int j = 0;
    while(sum(current_rates.subvec(0,j)) < w0_rxn){
      j = j + 1;
    }
    
    current_state = current_state + trans.row(j).t(); //update the current state
    trace.insert_rows(i,current_state.t());
    
    time = time + tau; //update time
    i++; //update iterator
  }
  
  return(trace);
}
  
