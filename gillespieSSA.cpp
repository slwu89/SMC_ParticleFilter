#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <boost/function.hpp>
using namespace Rcpp;


/*
 * This file gives functions to run exact discrete event (integer valued) simulations of stochastic processes.
 * These processes are three variants on the Gillespie algorithm (SSAs):
 * Gillespie's Direct Method
 * Gillespie's First Reaction Method
 * Gibson & Bruck's Next Reaction Method
 * 
 * You will need to install Rcpp, RcppArmadillo, RcppProgress, BH, and have a working C++ compiler to use them.
 */


/*
 * sir_demography_rate is function to return event rates based on current system state for SIR model with demography
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


/*
 * sir_rate is function to return event rates based on current system state for SIR model without demography
 * theta: vector of parameters
 * current_state: vector of state variables
 */
// [[Rcpp::export]]
arma::vec sir_rate(arma::vec theta, arma::vec current_state){
  
  //define parameters
  double beta = theta(0) / theta(1);
  double gamma = 1 / theta(1);
  
  //define states
  int s = current_state(0);
  int i = current_state(1);
  int r = current_state(2);
  int n = s + i + r;
  
  //return rates
  arma::vec rates = arma::zeros(2);
  rates(0) = beta * s * i/n; //S to I
  rates(1) = gamma * i; //I to R
  
  return(rates);
}


/*
 * seirs_demography_rate is a function to return event rates based on current system state for SEIRS model with demography
 * theta: vector of parameters
 * current_state: vector of state variables
 */
// [[Rcpp::export]]
arma::vec seirs_demography_rate(arma::vec theta, arma::vec current_state){
  
  //define parameters
  double beta = theta(0) / theta(2);
  double tau = 1 / theta(1);
  double rec = 1 / theta(2);
  double gamma = 1 / theta(3);
  double mu = 1 / theta(4);
  
  //define states
  int s = current_state(0);
  int e = current_state(1);
  int i = current_state(2);
  int r = current_state(3);
  int n = s + e + i + r;
  
  //return rates
  arma::vec rates = arma::zeros(9);
  rates(0) = beta * s * i/n; //S to E
  rates(1) = tau * e; //E to I
  rates(2) = rec * i; //I to R
  rates(3) = gamma * r; //R to S
  rates(4) = mu * n; //birth into S
  rates(5) = mu * s;
  rates(6) = mu * e;
  rates(7) = mu * i;
  rates(8) = mu * r;
  
  return(rates);
}


// define wrappers to allow any user specified rate function to be used in the SSA algorithm
typedef boost::function<arma::vec(arma::vec,arma::vec)> rate_arguments;

arma::vec rate_wrapper(arma::vec theta, arma::vec current_state, rate_arguments function){
  return(function(theta,current_state));
}


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
  
  Progress p(0, false); //progress to check for user interrupt
  
  //initialize trace of Monte Carlo simulation
  arma::mat trace = arma::zeros(1,init_state.n_elem);
  trace.row(0) = init_state.t();
  arma::vec current_state;
  
  //main simulation loop
  double time = 0.0;
  int i = 1;
  
  while(time <= t_end){
    
    //check for user abort
    if(Progress::check_abort()){
      Rcout << "User abort detected at time: " << time << ", i: " << i << std::endl;
      return(arma::zeros(2,2));
    }
    
    if(info){ //print simulation details
      Rcout << "time: " << time << ", i: " << i << std::endl;  
    }
    
    current_state = trace.row(i-1).t(); //extract state at beginning of time jump
    arma::vec current_rates = rate_wrapper(theta,current_state,seirs_demography_rate); //calculate current rates
    
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
# #specify stoichiometry/events
# #SIR with demography
# trans_sir <- matrix(c(1,0,0, #birth into S
#                   -1,1,0, #infection from S into I
#                   0,-1,1, #recovery from I into R
#                   -1,0,0, #death from S
#                   0,-1,0, #death from I
#                   0,0,-1), #death from R
#                 nrow=6,ncol=3,byrow=TRUE)
# 
# #specify theta (R0, infectious duration, lifespan)
# theta_sir <- c(2.5,5,365*65)
# 
# #specify initial state vector
# init_state_sir <- c(1e3,1,0)
#
# sim_firstRxn <- gillespie_first(theta_sir,init_state_sir,trans_sir,20,TRUE)

#SEIRS with demography
trans_seirs <- matrix(c(-1,1,0,0, #infection from S to E
                  0,-1,1,0, #progression from E to I
                  0,0,-1,1, #recovery from I to R
                  1,0,0,-1, #loss of immunity from R to S
                  1,0,0,0, #birth to S
                  -1,0,0,0, #death from S
                  0,-1,0,0, #death from E
                  0,0,-1,0, #death from I
                  0,0,0,-1), #death from R
                  nrow=9,ncol=4,byrow=TRUE)

#specify theta (R0, latent duration, infectious duration, immune duration, lifespan)
theta_seirs <- c(5,2,3,365,365*65)

#specify initial state vector
init_state_seirs <- c(1e3,1,0,0)

sim_firstRxn <- gillespie_first(theta_seirs,init_state_seirs,trans_seirs,20,TRUE)
*/


/*
 * gillespie_direct is a function to simulate a realization of a continuous time Markov process via the direct reaction method
 * theta: vector of parameters
 * init_state: vector of state variables
 * trans: matrix describing transitions
 * t_end: integer value of simulation length
 * rates: function to evaluate rates given current state and theta
 */
// [[Rcpp::export]]
arma::mat gillespie_direct(arma::vec theta, arma::vec init_state, arma::mat trans, int t_end, bool info = false){
  
  Progress p(0, false); //progress to check for user interrupt
  
  //initialize trace of Monte Carlo simulation
  arma::mat trace = arma::zeros(1,init_state.n_elem);
  trace.row(0) = init_state.t();
  arma::vec current_state;
  
  //main simulation loop
  double time = 0.0;
  int i = 1;
  
  while(time <= t_end){
    
    //check for user abort
    if(Progress::check_abort()){
      Rcout << "User abort detected at time: " << time << ", i: " << i << std::endl;
      return(arma::zeros(2,2));
    }
    
    if(info){ //print simulation details
      Rcout << "time: " << time << ", i: " << i << std::endl;  
    }
    
    current_state = trace.row(i-1).t(); //extract state at beginning of time jump
    arma::vec current_rates = rate_wrapper(theta,current_state,seirs_demography_rate); //calculate current rates
    
    double w0 = sum(current_rates); //sum of rate (propensity) functions
    double tau = 1/w0 * log(1/R::runif(0,1)); //calculate time to next reaction
    
    double r_num = R::runif(0,1);
    double w0_rxn = r_num * w0; //instantaneous event probabilities
    
    //decide which event j occured
    int j = 0;
    while(sum(current_rates.subvec(0,j)) < w0_rxn){
      j++;
    }
    
    current_state = current_state + trans.row(j).t(); //update the current state
    trace.insert_rows(i,current_state.t());
    
    time = time + tau; //update time
    i++; //update iterator
  }
  
  return(trace);
}
  
/***R
# #specify stoichiometry/events
# #SIRS with demography
# trans_sir <- matrix(c(1,0,0, #birth into S
#                 -1,1,0, #infection from S into I
#                 0,-1,1, #recovery from I into R
#                 -1,0,0, #death from S
#                 0,-1,0, #death from I
#                 0,0,-1), #death from R
#                 nrow=6,ncol=3,byrow=TRUE)
#   
# #specify theta (R0, infectious duration, lifespan)
# theta_sir <- c(2.5,5,365*65)
#   
# #specify initial state vector
# init_state_sir <- c(1e3,1,0)
#   
# sim_direct <- gillespie_direct(theta_sir,init_state_sir,trans_sir,20,TRUE)

#specify stoichiometry/events
#SEIRS with demography
trans_seirs <- matrix(c(-1,1,0,0, #infection from S to E
                        0,-1,1,0, #progression from E to I
                        0,0,-1,1, #recovery from I to R
                        1,0,0,-1, #loss of immunity from R to S
                        1,0,0,0, #birth to S
                        -1,0,0,0, #death from S
                        0,-1,0,0, #death from E
                        0,0,-1,0, #death from I
                        0,0,0,-1), #death from R
                      nrow=9,ncol=4,byrow=TRUE)

#specify theta (R0, latent duration, infectious duration, immune duration, lifespan)
theta_seirs <- c(5,2,3,365,365*65)

#specify initial state vector
init_state_seirs <- c(1e3,1,0,0)

sim_direct <- gillespie_direct(theta_seirs,init_state_seirs,trans_seirs,20,TRUE)
*/


/*
 * individual rate calculations used for the Gillespie Next-Reaction Method
 */
double s_to_e(arma::vec theta, arma::vec current_state){
  double beta = theta(0) / theta(2);
  int s = current_state(0);
  int i = current_state(2);
  int n = sum(current_state);
  double out = beta * s * i/n;
  return(out);
}

double e_to_i(arma::vec theta, arma::vec current_state){
  double tau = 1 / theta(1);
  int e = current_state(1);
  double out = tau * e;
  return(out);
}

double i_to_r(arma::vec theta, arma::vec current_state){
  double rec = 1 / theta(2);
  int i = current_state(2);
  double out = rec * i;
  return(out);
}

double r_to_s(arma::vec theta, arma::vec current_state){
  double gamma = 1 / theta(3);
  int r = current_state(3);
  double out = gamma * r;
  return(out);
}

double birth_s(arma::vec theta, arma::vec current_state){
  double mu = 1 / theta(4);
  int n = sum(current_state);
  double out = mu * n;
  return(out);
}

double death_s(arma::vec theta, arma::vec current_state){
  double mu = 1 / theta(4);
  int s = current_state(0);
  double out = mu * s;
  return(out);
}

double death_e(arma::vec theta, arma::vec current_state){
  double mu = 1 / theta(4);
  int e = current_state(1);
  double out = mu * e;
  return(out);
}

double death_i(arma::vec theta, arma::vec current_state){
  double mu = 1 / theta(4);
  int i = current_state(2);
  double out = mu * i;
  return(out);
}

double death_r(arma::vec theta, arma::vec current_state){
  double mu = 1 / theta(4);
  int r = current_state(3);
  double out = mu * r;
  return(out);
}


/*
 * next_rate_dynamic updates rate vector based on dependency graph
 * theta: vector of parameters
 * current_rate: vector of event rates
 * current_state: vector of state variables
 * depend: dependency graph (matrix of booleans)
 * j_event: which event fired last
 */
// [[Rcpp::export]]
arma::vec update_rates_seirs_demography(arma::vec theta, arma::vec current_rate, arma::vec current_state, LogicalMatrix depend, int j_event){
  
  arma::vec updated_rate = current_rate;
  
  LogicalVector j_update = depend.row(j_event); //boolean vector of which rates to update
  
  for(int i=0; i<j_update.length(); i++){ //iterate through j_update to choose which rates to update
    
    if(j_update[i]==true){ 
      if(i==0){ //S to E
        updated_rate(i) = s_to_e(theta,current_state);
      }
      if(i==1){ //E to I
        updated_rate(i) = e_to_i(theta,current_state);
      } 
      if(i==2){ //I to R
        updated_rate(i) = i_to_r(theta,current_state);
      } 
      if(i==3){ //R to S
        updated_rate(i) = r_to_s(theta,current_state);
      } 
      if(i==4){ //birth S
        updated_rate(i) = birth_s(theta,current_state);
      } 
      if(i==5){ //death S
        updated_rate(i) = death_s(theta,current_state);
      } 
      if(i==6){ //death E
        updated_rate(i) = death_e(theta,current_state);
      } 
      if(i==7){ //death I
        updated_rate(i) = death_i(theta,current_state);
      } 
      if(i==8){ //death R
        updated_rate(i) = death_r(theta,current_state);
      }
    } //end if
    
  } //end loop
  
  return(updated_rate);
}


/*
 * gillespie_next is a function to simulate a realization of a continuous time Markov process via the next reaction method
 * theta: vector of parameters
 * init_state: vector of state variables
 * trans: matrix describing transitions
 * depend: dependency graph (matrix of booleans)
 * t_end: integer value of simulation length
 * rates: function to evaluate rates given current state and theta
 */
// [[Rcpp::export]]
arma::mat gillespie_next(arma::vec theta, arma::vec init_state, arma::mat trans, LogicalMatrix depend, int t_end, bool info = false){
  
  int n_event = trans.n_rows;
  int j_event;
  
  Progress p(0, false); //progress to check for user interrupt
  
  //initialize trace of Monte Carlo simulation
  arma::mat trace = arma::zeros(1,init_state.n_elem);
  trace.row(0) = init_state.t();
  arma::vec current_state;
  arma::vec current_rates;
  
  //main simulation loop
  double time = 0.0;
  int i = 1;
  
  while(time <= t_end){
    
    //check for user abort
    if(Progress::check_abort()){
      Rcout << "User abort detected at time: " << time << ", i: " << i << std::endl;
      return(arma::zeros(2,2));
    }
    
    if(info){ //print simulation details
      Rcout << "time: " << time << ", i: " << i << std::endl;  
    }
    
    current_state = trace.row(i-1).t(); //extract state at beginning of time jump
    
    if(i == 1){
      current_rates = rate_wrapper(theta,current_state,seirs_demography_rate); //calculate initial rates   
    } else {
      current_rates = update_rates_seirs_demography(theta,current_rates,current_state,depend,j_event); //update rates based on dependency graph
    }
    
    arma::vec tau(n_event); //calculate vector of times to next event
    for(int j=0; j<n_event; j++){
      tau(j) = (-1 / current_rates(j)) * log(R::runif(0,1));
    }

    j_event = tau.index_min(); //fire j_event
    current_state = current_state + trans.row(j_event).t(); //update the current state
    trace.insert_rows(i,current_state.t());
    
    time = time + tau(j_event); //update time
    i++; //update iterator
  }
  
  return(trace);
}


/***R
#SEIRS with demography
trans_seirs <- matrix(c(-1,1,0,0, #infection from S to E
                      0,-1,1,0, #progression from E to I
                      0,0,-1,1, #recovery from I to R
                      1,0,0,-1, #loss of immunity from R to S
                      1,0,0,0, #birth to S
                      -1,0,0,0, #death from S
                      0,-1,0,0, #death from E
                      0,0,-1,0, #death from I
                      0,0,0,-1), #death from R
                      nrow=9,ncol=4,byrow=TRUE)
  
#specify theta (R0, latent duration, infectious duration, immune duration, lifespan)
theta_seirs <- c(5,2,3,365,365*65)

#specify initial state vector
init_state_seirs <- c(1e3,1,0,0)
  
#specify dependency graph
depend_seirs <- matrix(c(T,T,T,T,T,T,T,T,T, #S to E
                       T,T,F,F,F,F,T,F,F, #E to I
                       F,T,T,F,F,F,F,T,F, #I to R
                       F,F,T,T,F,F,F,F,T, #R to S
                       F,F,F,F,F,T,T,T,T, #birth to S
                       T,F,F,T,F,F,F,F,F, #death from S
                       T,T,F,F,F,F,F,F,F, #death from E
                       F,T,T,F,F,F,F,F,F, #death from I
                       F,F,T,T,F,F,F,F,F), #death from R
                       nrow=9,ncol=9,byrow=TRUE)

sim_nextRxn <- gillespie_next(theta_seirs,init_state_seirs,trans_seirs,depend_seirs,20,TRUE)
*/


/***R
#plot output
par(mfrow=c(1,3))
matplot(sim_firstRxn,type="l",ylab="",main="First Reaction Method")
matplot(sim_direct,type="l",ylab="",main="Direct Method")
matplot(sim_nextRxn,type="l",ylab="",main="Next Reaction Method")
par(mfrow=c(1,1))

#plot dependency graph
library(igraph)
depend_seirs_graph <- graph_from_adjacency_matrix(depend_seirs*1)
V(depend_seirs_graph)$name <- c("S->E","E->I","I->R","R->S","->S","S->","E->","I->","R->")
plot(depend_seirs_graph)

#benchmarking flavors of Gillespie SSAs
library(microbenchmark)
microbenchmark(
  gillespie_first(theta_seirs,init_state_seirs,trans_seirs,100,FALSE),
  gillespie_direct(theta_seirs,init_state_seirs,trans_seirs,100,FALSE),
  gillespie_next(theta_seirs,init_state_seirs,trans_seirs,depend_seirs,100,FALSE),
times=500)
*/
