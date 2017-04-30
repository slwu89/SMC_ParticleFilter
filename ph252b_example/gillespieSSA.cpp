#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp11)]]


/*
* This file gives functions to run exact discrete event (integer valued) simulations of stochastic processes.
* Gillespie's Direct Method
* 
* Sean Wu
* April 20, 2017
* 
* You will need to install Rcpp, RcppGSL, and have a working C++ compiler to use it.
*/

/*
 * hazardSEIRL: return instantaneous hazard rates for each event
 * theta: vector of parameters
 * currentState: vector of state variables
 */
// [[Rcpp::export]]
std::vector<double> hazardSEIRL(NumericVector theta, NumericVector currentState){
  
  //define parameters
  double beta = theta["R0"] / theta["inf_dur"];
  double f = 1.0 / theta["lat_dur"];
  double rec = 1.0 / theta["inf_dur"];
  double tau = 1.0 / theta["lat_dur"];
  double alpha = theta["alpha"];
  
  //define states
  double S = currentState["S"];
  double E1 = currentState["E1"];
  double E2 = currentState["E2"];
  double I1 = currentState["I1"];
  double I2 = currentState["I2"];
  double R1 = currentState["R1"];
  double R2 = currentState["R2"];
  double R3 = currentState["R3"];
  double R4 = currentState["R4"];
  double L = currentState["L"];
  double I = I1 + I2;
  double N = S + E1 + E2 + I1 + I2 + R1 + R2 + R3 + R4 + L;
  
  //return rates
  std::vector<double> rates(11);
  rates[0] = beta * S * I/N; // susceptible to exposed 1
  rates[1] = 2.0 * f * E1; // exposed 1 to exposed 2
  rates[2] = 2.0 * f * E2; // exposed 2 to infectious 1
  rates[3] = 2.0 * rec * I1; // infectious 1 to infectious 2
  rates[4] = 2.0 * rec * I2; // infectious 2 to recovered 1
  rates[5] = 4.0 * tau * R1; // recovered 1 to recovered 2
  rates[6] = 4.0 * tau * R2; // recovered 2 to recovered 3
  rates[7] = 4.0 * tau * R3; // recovered 3 to recovered 4
  rates[8] = 4.0 * tau * R4; // out of recovered 4
  rates[9] = 4.0 * (1-alpha) * tau * R4; // recovered 4 to susceptible
  rates[10] = 4.0 * alpha * tau * R4; // recovered 4 to immune

  return(rates);
}

/*
 * hazardSIR: return instantaneous hazard rates for each event
 * theta: vector of parameters
 * currentState: vector of state variables
 */
// [[Rcpp::export]]
std::vector<double> hazardSIR(NumericVector theta, NumericVector currentState){
  
  //define parameters
  double beta = theta["R0"] / theta["infDur"];
  double gamma = 1.0 / theta["infDur"];
  
  //define states
  int s = currentState["S"];
  int i = currentState["I"];
  int r = currentState["R"];
  int n = s + i + r;
  
  //return rates
  std::vector<double> rates(2);
  rates.at(0) = beta * s * i/n; //S to I
  rates.at(1) = gamma * i; //I to R
  return(rates);
}


/*
 * gillespieDirect is a function to simulate a realization of a continuous time Markov process via the direct reaction method
 * theta: vector of parameters
 * initState: vector of state variables
 * trans: matrix describing transitions
 * tEnd: integer value of simulation length
 */
// [[Rcpp::export]]
Rcpp::List gillespieDirect(NumericVector theta, NumericVector initState, NumericMatrix trans, int tEnd, int seed, bool info = false){
  
  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  
  // initialize trace of Monte Carlo simulation
  std::vector<NumericVector> trace;
  trace.push_back(clone(initState));
  std::vector<double> times;

  // main simulation loop
  double time = 0.0;
  int i = 0;
  
  times.push_back(time);
  
  while(time <= tEnd){
    
    // print info
    if(info){
      Rcout << "time: " << time << ", i: " << i << std::endl;
    }
    
    NumericVector currentState = trace.at(i); // extract state at beginning of time jump
    std::vector<double> currentRates = hazardSEIRL(theta, currentState); // calculate event rates
    // std::vector<double> currentRates = hazardSIR(theta, currentState); // calculate event rates
    
    double w0 = std::accumulate(currentRates.begin(), currentRates.end(), 0.0); // sum of rate (propensity) functions
    double tau = (1.0/w0) * log(1.0/gsl_ran_flat(r,0.0,1.0)); // calculate time to next reaction
    
    double rNum = gsl_ran_flat(r,0.0,1.0);
    double w0rxn = rNum * w0; // instantaneous event probability
    
    // decide which event j occured
    int j = 1;
    while(std::accumulate(currentRates.begin(), currentRates.begin()+j, 0.0) < w0rxn){
      j++;
    }
    
    currentState = currentState + trans.row(j-1);
    trace.push_back(clone(currentState));
    
    time = time + tau; // update time
    times.push_back(time);
    i++; // update iterator
  }
  
  return(List::create(
      _["trace"] = trace,
      _["times"] = times
  ));
}