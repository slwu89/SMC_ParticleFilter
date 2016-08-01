#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
using namespace Rcpp;

/*
 * This file includes functions to run a standard bootstrap particle filter in C++ through Rcpp
 * 
 */

// [[Rcpp::export]]
List particle_filter(Function model, NumericVector theta, NumericVector init_state, NumericMatrix data, int n_particles, bool info){
  
  //initialize system
  int n_iter = data.nrow(); //number of iterations for main particle filter
  NumericVector time_points = data(_,0); //extract time points to run model over
  double loglike = 0.0;
  NumericVector particle_current_state(Dimension(1,init_state.length(),n_particles));
  NumericVector particle_traj(Dimension(n_iter,init_state.length(),n_particles));
  NumericVector particle_weight = NumericVector(1/n_particles,n_particles);
  
  //progress bar
  Progress p(n_iter,info);
  
  //run the particle filter REMEMBER TO SUBSET WITH i-1!!!
  for(int i=1; i<=n_iter; i++){
    
    NumericVector times(Dimension(2));
    int begin;
    int end;
    if(i == 1){
      begin = 0;
    } else {
      begin = time_points(i-1);
    }
    end = time_points(i);
    
  }
  
  //return estimated marginal log likelihood, particle trajectory, and particle weights
  return(List::create(Named("loglike")=loglike,Named("trajectory")=particle_traj,Named("traj_weight")=particle_weight/sum(particle_weight)));
}


/*** R

*/
