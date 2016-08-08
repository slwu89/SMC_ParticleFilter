#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

/*
 * This file includes functions to run a standard bootstrap particle filter in C++ through Rcpp
 * 
 */


// whichNA will return NA indicies
// [[Rcpp::export]]
IntegerVector whichNA(NumericVector input){
  IntegerVector out;
  for(int i=0; i<input.size(); i++){
    if(R_IsNA(input(i))){
      out.push_back(i);
    }
  }
  return(out);
}


// setZeroIndex will replace indexed elements of input vector with zero
// [[Rcpp::export]]
NumericVector setZeroIndex(NumericVector input, IntegerVector index){
  NumericVector out = input;
  if(index.isNULL()){
    return(out);
  } else {
    for(int i=0; i<index.size(); i++){
      out(index(i)) = 0;
    }
    return(out); 
  }
}


// [[Rcpp::export]]
/*
 * 
 */
List particle_filter(Function model, NumericVector theta, NumericVector init_state, NumericMatrix data, int n_particles, bool info){
  
  //initialize system
  int n_iter = data.nrow(); //number of iterations for main particle filter
  NumericVector time_points = data(_,0); //extract time points to run model over
  double loglike = 0.0;
  arma::cube particle_current_state = arma::zeros(1,init_state.length(),n_particles);
  arma::cube particle_traj = arma::zeros(n_iter,init_state.length(),n_particles);
  //NumericVector particle_current_state(Dimension(1,init_state.length(),n_particles));
  //NumericVector particle_traj(Dimension(n_iter,init_state.length(),n_particles));
  double init_weight = 1 / Rcpp::as<double>(wrap(n_particles));
  NumericVector particle_weight = NumericVector(n_particles,init_weight);
  IntegerVector n_particles_vec = seq_len(n_particles)-1;
  
  //progress bar
  Progress p(n_iter,info);
  
  //run the particle filter REMEMBER TO SUBSET WITH i-1!!!
  for(int i=1; i<=n_iter; i++){
    
    //run inner loop only over single timestep indexed by times
    NumericVector times(Dimension(2));
    int begin;
    int end;
    if(i == 1){
      begin = 0;
    } else {
      begin = time_points(i-2);
    }
    end = time_points(i-1);
    times(0) = begin;
    times(1) = end;
    NumericVector data_point = data(i-1,_); //vector of the current observed data
    
    //resample particles
    IntegerVector na_weight = whichNA(particle_weight); //find particles with NA weight
    particle_weight = setZeroIndex(particle_weight,na_weight); //replace NA weight particles with 0
    if(!all(particle_weight==0).is_true()){
      IntegerVector index_resampled = RcppArmadillo::sample(n_particles_vec,n_particles,true, particle_weight);
    } else {
      Rcout << "All particles depleted at step " << i << " of SMC. Return margLogLike = -Inf for theta: " << theta << std::endl;
      return(List::create(Named("margLogLike")=R_NegInf,Named("trajectory")=R_NilValue,Named("traj_weight")=R_NilValue));
    }
    
    //update particles after resampling
    particle_traj ;
    particle_current_state ;

    
    
  }
  
  //return estimated marginal log likelihood, particle trajectory, and particle weights
  return(List::create(Named("loglike")=loglike,Named("trajectory")=particle_traj,Named("traj_weight")=particle_weight/sum(particle_weight)));
}


/*** R

*/
