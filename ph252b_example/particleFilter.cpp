#include <RcppGSL.h>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppGSL)]]

/*
 * This file includes functions to run a standard bootstrap particle filter in C++ through Rcpp
 * The function particle_filter will run a bootstrap resample particle filter for count data (infections),
 * and it is intended to be used for stochastic models such as Gillespie algorithm style models,
 * implemented by the adaptivetau package. It can be modified for multiple data sources and more complex
 * likelihood calculations easily, but is currently set up to estimated the log-likelihood of a stochastic
 * SIR model with incidence modeled as a Poisson process.
 * 
 * You will need to install Rcpp, RcppArmadillo, and RcppProgress packages and have a working C++ compiler to use it
 * 
 * The function takes the following arguments:
 * model: this is an R function (usually a stochastic model; example here uses ssa.adaptivetau, a modification of the Gillespie algorithm)
 *        that returns a data frame of the state values indexed by time
 * theta: the vector of parameters that defines the model
 * init_state: the initial state values of the model
 * data: the data that is used for likelihood calculations
 * n_particles: the number of particles used to evaluate the estimated log-likelihood
 * info: boolean to print detailed info or not
 */


// Particle Filter (SMC) utility functions

// whichNA will return NA indicies of a NumericVector
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


// // subsetCube will subset/sample an Armadillo cube by an index of same outer dimension
// // [[Rcpp::export]]
// arma::cube subsetCube(arma::cube data, IntegerVector index){
//   if(data.n_slices != index.size()){ //informative error message
//     Rcout << "subsetCube requires an array and index of the same outer dimension!" << std::endl;
//   }
//   arma::cube out = arma::zeros(data.n_rows,data.n_cols,data.n_slices);
//   for(int i=0; i<data.n_slices; i++){
//     out.slice(i) = data.slice(index(i));
//   }
//   return(out);
// }
// 
// 
// // subsetMat will subset/sample an armadillo matrix by an index of same row dimension
// // [[Rcpp::export]]
// arma::mat subsetMat(arma::mat data, IntegerVector index){
//   if(data.n_rows != index.size()){
//     Rcout << "subsetMat requires a matrix and index of the same row dimension!" << std::endl;
//   }
//   arma::mat out = data;
//   for(int i=0; i<data.n_rows; i++){
//     out.row(i) = data.row(index(i));
//   }
//   return(out);
// }
// 
// 
// // makeCube will create a cube where the first row of each slice is given as an argument
// // [[Rcpp::export]]
// arma::cube makeCube(int nrow, int nslice ,arma::rowvec init_state){
//   arma::cube out = arma::zeros(nrow,init_state.n_elem,nslice);
//   arma::mat data = arma::zeros(nrow,init_state.n_elem);
//   data.row(0) = init_state;
//   for(int i=0; i<out.n_slices; i++){
//     out.slice(i) = data;
//   }  
//   return(out);
// }


/*
 * Bootstrap Particle Filter (SMC)
 * particle_filter returns an R list:
 * loglike is the estimate of the marginal log likelihood
 * traj is the filtered particle trajectories
 * traj_weight is the weight of each particle trajectory
 */

// imitate R's sample
IntegerVector sampleGSL(const gsl_rng *r, int x, NumericVector probs){
  
  size_t K = int(x);
  
  gsl_ran_discrete_t * sampler;
  sampler = gsl_ran_discrete_preproc(K, probs.begin());
  
  size_t ii;
  IntegerVector out(x);
  
  for(int i=0; i<x; i++){
    ii = gsl_ran_discrete(r,sampler);
    out[i] = ii;
  }
  
  gsl_ran_discrete_free(sampler);
  
  return(out);
}

typedef std::vector<NumericVector> particleCurrentStateDefn;
typedef std::vector<particleCurrentStateDefn> particleTrajDefn;

// [[Rcpp::export]]
Rcpp::List particleFilter(Function model, NumericVector theta, NumericVector initState, NumericMatrix data, int nParticles, bool progress, int seed){
  
  // initialize random number generator
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  
  // initialize particle filter
  int nIter = data.nrow(); // number of times to bootstrap weighted particle trajectories
  NumericVector timePoints = data.column(0); // first column must be integer times
  double logLike = 0.0; // marginal log-likelihood of theta
  
  particleCurrentStateDefn particleCurrentState(nParticles, initState); // current state of each particle
  particleTrajDefn particleTrajectory(nIter+1, particleCurrentState); // all filtered trajectories
  
  double initWeight = 1.0 / double(nParticles); // uniform initial weights of particles
  NumericVector particleWeight = NumericVector(nParticles, initWeight);
  
  // interruptable progress bar
  Progress prog(nIter, progress);
  
  // run particle filter
  for(int i = 0; i < nIter; i++){
    
    //check for user abort
    if(Progress::check_abort()){
      Rcout << "User abort detected at step " << i << std::endl;
      return(List::create(Named("logLike")=R_NegInf,Named("trajectory")=R_NilValue,Named("weights")=R_NilValue));
    }
    
    // time step to propagate particles over
    IntegerVector times(2);
    if(i == 0){
      times[0] = 0;
      times[1] = timePoints[i];
    } else {
      times[0] = timePoints[i-1];
      times[1] = timePoints[i];
    }
    double dataPoint = data(i,1); // observed data point
    
    // bootstrap resample weighted particles
    IntegerVector naWeight = whichNA(particleWeight); // find particles with NA weight
    particleWeight = setZeroIndex(particleWeight, naWeight); // replace NA weight particles with 0
    IntegerVector indexResample;
    if(!all(particleWeight==0).is_true()){
      indexResample = sampleGSL(r, nParticles, particleWeight);
    } else {
      Rcout << "All particles depleted at step " << i << " of SMC. Return margLogLike = -Inf for theta: " << theta << std::endl;
      return(List::create(Named("logLike")=R_NegInf,Named("trajectory")=R_NilValue,Named("weights")=R_NilValue));
    }
    
    // update filtered trajectories of particles after bootstrap resampling
    for(int i=0; i<particleTrajectory.size(); i++){
      for(int ii=0; ii<indexResample.size(); ii++){
        particleTrajectory[i][ii] = particleTrajectory[i][indexResample[ii]];
      }
    }
    
    // update current state of particles after bootstrap resampling
    for(int i=0; i<particleCurrentState.size(); i++){
      particleCurrentState[i] = particleCurrentState[indexResample[i]];
    }

    // propagate particles over particleCurrentState
    particleCurrentStateDefn particlesPropagateState(nParticles); // states for each of the j particles after being propagated to next time point
    NumericVector propagateWeight(nParticles); // weights for each of the j particles after being propagated to next time point

    // propagate particles
    for(int i=0; i<nParticles; i++){
      
    }
    
  }
  
  return(List::create(
      _["particleCurrentState"] = particleCurrentState,
      _["particleTrajectory"] = particleTrajectory
      // _["indexResample"] = indexResample
  ));
}

/***R
out = particleFilter(model = function(){1},theta = c(1),initState = c(1,10,100),data = matrix(rnorm(9)),nParticles = 5,progress = FALSE,seed = 42)

*/