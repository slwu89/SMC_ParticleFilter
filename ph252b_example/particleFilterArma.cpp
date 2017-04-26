#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

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


// subsetCube will subset/sample an Armadillo cube by an index of same outer dimension
// [[Rcpp::export]]
arma::cube subsetCube(arma::cube data, IntegerVector index){
  if(data.n_slices != index.size()){ //informative error message
    Rcout << "subsetCube requires an array and index of the same outer dimension!" << std::endl;
  }
  arma::cube out = arma::zeros(data.n_rows,data.n_cols,data.n_slices);
  for(int i=0; i<data.n_slices; i++){
    out.slice(i) = data.slice(index(i));
  }
  return(out);
}


// subsetMat will subset/sample an armadillo matrix by an index of same row dimension
// [[Rcpp::export]]
arma::mat subsetMat(arma::mat data, IntegerVector index){
  if(data.n_rows != index.size()){
    Rcout << "subsetMat requires a matrix and index of the same row dimension!" << std::endl;
  }
  arma::mat out = data;
  for(int i=0; i<data.n_rows; i++){
    out.row(i) = data.row(index(i));
  }
  return(out);
}


// makeCube will create a cube where the first row of each slice is given as an argument
// [[Rcpp::export]]
arma::cube makeCube(int nrow, int nslice ,arma::rowvec init_state){
  arma::cube out = arma::zeros(nrow,init_state.n_elem,nslice);
  arma::mat data = arma::zeros(nrow,init_state.n_elem);
  data.row(0) = init_state;
  for(int i=0; i<out.n_slices; i++){
    out.slice(i) = data;
  }  
  return(out);
}


/*
 * Bootstrap Particle Filter (SMC)
 * particle_filter returns an R list:
 * loglike is the estimate of the marginal log likelihood
 * traj is the filtered particle trajectories
 * traj_weight is the weight of each particle trajectory
 */
// [[Rcpp::export]]
List particleFilterArmadillo(Function model, NumericVector theta, arma::vec init_state, NumericMatrix data, int n_particles, bool progress){
  
  //initialize system
  int n_iter = data.nrow(); //number of iterations for main particle filter
  NumericVector time_points = data(_,0); //extract time points to run model over
  double loglike = 0.0;
  arma::mat particle_current_state = arma::zeros(n_particles,init_state.n_elem);
  particle_current_state.each_row() = init_state.t();
  arma::cube particle_traj = makeCube(n_iter+1,n_particles,init_state.t());
  double init_weight = 1 / Rcpp::as<double>(wrap(n_particles));
  NumericVector particle_weight = NumericVector(n_particles,init_weight);
  IntegerVector n_particles_vec = seq_len(n_particles)-1; //0:n_particles needed for armadillo bootstrap resample
  
  //progress bar
  Progress p(n_iter,progress);
  
  //run the particle filter
  for(int i=0; i<n_iter; i++){

    //check for user abort
    if(Progress::check_abort()){
      Rcout << "User abort detected at step " << i << std::endl;
      return(List::create(Named("margLogLike")=R_NegInf,Named("trajectory")=R_NilValue,Named("traj_weight")=R_NilValue));
    }
    
    //run inner loop only over single timestep indexed by times
    NumericVector times(Dimension(2));
    int begin;
    int end;
    if(i == 0){
      begin = 0;
    } else {
      begin = time_points(i-1);
    }
    end = time_points(i);
    times(0) = begin;
    times(1) = end;
    double data_point = data(i,1); //current observed data point
    
    //resample particles
    IntegerVector index_resampled;
    IntegerVector na_weight = whichNA(particle_weight); //find particles with NA weight
    particle_weight = setZeroIndex(particle_weight,na_weight); //replace NA weight particles with 0
    if(!all(particle_weight==0).is_true()){
      index_resampled = RcppArmadillo::sample(n_particles_vec,n_particles,true,particle_weight);
    } else {
      Rcout << "All particles depleted at step " << i << " of SMC. Return margLogLike = -Inf for theta: " << theta << std::endl;
      return(List::create(Named("margLogLike")=R_NegInf,Named("trajectory")=R_NilValue,Named("traj_weight")=R_NilValue));
    }
    
    //update particles after resampling
    particle_traj = subsetCube(particle_traj,index_resampled);
    particle_current_state = subsetMat(particle_current_state,index_resampled);
    
    //propagate particles over particle_current_state
    arma::mat propagate_model_point = arma::zeros(n_particles,init_state.n_elem); //states for each of the j particles after being propagated to next time point
    NumericVector propagate_weight = NumericVector(n_particles); //weights for each of the j particles after being propagated to next time point
    
    for(int j=0; j<n_particles; j++){ //run model for each particle
      
      //run model over current time step
      arma::rowvec current_state_j = particle_current_state.row(j);
      Rcpp::DataFrame particle_j_out = model(theta,current_state_j,times);
      
      //extract trajectories at latest time step
      IntegerVector j_times = particle_j_out["time"];
      NumericVector j_S = particle_j_out["S"];
      NumericVector j_E1 = particle_j_out["E1"];
      NumericVector j_E2 = particle_j_out["E2"];
      NumericVector j_I1 = particle_j_out["I1"];
      NumericVector j_I2 = particle_j_out["I2"];
      NumericVector j_R1 = particle_j_out["R1"];
      NumericVector j_R2 = particle_j_out["R2"];
      NumericVector j_R3 = particle_j_out["R3"];
      NumericVector j_R4 = particle_j_out["R4"];
      NumericVector j_L = particle_j_out["L"];
      propagate_model_point(j,0) = j_S(j_S.size()-1);
      propagate_model_point(j,1) = j_E1(j_E1.size()-1);
      propagate_model_point(j,2) = j_E2(j_E2.size()-1);
      propagate_model_point(j,3) = j_I1(j_I1.size()-1);
      propagate_model_point(j,4) = j_I2(j_I2.size()-1);
      propagate_model_point(j,5) = j_R1(j_R1.size()-1);
      propagate_model_point(j,6) = j_R2(j_R2.size()-1);
      propagate_model_point(j,7) = j_R3(j_R3.size()-1);
      propagate_model_point(j,8) = j_R4(j_R4.size()-1);
      propagate_model_point(j,9) = j_L(j_L.size()-1);
      
      //compute particle weight
      double f = 1.0/theta["lat_dur"];
      double rho = 0.72;
      double lambda = f*rho*(propagate_model_point(j,1) + propagate_model_point(j,2));
      double j_weight = R::dpois(data_point,lambda,false); //incidence is modeled as a Poisson process
      
      //collect parallel trajectories of each particle
      propagate_weight(j) = j_weight;
      particle_traj(i+1,0,j) = j_S(j_S.size()-1);
      particle_traj(i+1,1,j) = j_E1(j_E1.size()-1);
      particle_traj(i+1,2,j) = j_E2(j_E2.size()-1);
      particle_traj(i+1,3,j) = j_I1(j_I1.size()-1);
      particle_traj(i+1,4,j) = j_I2(j_I2.size()-1);
      particle_traj(i+1,5,j) = j_R1(j_R1.size()-1);
      particle_traj(i+1,6,j) = j_R2(j_R2.size()-1);
      particle_traj(i+1,7,j) = j_R3(j_R3.size()-1);
      particle_traj(i+1,8,j) = j_R4(j_R4.size()-1);
      particle_traj(i+1,9,j) = j_L(j_L.size()-1);
    }
    
    //collect parallel trajectories of each particle
    particle_current_state = propagate_model_point;
    particle_weight = propagate_weight;
    
    //update marginal log-likelihood of theta
    loglike += log(mean(particle_weight));
    
    //advance progress bar
    p.increment();
    
  }
  
  //return estimated marginal log likelihood, particle trajectory, and particle weights
  return(List::create(Named("loglike")=loglike,Named("trajectory")=particle_traj,Named("traj_weight")=particle_weight/sum(particle_weight)));
}