#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


/*
 * C++ samples and snippets from coding simple SMC/Particle Filter algorithm
 */

//function to create null vector and fill it manually
// [[Rcpp::export]]
IntegerVector fillVectorMan(){
  IntegerVector out(Dimension(2));
  out(0) = 0;
  out(1) = 1;
  return(out);
}
/*** R
message("Running fillVectorMan")
fillVectorMan()
*/


//testing initial system for particle filter
// [[Rcpp::export]]
List initPF(NumericMatrix data, NumericVector init_state, int n_particles){
  //initialize system
  int n_iter = data.nrow(); //number of iterations for main particle filter
  NumericVector time_points = data(_,0); //extract time points to run model over
  double loglike = 0.0;
  NumericVector particle_current_state(Dimension(1,init_state.length(),n_particles));
  NumericVector particle_traj(Dimension(n_iter,init_state.length(),n_particles));
  double init_weight = 1 / Rcpp::as<double>(wrap(n_particles));
  NumericVector particle_weight = NumericVector(n_particles,init_weight);
  return(List::create(Named("n_iter")=n_iter,
                      Named("time_points")=time_points,
                      Named("loglike")=loglike,
                      Named("particle_current_state")=particle_current_state,
                      Named("particle_traj")=particle_traj,
                      Named("particle_weight")=particle_weight));
}
/***R
message("Running initPF")
dataMat <- matrix(c(1,2,3,4,5,6,7,8,5,4,2,6,7,5,3,5),ncol=2,byrow=FALSE)
initPF(data=dataMat,init_state=c(3,4,2),n_particles=10)
*/


//testing setting up data and inner loop time iterator
// [[Rcpp::export]]
List innerTimes(NumericMatrix data, int i){
  
  NumericVector time_points = data(_,0); //extract time points to run model over
  
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
  
  return(List::create(Named("times")=times,Named("data_point")=data_point));
}
/***R
message("Running innerTimes")
dataMat <- matrix(c(1,2,3,4,5,6,7,8,5,4,2,6,7,5,3,5),ncol=2,byrow=FALSE)
innerTimes(dataMat,1)
innerTimes(dataMat,2)
*/


//function to return indices of NA values
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
/***R
message("Running whichNA")
whichNA(c(NA,1,4,2,5,NA,4,2,NA,5))
whichNA(rep(10,10))
*/


//test subsetting a NumericVector by an IntegerVector (subset by index)
// [[Rcpp::export]]
NumericVector subsetNumVec(NumericVector input, IntegerVector index){
  NumericVector out(Dimension(index.size()));
  for(int i=0; i<index.size(); i++){
    out(i) = input(index(i));
  }
  return(out);
}
/***R
message("Running subsetNumVec")
subsetNumVec(1:10,c(0,5))
subsetNumVec(1:10,integer(0))
*/


// Replace selected indicies of a NumericVector input with zeros
// [[Rcpp::export]]
NumericVector setZero(NumericVector input, IntegerVector index){
  NumericVector out = input;
  for(int i=0; i<index.size(); i++){
    out(index(i)) = 0;
  }
  return(out);
}
/***R
message("Running setZero")
setZero(1:10,c(0,5))
setZero(1:10,integer(0))
*/


// Replace selected indicies of a NumericVector input with zeros
// [[Rcpp::export]]
NumericVector setZero1(NumericVector input, IntegerVector index){
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
/***R
message("Running setZero1")
setZero1(1:10,c(0,5))
setZero1(1:10,integer(0))
  */


// Test returning SEXP (R) stuff from within C++
// [[Rcpp::export]]
List returnSEXP(){
  return(List::create(Named("margLogLike")=R_NegInf,Named("trajectory")=R_NilValue,Named("traj_weight")=R_NilValue));
}
/***R
message("Running returnSEXP")
returnSEXP()
*/


//testing RcppArmadillo sample
// [[Rcpp::export]]
NumericVector sampleTest(NumericVector x, int size, bool replace, NumericVector probs){
  return(RcppArmadillo::sample(x,size,replace,probs));
}
/***R
message("Running sampleTest")
sampleTest(1:10,10,TRUE,rep(1/10,10))
*/


// testing what the hell Rcpp sugar functions return
// [[Rcpp::export]]
List rcppSugarTest(){
  IntegerVector seqOut = seq_len(10)-1;
  return(List::create(Named("seq")=seqOut));
}
/***R
message("Running seqOut")
rcppSugarTest()
*/