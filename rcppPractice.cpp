#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(BH)]]
#include <boost/function.hpp>
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


// testing returning -Inf and NULL
// [[Rcpp::export]]
List rcppNA(){
  return(List::create(Named("NegInf")=R_NegInf,
                      Named("NULL")=R_NilValue,
                      Named("NaN")=R_NaN,
                      Named("NA")=NumericVector::get_na()));
}
/***R
message("Running rcppNA")
rcppNA()
*/


// how to subset outer dimension of arma cube by IntegerVector
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
/***R
message("Running subsetCube")
subsetCube(array(rep(1:10,each=9),c(3,3,10)),c(0,0,0,5,5,6,6,7,8,9))
*/


// how to subset rows of NumericMatrix by IntegerVector
// [[Rcpp::export]]
NumericMatrix subsetMat(NumericMatrix data, IntegerVector index){
  if(data.nrow() != index.size()){
    Rcout << "subsetMat requires a matrix and index of the same row dimension!" << std::endl;
  }
  NumericMatrix out = NumericMatrix(Dimension(data.nrow(),data.ncol()));
  for(int i=0; i<data.nrow(); i++){
    out(i,_) = data(index(i),_);
  }
  return(out);
}
/***R
message("Running subsetMat")
subsetMat(matrix(rep(1:10,each=3),ncol=3,byrow=T),c(0,0,0,5,5,6,6,7,8,9))
*/


// how to subset rows of armadillo mat by IntegerVector
// [[Rcpp::export]]
arma::mat subsetArmaMat(arma::mat data, IntegerVector index){
  if(data.n_rows != index.size()){
    Rcout << "subsetMat requires a matrix and index of the same row dimension!" << std::endl;
  }
  arma::mat out = data;
  for(int i=0; i<data.n_rows; i++){
    out.row(i) = data.row(index(i));
  }
  return(out);
}
/***R
message("Running subsetMat")
subsetArmaMat(matrix(rep(1:10,each=3),ncol=3,byrow=T),c(0,0,0,5,5,6,6,7,8,9))
*/

// evaluating dpois through Rcpp
// [[Rcpp::export]]
double poisPDF(double x, double lambda){
  double out = R::dpois(x,lambda,true);
  return(out);
}
/***R
message("Running poisPDF")
poisPDF(5,6)
*/


// armadillo cube dimensions
// [[Rcpp::export]]
arma::cube cubeMake(int n_iter, NumericVector init_state, int n_particles){
  return(arma::zeros(n_iter,init_state.length(),n_particles));
}
/***R
message("Running cubeMake")
cubeMake(10,c(1,2,3),3)
*/


// armadillo subcube testing
// [[Rcpp::export]]
arma::cube subCube(arma::cube cubeDat, int first_row, int first_col, int first_slice, int last_row, int last_col, int last_slice){
  return(cubeDat.subcube(first_row, first_col, first_slice, last_row, last_col, last_slice));
}
/***R
message("Running subCube")
subCube(array(1:300,c(10,3,10)),0,0,0,2,2,9)
*/


// making empty cube except for first row in each sub mat
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
/***R
message("Running makeCube")
makeCube(10,10,c(1,2,3))
*/


// dealing with data frame output from ssa.adaptivetau in C++
// [[Rcpp::export]]
List dataframeAdaptTau(Function sim_algorithm,arma::vec theta, arma::vec init_state, IntegerVector times){
  Rcpp::DataFrame sim_out = sim_algorithm(theta,init_state,times);
  IntegerVector sim_times = sim_out["time"];
  NumericMatrix sim_trace = Dimension(times.size(),theta.n_elem+1);
  NumericVector outS = sim_out["S"];
  NumericVector outI = sim_out["I"];
  NumericVector outR = sim_out["R"];
  sim_trace(_,0) = outS;
  sim_trace(_,1) = outI;
  sim_trace(_,2) = outR;
  return(List::create(Named("sim_times")=sim_times,Named("sim_trace")=sim_trace));
}
/***R
message("Running dataframeAdaptTau")
library(adaptivetau)
#stochastic SIR simulation using ssa.adaptivetau
SIR_stoch <- function(theta,init_state,times){

#fix names to deal with armadillo input vectors
names(theta) <- c("R0","inf_dur")
names(init_state) <- c("S","I","R")
  
#create transition matrix for Markov process
  trans <- list(
      c(S=-1,I=1), #infection
        c(I=-1,R=1) #recovery
  )
  
#create transition rates for Markov process at time t
  rates <- function(current_state,theta,t){
    
#define parameters
    beta <- theta[["R0"]] / theta[["inf_dur"]]
    gamma <- 1/theta[["inf_dur"]]
    
#define states
    S <- current_state[["S"]]
    I <- current_state[["I"]]
    R <- current_state[["R"]]
    N <- S + I + R
    
#return rates
    rates <- c(
        beta * S * I/N, #infection rate
        gamma * I #recovery rate
    )
    
    return(rates)
  }
  
#run the simulation
  mod_out <- as.data.frame(ssa.adaptivetau(init.values=init_state,transitions=trans,rateFunc=rates,params=theta,tf=diff(range(times))))
  
#need to interpolate continuous time Markov process output to discrete time points
  mod_out$time <- mod_out$time + min(times)
  
#interpolate output to integer values
  int_out <- apply(mod_out[,-1],2,function(i){
    approx(x=mod_out[,1],y=i,xout=times,method="constant")$y
  })
  
  mod_traj <- cbind(time=times,int_out)
  return(as.data.frame(mod_traj))
}

# init_theta_sir <- c(R0=1.5,inf_dur=5)
# init_state_sir <- c(S=999,I=1,R=0)
# times <- 1:5
# 
# dataframeAdaptTau(SIR_stoch,init_theta_sir,init_state_sir,times)
dataframeAdaptTau(SIR_stoch,c(1.5,5),c(999,1,0),1:5)
*/


/*
 * How to pass functions as arguments
 */
typedef boost::function<double(double,double)> funcArg;

double funcArg_function(double x, double y){
  return(x*y);
}

double funcArg_wrapper(double x, double y, funcArg somefunc){
  return(somefunc(x,y));
}

// [[Rcpp::export]]
void funcArg_run(double x, double y){
  Rcout << "testing " << funcArg_wrapper(x,y,funcArg_function) << std::endl;
}

/***R
funcArg_run(4.23,2.3)
*/


// testing inserting rows
// [[Rcpp::export]]
arma::mat insertRow(){
  arma::mat A = arma::randu(5,4);
  arma::mat B = arma::ones(2,4);
  
  A.insert_rows(5,B);
  return(A);
}
/***R
insertRow()
#looks like you insert it into the empty space.
*/
