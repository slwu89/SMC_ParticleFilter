#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppProgress)]]
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <algorithm> //for sort and unique
#include <vector>
#include <iterator>
#include <ctime>
#include <progress.hpp>
#include <boost/function.hpp>
using namespace Rcpp;

typedef unsigned int COUNTER;
typedef unsigned int NODE;
typedef std::vector<NODE> NODES; // list of nodes
typedef std::vector<bool> BOOLS;
struct CONTACT {NODE i; NODE j;}; // contact (i,j)
typedef std::vector<CONTACT> CONTACTS; // contacts in a single time-frame
typedef std::vector<CONTACTS> CONTACTS_LIST; // list of contact lists
typedef std::vector<double> CHARACTERISTICS; // characteristics of nodes, e.g., susceptibilities
// Random number generators:
typedef boost::mt19937_64 ENG; // use Mersenne Twister 19937 as PRNG engine
typedef boost::uniform_int<> DIST_INT; // define uniform distribution of integers
typedef boost::uniform_real<> DIST_REAL; // define uniform distribution of reals on [0,1)
typedef boost::exponential_distribution<> DIST_EXP; // define exponential distribution





/*
 * This file includes functions to run stochastic simulation algorithms on static networks.
 */

// [[Rcpp::export]]
arma::mat gillespie_direct_network(){
  return(arma::zeros(3,3));
}


//======================================================================
// Global parameters
//======================================================================
COUNTER N; //number of nodes in network
COUNTER T_data; //length of dataset (time steps)
COUNTER dt; //input time resolution (integer)
CHARACTERISTICS susceptibilities;
CHARACTERISTICS infectivities;
CHARACTERISTICS recoverabilities;
char inputname[200], outputname[200];
// Input/Output streams:
std::ofstream output;
std::ifstream input;
// Initialize PRNG:
ENG eng(9071982);
DIST_INT dist_int;
DIST_REAL dist_rand(0,1);
DIST_EXP dist_exp;

//======================================================================
// Function for importing list of contacts from tij format file:


//======================================================================
// Main:
//======================================================================
int main(int argc, char *argv[])
{
  //-------------------------------------------------------------------------------------
  // Load contact data set:
  //-------------------------------------------------------------------------------------
  // Check if correct number of parameters was passed to program:
  if(argc<3){ std::cout << "Error! Data file not defined.\n"; return 0; }
  else
  {
    if(argc<5){ std::cout << "Error! Epidemic parameters not defined.\n"; return 0; }
    else
    {
      if(argc<6){ std::cout << "Error! Simulation time not specified.\n"; return 0; }
      else
      {
        if(argc<7){ std::cout << "Error! Ensemble size not specified.\n"; return 0; }
        else
        {
          if(argc<8){ std::cout << "Error! Output time-resolution not specified.\n"; return 0;}
        }
      }
    }
  
  COUNTER dt; //input time resolution (integer)}
  
  // Set parameter values as specified:
  char *datafile=argv[1]; //dataset file name
  dt=atoi(argv[2]); //dataset time resolution (integer)
  char *beta_str = argv[3]; // base infection rate
  double beta = atof(beta_str);
  char *mu_str = argv[4]; // base recovery rate
  double mu = atof(mu_str);
  COUNTER betaprec=strlen(beta_str)-2;
  COUNTER muprec=strlen(mu_str)-2;
  COUNTER T_simulation = atoi(argv[5]); //simulation time
  COUNTER ensembleSize = atoi(argv[6]); //ensemble size (number of realizations)
  COUNTER outputTimeResolution = atoi(argv[7]); //output time-resolution
  
  // Open input file and load contact_lists:
  sprintf(inputname,"%s",datafile);
  CONTACTS_LIST contactListList=loadContactListList(inputname);
  // Check if length of contacts_list > 0 and end program if it is not:
  if(contactListList.size()==0){ std::cout << "Error! Dataset empty.\n"; return 0; }
  
  //-------------------------------------------------------------------------------------
  // Define variables:
  //-------------------------------------------------------------------------------------
  NODES infected; //list of infected nodes
  CHARACTERISTICS mus; //list of cumulative sums of their recovery rates
  double Mu; //total recovery rate
  BOOLS isInfected; //list which nodes are susceptible/infected, respectively
  COUNTER I; //number of infected and recovered nodes
  COUNTER SI; //number of susceptible nodes in contact with infectious nodes
  NODES si_s; //list of susceptible nodes in contact with infected nodes
  CHARACTERISTICS betas; //list of cumulative sums of their infection rates
  double Beta; //total infection rate
  double Lambda; //total transition rate
  double xi;
  COUNTER t; //time counter
  COUNTER t_infectionStart; //starting time of infection
  NODE root; //root node of infection
  double tau; //renormalized waiting time until next event
  NODE i,j;
  CONTACTS_LIST::iterator contact_iterator; //iterator over list of contacts
  CONTACTS::iterator contactList_iterator; //iterator over contacts
  double r_transition; //random variable for choosing which transition happens
  CHARACTERISTICS::iterator r_weightedSampling; //iterator for drawing the transition
  COUNTER m; //number of this transition
  COUNTER n; //time counter
  double lambda_draw; //rate of drawn transition
  double lambda_m_new;
  CHARACTERISTICS::iterator mu_iterator; //iterator over list of recovery rates
  NODES::iterator node_iterator; //iterator over list of nodes
  NODES::iterator last; //iterator for use when generating unique list of new infected nodes
  // Containers for output data:
  NODES sumI_t(T_simulation/outputTimeResolution); //list of number of infected nodes in each recorded frame
  NODES hist_I(N+1); //histogram of I values at end of simulations
  // Random number generators:
  boost::variate_generator<ENG,DIST_INT> randint(eng,dist_int); //random integer
  boost::variate_generator<ENG,DIST_REAL> rand(eng,dist_rand); //random float on [0,1[
  boost::variate_generator<ENG,DIST_EXP> randexp(eng,dist_exp); //random exponentially distributed float
  // Assign individual values of population:
  recoverabilities.assign(N,1.);
  susceptibilities.assign(N,1.);
  infectivities.assign(N,1.);
  
  //-------------------------------------------------------------------------------------
  // Simulate:
  //-------------------------------------------------------------------------------------
  std::clock_t clockStart = std::clock();     //timer
  COUNTER stopped=0; //counter of number of simulations that stopped (I=0) during T_simulation
  for(int q=0; q<ensembleSize; q++)
  {
    std::cout << q << "/" << ensembleSize << std::endl; //print realization # to screen
    // Choose at random infectious root node and run SIR process starting from root:
    root=randint(N);
    // Clear parameters:
    infected.clear();
    // Initialize lists of infected nodes and infected node IDs:
    infected.push_back(root);
    I=1;
    mus.assign(1,0.);
    mus.push_back(mu*recoverabilities[root]);
    Mu=mus.back();
    isInfected.assign(N,false);
    isInfected[root]=true;
    // First waiting time:
    tau=randexp(1);
    // random starting time of infection:
    t_infectionStart=randint(T_data);
    // set simulation time to zero:
    t=0;
    
    //--- Loop over list of contact lists: ---
    while(I>0 && t<T_simulation) //loop until either I=0 or t>=T_simulation
    {
      for(contact_iterator=contactListList.begin()+t_infectionStart; contact_iterator!=contactListList.end(); contact_iterator++)
      {
        // Create list of susceptible nodes in contact with infected nodes:
        si_s.clear();
        betas.assign(1,0.);
        for(contactList_iterator=(*contact_iterator).begin(); contactList_iterator!=(*contact_iterator).end(); contactList_iterator++)
        {
          i=(*contactList_iterator).i;
          j=(*contactList_iterator).j;
          if(isInfected[i])
          {
            if(!isInfected[j])
            {
              si_s.push_back(j);
              betas.push_back(beta*susceptibilities[j]*infectivities[i]+betas.back());
            }
          }
          else
          {
            if(isInfected[j])
            {
              si_s.push_back(i);
              betas.push_back(beta*susceptibilities[i]*infectivities[j]+betas.back());
            }
          }
        }
        SI=si_s.size();
        Beta=betas.back();
        Lambda=Beta+Mu;
        
        // Check if transition takes place during time-step:
        if(tau>=Lambda) //no transition takes place
        {
          tau-=Lambda;
        }
        else //at least one transition takes place
        {
          xi=1.;
          // Sampling step:
          while(tau<xi*Lambda) //repeat if next tau is smaller than ~ Lambda-tau
          {
            xi-=tau/Lambda;
            r_transition=Lambda*rand(); //random variable for weighted sampling of transitions
            if(r_transition<Beta) //S->I
            {
              r_weightedSampling=std::upper_bound(betas.begin(),betas.end(),r_transition)-1; //find beta^(m) interval corresponding to x_choice
              m=r_weightedSampling-betas.begin(); //corresponding m
              // Add infected node to lists:
              isInfected[si_s[m]]=true;
              infected.push_back(si_s[m]);
              mus.push_back(mu*recoverabilities[si_s[m]]+mus.back());
              Mu=mus.back();
            }
            else //I->S
            {
              r_weightedSampling=std::upper_bound(mus.begin(),mus.end(),r_transition-Beta)-1; //find mu^(m) interval corresponding to x_choice
              m=r_weightedSampling-mus.begin(); //corresponding m
              lambda_draw=mu*recoverabilities[infected[m]]; //extract rate for transition m
              isInfected[infected[m]]=false;
              infected[m]=infected.back(); //remove drawn element from infected
              infected.pop_back();
              lambda_m_new=mu*infectivities[infected[m]];
              r_weightedSampling=r_weightedSampling-1+lambda_m_new; //replace mus[m]
              mus.pop_back();
              // Subtract lambda_draw from and add lambda_m_new to all mus[n] for n>m
              for(mu_iterator=r_weightedSampling+1; mu_iterator!=mus.end(); mu_iterator++) *mu_iterator+=lambda_m_new-lambda_draw;
              Mu=mus.back();
            }
            // Redo list of susceptible nodes in contact with infected nodes to update betas:
            si_s.clear();
            betas.assign(1,0.);
            for(contactList_iterator=(*contact_iterator).begin(); contactList_iterator!=(*contact_iterator).end(); contactList_iterator++)
            {
              i=(*contactList_iterator).i;
              j=(*contactList_iterator).j;
              if(isInfected[i])
              {
                if(!isInfected[j])
                {
                  si_s.push_back(j);
                  betas.push_back(beta*susceptibilities[j]*infectivities[i]+betas.back());
                }
              }
              else
              {
                if(isInfected[j])
                {
                  si_s.push_back(i);
                  betas.push_back(beta*susceptibilities[i]*infectivities[j]+betas.back());
                }
              }
            }
            SI=si_s.size();
            Beta=betas.back();
            Lambda=Beta+Mu; //new cumulative transition rate
            // Draw new renormalized waiting time:
            tau=randexp(1);
          }
          // Update list of number of infected nodes:
          I=infected.size();
        }
        // Stop if I=0:
        if(I==0)
        {
          stopped++;
          break;
        }
        // read out I and R if t is divisible by outputTimeResolution
        if(t % outputTimeResolution ==0)
        {
          if(t>=T_simulation)
          {
            break;
          }
          else
          {
            sumI_t[t/outputTimeResolution]+=I;
          }
        }
        t++;
      }
      t_infectionStart=0;
    }
    hist_I[I]++;
  }

}
}

/***R
#generate the network medium the simulation will be run on
library(igraph)
erdos <- erdos.renyi.game(n=100,p=0.025,directed=TRUE,loops=FALSE)
erdos_edge <- as_edgelist(erdos)
erdos_edge <- erdos_edge[order(erdos_edge[,1]),]


*/
