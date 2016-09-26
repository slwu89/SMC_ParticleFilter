#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <boost/function.hpp>
using namespace Rcpp;



/*
 * This file includes functions to run stochastic simulation algorithms on static networks.
 */

// node_tuple is a structure that defines the state of a pair of nodes
struct node_tuple {int i; int j; char i_state; char j_state;};
// node_tuple_list
typedef std::vector<node_tuple> node_tuple_list;




/*
 * update_si updates m_SI; the takes a node_tuple_list at time t and returns a vector of 
 * indicies where either the i_state or j_state of the node tuple at that index has one infected. 
 */

/***R
#generate the network medium the simulation will be run on
library(igraph)
erdos <- erdos.renyi.game(n=100,p=0.025,directed=TRUE,loops=FALSE)
erdos_edge <- as_edgelist(erdos)
erdos_edge <- erdos_edge[order(erdos_edge[,1]),]


*/
