#include <RcppArmadillo.h>
#include "functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// a function to implement in for uvec index of a%in%b
arma::uvec findin(arma::uvec &a, arma::uvec &b) 
{
  int asize = a.size();
  arma::uvec indexa = zeros<uvec>(asize);
  arma::uvec bi;
  
  for(int i=0; i<asize; i++)
  {
    bi = find( b == a(i));
    if(bi.size() >0){indexa(i) = 1;}
  }
  return(indexa);
}
