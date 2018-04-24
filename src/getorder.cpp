#include <RcppArmadillo.h>
#include "functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// This is a function to give the order vector 
// [[Rcpp::export]]
arma::uvec getorder(arma::sp_mat &Cmat)
{
  arma::mat Cm(Cmat);
  int n = Cm.n_cols;
  arma::uvec ordervec = zeros<uvec>(n*(n-1)/2);
  arma::uvec orderj = zeros<uvec>(n);
  arma::uvec indexknd;
  arma::uvec locj;
  arma::uvec order0;
  arma::uvec orderl;
  
  int indexj1;
  int indexj2;
  int k = 0;
  
  for(int j =0; j < n - 1; j++)
  {
    indexj1 = (2*n -j- 1)*j/2;
    indexj2 = indexj1 + n - j - 2;
    
    orderj = arma::conv_to<arma::uvec>::from((Cm.col(j)));
    orderj(j) = 1;
    
    k = 0;
    do{
      k = k + 1;
      indexknd = find(sum(Cm.cols(find(orderj == k)),1) >0);
      order0 = find(orderj==0);
      locj  = indexknd( find( findin(indexknd, order0) == 1));
      orderj(locj) =  (k+1)*ones<uvec>(locj.size());
      orderl = find(orderj == 0);
      
    } while (orderl.size()>0);
    
    ordervec.subvec(indexj1,indexj2) = orderj.subvec(j + 1,n-1);
  }
  
  return(ordervec);
}