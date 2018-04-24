#include <RcppArmadillo.h>
#include "functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// obtain group estimate
// [[Rcpp::export]]
arma::uvec getgroup(const arma::mat &deltam, 
                    const int n, const double tol = 1e-4)
{
  int p = deltam.n_rows;
  arma::vec b2value  =  arma::conv_to<arma::vec>::from(sqrt(sum(square(deltam),0)/p));
  arma::uvec indexb = find(abs(b2value) <= tol);
  b2value(indexb) = zeros(size(indexb));
  
  arma::mat d2 = zeros(n,n);
  int indexj1;
  int indexj2;
  
  for(int j = 0; j < n - 1; j++)
  {
    indexj1 = (2*n -j- 1)*j/2;
    indexj2 = indexj1 + n - j - 2;
    d2(span(n - indexj2 + indexj1 - 1, n - 1),j) = b2value(span(indexj1,indexj2));
  }
  
  d2 = d2.t() + d2;
  
  int j = 0;
  arma::uvec ngj = linspace<uvec>(0, n - 1, n);
  arma::uvec gj;
  arma::uvec groupest = zeros<uvec>(n);
  
   do{
   j = j + 1;
   gj = find(d2.row(ngj(0))==0);
   gj = ngj(find(findin(ngj, gj)==1));
   ngj = ngj(find(findin(ngj,gj)==0));
   groupest(gj) = j*ones<uvec>(gj.size());
   } while (ngj.size()>0);
   
  return(groupest);
}