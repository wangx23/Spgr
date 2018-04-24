#include <RcppArmadillo.h>
#include "functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List Spgrx(arma::vec &y, arma::mat &x, 
                         arma::vec &weights, arma::mat &betam0,
                         double nu = 1, double gam =3 , double lam = 0.5 ,
                         int maxiter = 1000, double tolabs = 1e-4, double tolrel = 1e-2)
{
  int n = x.n_rows;
  int p = x.n_cols;
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);

  arma::sp_mat Dmat = Dfun(n);
  
  arma::vec Xty =  reshape(trans(x.each_col() % y), n0 , 1);
  
  arma::mat Xinv = inversex(x,nu);
  arma::mat reg1 = Xinv * Xty;
  
  // initial deltam 
  arma::mat deltam(p,npair);
  arma::mat deltamold(p, npair);
  arma::mat betadiff(p,npair);
  
  arma::mat vm = zeros(p,npair);
  deltamold = trans(Dmat * betam0);
  
  // define some varaibles 
  
  arma::vec temp = zeros<vec>(n0);
  arma::vec betanew = zeros<vec>(n0);
  arma::mat betam = zeros(n,p);
  arma::vec normbd(2);
  
  int flag = 0;
  double rm  = 1;
  double sm = 1;
  double tolpri;
  double toldual;
  
  
  int m = 0;
  
  for( m = 0; m < maxiter; m++)
  {
    
    temp =  reshape((deltamold - 1/nu * vm)*Dmat,n0,1);
    betanew =  reg1 + nu*Xinv * temp;
    betam = trans(reshape(betanew, p, n));
    betadiff = trans(Dmat * betam);
    
    deltam = betadiff + (1/nu) * vm;
    
    // update deltam 
    for(int i = 0; i < npair; i++)
    {
      deltam.col(i) = scad(deltam.col(i),weights(i)*lam,nu,gam);
    }
    
    vm =  vm + nu * (betadiff - deltam);
    
    normbd(0) = norm(betadiff,"fro");
    normbd(1) = norm(deltam,"fro");
    
    tolpri = tolabs*sqrt(npair*p) + tolrel*max(normbd);
    toldual = tolabs*sqrt(n * p) + tolrel * norm(vm * Dmat, "fro");
    
    rm = norm(betadiff - deltam, "fro");
    sm = nu * norm((deltam - deltamold)*Dmat, "fro");
    
    deltamold = deltam;
    
    if(rm <= tolpri & sm <= toldual)
      break;
  }
  
  if(m == maxiter) {flag = 1;}
  
  arma::vec yhat = sum(x % betam,1);

  arma::uvec group = getgroup(deltam,n,tolabs);
  arma::uvec ugroup = unique(group);
  int ng = ugroup.size();
  int nj; 
  
  arma::uvec indexj;
  arma::mat betaest(ng,p);
  
  
  for(int j = 0; j < ng; j++)
  {
    indexj = find(group == ugroup(j));
    nj = indexj.size();
    betaest.row(j) = sum(betam.rows(indexj),0)/nj;
  }
  
  double sig2 = sum(square(y  - yhat))/(n  - ng*p);
  
  
  
  return Rcpp::List::create(Named("beta") = betam,
                            Named("betaest") = betaest,
                            Named("sig2") = sig2,
                            Named("group") = group,
                            Named("deltam") = deltam,
                            Named("flag") = flag,
                            Named("rm") = rm,
                            Named("sm") = sm,
                            Named("tolpri") = tolpri,
                            Named("toldual") = toldual,
                            Named("niteration") = m);
  
}


//[[Rcpp::export]]
Rcpp::List selectlamx(arma::vec &y, arma::mat &x, arma::vec &weights,
                     arma::vec &lamv,arma::mat &betam0, double nu = 1, double gam = 3,
                     int maxiter = 1000, double tolabs = 1e-4, double tolrel = 1e-2)
{
  int nlam = lamv.size();
  int p = x.n_cols;
  int n = x.n_rows;
  
  // initial value
  
  arma::cube betacube(n,p,nlam);
  arma::vec sig2(nlam);
  arma::umat groupm(n, nlam);
  arma::vec bic(nlam);
  arma::vec flag(nlam);
  arma::uvec niterationv(nlam);
  
  Rcpp::List res;
  
  arma::uvec groupj(n);
  arma::mat betaestj = betam0;
  
  for(int j1 = 0; j1 < nlam; j1++)
  {
    res = Spgrx(y, x, weights, betaestj, nu, gam, lamv(j1), maxiter, tolabs, tolrel);
    betaestj = Rcpp::as<arma::mat>(res("beta"));
    groupj  = Rcpp::as<arma::uvec>(res("group"));
    
    betacube.slice(j1) = betaestj;
    sig2(j1) = Rcpp::as<double>(res("sig2"));
    groupm.col(j1) = groupj;
    flag(j1) = Rcpp::as<double>(res("flag"));
    bic(j1) = BIClogx(res,y,x, tolabs);
    niterationv(j1) = res("niteration");
  }
  
  uword index =  bic.index_min();
  
  arma::uvec group = groupm.col(index);
  arma::uvec ugroup = unique(group);
  int ng = ugroup.size();
  int nj; 
  
  arma::uvec indexj;
  arma::mat betaest(ng,p);
  arma::mat betam = betacube.slice(index);
  
  for(int j = 0; j < ng; j++)
  {
    indexj = find(group == ugroup(j));
    nj = indexj.size();
    betaest.row(j) = sum(betam.rows(indexj),0)/nj;
  }
  
  
  return Rcpp::List::create(Named("beta") = betam,
                            Named("betaest") = betaest,
                            Named("sig2") = sig2(index),
                            Named("group") = group,
                            Named("convergence") = flag(index),
                            Named("niteration") = niterationv(index),
                            Named("lambda") = lamv(index),
                            Named("BIC") = bic(index),
                            Named("bicv") = bic,
                            Named("betacube") = betacube,
                            Named("groupm") = groupm);
}


