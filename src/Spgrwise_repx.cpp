#include <RcppArmadillo.h>
#include "functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
using namespace Rcpp;
using namespace arma; // repeated measure

// [[Rcpp::export]]
Rcpp::List Spgrwise_repx_lasso(arma::vec indexy,arma::vec &y, arma::mat &x, 
                 arma::vec &weights, arma::mat &betam0,
                 double nu = 1 , double lam = 0.5 ,
                 int maxiter = 1000, double tolabs = 1e-4, double tolrel = 1e-2)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  //transformation for y,x 
  arma::uvec nJ = zeros<uvec>(n);
  arma::uvec indexi;
  int ni;
  
  arma::mat xm = zeros(nt,p);
  arma::vec ym = zeros(nt);
  for(int i = 0; i <n ; i++)
  {
    indexi = find(indexy == uindexy(i));
    ni = indexi.size();
    nJ(i) = ni;
    xm.rows(indexi) = x.rows(indexi)/sqrt(ni);
    ym(indexi) = y(indexi)/sqrt(ni);
  }
  

  arma::mat tempxy =  trans(xm.each_col() % ym);
  
  arma::vec Xty(n0);
  
  for(int i = 0 ; i < n; i++ )
  {
    indexi = find(indexy == uindexy(i));
    Xty(span(i*p,(i+1)*p - 1)) = sum(tempxy.cols(indexi),1);
  }
  
  arma::mat Xinv = inverserx(indexy,xm,nu);
  arma::mat reg1 = Xinv * Xty;
  
  //initial deltam
  arma::mat deltam(p,npair);
  arma::mat deltamold(p, npair);
  arma::mat betadiff(p,npair);
  
  arma::mat vm = zeros(p,npair);
  deltamold = trans(Dmat * betam0);
  
  // define some variables
  
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
      for(int l = 0; l < p; l++)
      {
        deltam(l,i) = nsfun(deltam(l,i),weights(i)*lam);
      }
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
  
  arma::vec yhat = zeros(nt);
  arma::vec yhatm = zeros(nt);
  for(int i = 0; i < n; i++)
  {
    indexi = find(indexy == uindexy(i));
    yhatm(indexi) = xm.rows(indexi)*trans(betam.row(i));
    yhat(indexi) = yhatm(indexi) * sqrt(nJ(i));
  }
  
 
  // change deltam
  deltam = deltam.t();
  
  arma::umat group(n,p);
  
  for(int l=0;l <p; l++)
  {
    group.col(l) = ngetgroup(deltam.col(l),n,tolabs);
  }

  double sig2 = sum(square(y  - yhat))/nt;
  
  
  return Rcpp::List::create(Named("beta") = betam,
                            Named("index") = uindexy,
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




// [[Rcpp::export]]
Rcpp::List Spgrwise_repx_scad(arma::vec indexy,arma::vec &y, arma::mat &x, 
                              arma::vec &weights, arma::mat &betam0,
                              double nu = 1, double gam = 3, double lam = 0.5 ,
                              int maxiter = 1000, double tolabs = 1e-4, double tolrel = 1e-2)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  //transformation for y,x 
  arma::uvec nJ = zeros<uvec>(n);
  arma::uvec indexi;
  int ni;
  
  arma::mat xm = zeros(nt,p);
  arma::vec ym = zeros(nt);
  for(int i = 0; i <n ; i++)
  {
    indexi = find(indexy == uindexy(i));
    ni = indexi.size();
    nJ(i) = ni;
    xm.rows(indexi) = x.rows(indexi)/sqrt(ni);
    ym(indexi) = y(indexi)/sqrt(ni);
  }
  
  arma::mat tempxy =  trans(xm.each_col() % ym);
  
  arma::vec Xty(n0);
  
  for(int i = 0 ; i < n; i++ )
  {
    indexi = find(indexy == uindexy(i));
    Xty(span(i*p,(i+1)*p - 1)) = sum(tempxy.cols(indexi),1);
  }
  
  arma::mat Xinv = inverserx(indexy,xm,nu);
  arma::mat reg1 = Xinv * Xty;
  
  //initial deltam
  arma::mat deltam(p,npair);
  arma::mat deltamold(p, npair);
  arma::mat betadiff(p,npair);
  
  arma::mat vm = zeros(p,npair);
  deltamold = trans(Dmat * betam0);
  
  // define some variables
  
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
      for(int l = 0; l < p; l++)
      {
        deltam(l,i) = nscad(deltam(l,i),weights(i)*lam);
      }
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
  
  arma::vec yhat = zeros(nt);
  arma::vec yhatm = zeros(nt);
  for(int i = 0; i < n; i++)
  {
    indexi = find(indexy == uindexy(i));
    yhatm(indexi) = xm.rows(indexi)*trans(betam.row(i));
    yhat(indexi) = yhatm(indexi) * sqrt(nJ(i));
  }
  
  
  // change deltam
  deltam = deltam.t();
  
  arma::umat group(n,p);
  
  for(int l=0;l <p; l++)
  {
    group.col(l) = ngetgroup(deltam.col(l),n,tolabs);
  }
  
  double sig2 = sum(square(y  - yhat))/nt;
  
  
  return Rcpp::List::create(Named("beta") = betam,
                            Named("index") = uindexy,
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