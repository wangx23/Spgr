#include <RcppArmadillo.h>
#include "functions.hpp"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]] 
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List Spgrwise_lasso(arma::vec &y, arma::mat &z, arma::mat &x,
                arma::vec &weights, arma::mat &betam0,
                double nu = 1 , double lam = 0.5 ,
                int maxiter = 1000, double tolabs = 1e-4, double tolrel = 1e-2)
{
  int n = x.n_rows;
  int p = x.n_cols;
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);


  arma::sp_mat Dmat = Dfun(n);

  arma::mat ztz = inv(z.t() * z) * z.t();
  arma::vec Xty =  reshape(trans(x.each_col() % y  -
    x.each_col() %(z*ztz *y)), n0 , 1);

  arma::mat Xinv = inverse(x,z,nu);
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
  arma::umat group(n,p);

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

  arma::vec yhat = sum(x % betam,1);
  arma::vec eta =  ztz * (y - yhat);

  // change deltam
  deltam = deltam.t();


  for(int l=0;l <p; l++)
  {
    group.col(l) = ngetgroup(deltam.col(l),n,tolabs);
  }


  double sig2 = sum(square(y - z * eta - yhat))/n;

  return Rcpp::List::create(Named("beta") = betam,
                            Named("eta") = eta,
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
Rcpp::List Spgrwise_scad(arma::vec &y, arma::mat &z, arma::mat &x,
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

  arma::mat ztz = inv(z.t() * z) * z.t();
  arma::vec Xty =  reshape(trans(x.each_col() % y  -
    x.each_col() %(z*ztz *y)), n0 , 1);

  arma::mat Xinv = inverse(x,z,nu);
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
  arma::umat group(n,p);

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
        deltam(l,i) = nscad(deltam(l,i),weights(i)*lam,nu,gam);
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

  arma::vec yhat = sum(x % betam,1);
  arma::vec eta =  ztz * (y - yhat);

  // change deltam
  deltam = deltam.t();


  for(int l=0;l <p; l++)
  {
    group.col(l) = ngetgroup(deltam.col(l),n,tolabs);
  }


  double sig2 = sum(square(y - z * eta - yhat))/n;

  return Rcpp::List::create(Named("beta") = betam,
                            Named("eta") = eta,
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
