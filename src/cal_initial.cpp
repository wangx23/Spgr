#include "functions.hpp"
#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;


// calculate initial value 
// [[Rcpp::export]]
arma::mat cal_initial(arma::vec &y,  arma::mat &z, arma::mat &x,
                      double lam0 = 0.0001)
{
  int n = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  
  int indexk = 0;
  int index1 = 0;
  
  arma::mat ztz = inv(z.t() * z) * z.t();
  arma::vec Xty =  reshape(trans(x.each_col() % y  - 
    x.each_col() %(z*ztz *y)), n0 , 1);
  
  arma::mat Xinv = inverse(x,z,lam0);
  arma::mat betam0 = reshape(Xinv * Xty,p,n);
  arma::vec betamed = zeros(n);
  for(int i=0; i <n; i++)
  {
    betamed(i) = median(betam0.col(i));
  }
  
  int K0 = floor(sqrt(n));
  arma::uvec indexm = sort_index(betamed);
  arma::vec rankm = zeros(n);
  rankm(indexm) = linspace<vec>(0,n-1,n);
  arma::vec breaks = linspace<vec>(0,n - 1, K0+1);
  breaks(K0) = breaks(K0) + 1;
  
  Function f = Environment::base_env()["findInterval"];
  Rcpp::IntegerVector group = f(rankm, breaks);
  arma::uvec group1= as<arma::uvec>(group);

  arma::mat X0 = zeros(n, K0*p);
  arma::uvec colindex = zeros<uvec>(p);
  for(int k=0; k < K0; k++)
  {
    colindex = linspace<uvec>(p*k, p*(k+1) - 1, p);
    X0(find(group1==k+1),colindex) = x.rows(find(group1 == k+1));
  }
  
  arma::mat W0 = zeros(n, q + K0*p);
  int colw = W0.n_cols;
  W0.cols(0,q-1) = z;
  W0.cols(q, colw - 1) = X0;
  arma::mat est0 = inv(W0.t() * W0)*W0.t()*y;
  arma::mat alphaest0 = trans(reshape(est0.rows(q,colw-1),p,K0));
  arma::mat betam = alphaest0.rows(group1 - 1);
  return(betam);
}


// z is not included 
// [[Rcpp::export]]
arma::mat cal_initialx(arma::vec &y, arma::mat &x, double lam0 = 0.0001)
{
  int n = x.n_rows;
  int p = x.n_cols;
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  
  int indexk = 0;
  int index1 = 0;
  
  arma::vec Xty =  reshape(trans(x.each_col() % y), n0 , 1);
  
  arma::mat Xinv = inversex(x,lam0);
  arma::mat betam0 = reshape(Xinv * Xty,p,n);
  arma::vec betamed = zeros(n);
  for(int i=0; i <n; i++)
  {
    betamed(i) = median(betam0.col(i));
  }
  
  int K0 = floor(sqrt(n));
  arma::uvec indexm = sort_index(betamed);
  arma::vec rankm = zeros(n);
  rankm(indexm) = linspace<vec>(0,n-1,n);
  arma::vec breaks = linspace<vec>(0,n - 1, K0+1);
  breaks(K0) = breaks(K0) + 1;
  
  Function f = Environment::base_env()["findInterval"];
  Rcpp::IntegerVector group = f(rankm, breaks);
  arma::uvec group1= as<arma::uvec>(group);
  
  arma::mat X0 = zeros(n, K0*p);
  arma::uvec colindex = zeros<uvec>(p);
  for(int k=0; k < K0; k++)
  {
    colindex = linspace<uvec>(p*k, p*(k+1) - 1, p);
    X0(find(group1==k+1),colindex) = x.rows(find(group1 == k+1));
  }
  

  arma::mat est0 = inv(X0.t() * X0)*X0.t()*y;
  arma::mat alphaest0 = trans(reshape(est0,p,K0));
  arma::mat betam = alphaest0.rows(group1 - 1);
  return(betam);
}


// [[Rcpp::export]]
arma::mat cal_initialr(arma::vec indexy, arma::vec &y,arma::mat &z, arma::mat &x,
                       double lam0 = 0.0001)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  //transformation for y,z,x 
  arma::uvec nJ = zeros<uvec>(n);
  arma::uvec indexi;
  int ni;
  
  arma::mat xm = zeros(nt,p);
  arma::mat zm = zeros(nt,q);
  arma::vec ym = zeros(nt);
  for(int i = 0; i <n ; i++)
  {
    indexi = find(indexy == uindexy(i));
    ni = indexi.size();
    nJ(i) = ni;
    xm.rows(indexi) = x.rows(indexi)/sqrt(ni);
    zm.rows(indexi) = z.rows(indexi)/sqrt(ni);
    ym(indexi) = y(indexi)/sqrt(ni);
  }
  
  arma::mat ztz = inv(zm.t() * zm) * zm.t();
  arma::mat tempxy =  trans(xm.each_col() % ym  -
    xm.each_col() %(zm*ztz *ym));
  
  arma::vec Xty(n0);
  
  for(int i = 0 ; i < n; i++ )
  {
    indexi = find(indexy == uindexy(i));
    Xty(span(i*p,(i+1)*p - 1)) = sum(tempxy.cols(indexi),1);
  }
  
  arma::mat Xinv = inverser(indexy,xm,zm,lam0);
  arma::mat reg1 = Xinv * Xty;
  arma::mat betam = trans(reshape(reg1, p, n));
  return(betam);
}


// [[Rcpp::export]]
arma::mat cal_initialrx(arma::vec indexy, arma::vec &y, arma::mat &x,
                        double lam0 = 0.0001)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  //transformation for y,z,x 
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
  
  arma::mat Xinv = inverserx(indexy,xm,lam0);
  arma::mat reg1 = Xinv * Xty;
  arma::mat betam = trans(reshape(reg1, p, n));
  return(betam);
}


// [[Rcpp::export]]
arma::mat cal_initialr2(arma::vec indexy, arma::vec &y,arma::mat &z, arma::mat &x, int K0,
                        double lam0 = 0.0001)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  //transformation for y,z,x 
  arma::uvec nJ = zeros<uvec>(n);
  arma::uvec indexi;
  int ni;
  
  arma::mat xm = zeros(nt,p);
  arma::mat zm = zeros(nt,q);
  arma::vec ym = zeros(nt);
  for(int i = 0; i <n ; i++)
  {
    indexi = find(indexy == uindexy(i));
    ni = indexi.size();
    nJ(i) = ni;
    xm.rows(indexi) = x.rows(indexi)/sqrt(ni);
    zm.rows(indexi) = z.rows(indexi)/sqrt(ni);
    ym(indexi) = y(indexi)/sqrt(ni);
  }
  
  arma::mat ztz = inv(zm.t() * zm) * zm.t();
  arma::mat tempxy =  trans(xm.each_col() % ym  -
    xm.each_col() %(zm*ztz *ym));
  
  arma::vec Xty(n0);
  
  for(int i = 0 ; i < n; i++ )
  {
    indexi = find(indexy == uindexy(i));
    Xty(span(i*p,(i+1)*p - 1)) = sum(tempxy.cols(indexi),1);
  }
  
  arma::mat Xinv = inverser(indexy,xm,zm,lam0);
  arma::mat reg1 = Xinv * Xty;
  arma::mat betam0 = reshape(reg1, p, n);
  
  arma::vec betamed = zeros(n);
  for(int i=0; i <n; i++)
  {
    betamed(i) = median(betam0.col(i));
  }
  
  arma::uvec indexm = sort_index(betamed);
  arma::vec rankm = zeros(n);
  rankm(indexm) = linspace<vec>(0,n-1,n);
  arma::vec breaks = linspace<vec>(0,n - 1, K0+1);
  breaks(K0) = breaks(K0) + 1;
  
  Function f = Environment::base_env()["findInterval"];
  Rcpp::IntegerVector group = f(rankm, breaks);
  arma::uvec group1= as<arma::uvec>(group);
  
  arma::uvec groupall = zeros<arma::uvec>(nt);
  for(int i=0; i < n; i++)
  {
    indexi = find(indexy == uindexy(i));
    groupall(indexi) = group1(i)*ones<arma::uvec>(nJ(i));
  }
  
  arma::mat X0 = zeros(nt, K0*p);
  arma::uvec colindex = zeros<uvec>(p);
  for(int k=0; k < K0; k++)
  {
    colindex = linspace<uvec>(p*k, p*(k+1) - 1, p);
    X0(find(groupall==k+1),colindex) = xm.rows(find(groupall == k+1));
  }
  
  arma::mat W0 = zeros(nt, q + K0*p);
  int colw = W0.n_cols;
  W0.cols(0,q-1) = zm;
  W0.cols(q, colw - 1) = X0;
  arma::mat est0 = inv(W0.t() * W0)*W0.t()*ym;
  arma::mat alphaest0 = trans(reshape(est0.rows(q,colw-1),p,K0));
  arma::mat betam = alphaest0.rows(group1 - 1);
  return(betam);
}


// [[Rcpp::export]]
arma::mat cal_initialrx2(arma::vec indexy, arma::vec &y, arma::mat &x, int K0,
                         double lam0 = 0.0001)
{
  int nt = x.n_rows;
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  int npair = n*(n-1)/2;
  arma::mat Ip =  eye(p,p);
  
  arma::sp_mat Dmat = Dfun(n);
  
  //transformation for y,z,x 
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
  
  arma::mat Xinv = inverserx(indexy,xm,lam0);
  arma::mat reg1 = Xinv * Xty;
  arma::mat betam0 = reshape(reg1, p, n);
  
  arma::vec betamed = zeros(n);
  for(int i=0; i <n; i++)
  {
    betamed(i) = median(betam0.col(i));
  }
  
  arma::uvec indexm = sort_index(betamed);
  arma::vec rankm = zeros(n);
  rankm(indexm) = linspace<vec>(0,n-1,n);
  arma::vec breaks = linspace<vec>(0,n - 1, K0+1);
  breaks(K0) = breaks(K0) + 1;
  
  Function f = Environment::base_env()["findInterval"];
  Rcpp::IntegerVector group = f(rankm, breaks);
  arma::uvec group1= as<arma::uvec>(group);
  
  arma::uvec groupall = zeros<arma::uvec>(nt);
  for(int i=0; i < n; i++)
  {
    indexi = find(indexy == uindexy(i));
    groupall(indexi) = group1(i)*ones<arma::uvec>(nJ(i));
  }
  
  arma::mat X0 = zeros(nt, K0*p);
  arma::uvec colindex = zeros<uvec>(p);
  for(int k=0; k < K0; k++)
  {
    colindex = linspace<uvec>(p*k, p*(k+1) - 1, p);
    X0(find(groupall==k+1),colindex) = xm.rows(find(groupall == k+1));
  }
  
  arma::mat est0 = inv(X0.t() * X0)*X0.t()*ym;
  arma::mat alphaest0 = trans(reshape(est0,p,K0));
  arma::mat betam = alphaest0.rows(group1 - 1);
  return(betam);
}