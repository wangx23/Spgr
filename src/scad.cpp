#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// sfun
arma::vec sfun(arma::vec x, double th)
{
  double xn = norm(x,2);
  double thval = 1 - th/xn;
  x = thval*(thval >0)*x;
  return(x);
}


//scad
arma::vec scad(arma::vec v, double lam, const double nu =1 ,
               const double gam  = 3)
{

  double temp1 = lam/nu;
  double temp2 = gam * lam;
  double xn = norm(v,2);
  arma::vec vnew;

  if(xn <= lam + temp1)
  {
    vnew = sfun(v, temp1);
  }else if(xn <= temp2 & xn >= lam + temp1)
  {
    vnew = sfun(v, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu));
  }else{
    vnew = v;
  }
  return(vnew);
}
