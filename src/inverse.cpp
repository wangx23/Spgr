#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace arma;
using namespace Rcpp;

//inverse
arma::mat inverse(const arma::mat &x, const arma::mat &z,
                  const double nu)
{
  int n = x.n_rows;
  int p = x.n_cols;
  int q = z.n_cols;

  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat zt = z.t();
  arma::mat matinv = zeros(n0,n0);


  arma::mat DB = zeros(p,p);
  arma::mat zAz = zeros(q, q);
  arma::mat Az = zeros(p,q);
  arma::mat AB = zeros(n0, p);
  arma::mat AZ = zeros(n0,q);

  arma::mat mati = zeros(p,p);

  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){

    mati = inv(xt.col(i) * x.row(i) + nIp);
    DB = DB + mati;
    zAz = zAz + zt.col(i)*z.row(i) -  zt.col(i)*x.row(i) *mati * xt.col(i)*z.row(i);
    Az = Az + mati * xt.col(i)*z.row(i);
    AB.rows(idp1,idp2) = mati;
    AZ.rows(idp1,idp2) = mati*xt.col(i)* z.row(i);
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output

    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }

  arma::mat zAzinv = zAz.i();
  arma::mat Azt = Az.t();
  arma::mat temp1 = inv(Ip - DB - Az * zAzinv* Azt);

  arma::mat A1 = AB * temp1 * AB.t();
  arma::mat A2 = AZ * zAzinv * Az.t();
  arma::mat A3 = A2 * temp1;
  arma::mat A22 = A3 * A2.t();
  arma::mat A12 = A3 * AB.t();

  matinv = matinv + AZ * zAzinv * AZ.t() + A1 + A22 + A12 + A12.t();
  return(matinv);

}

// z is not included
arma::mat inversex(const arma::mat &x, const double nu)
{
  int n = x.n_rows;
  int p = x.n_cols;
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat matinv = zeros(n0,n0); 
  
  
  arma::mat DB = zeros(p,p);
  arma::mat AB = zeros(n0, p);
  arma::mat mati = zeros(p,p);
  
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    mati = inv(xt.col(i) * x.row(i) + nIp);
    DB = DB + mati;
    AB.rows(idp1,idp2) = mati;
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  arma::mat IB = inv(Ip - DB);
  
  matinv = matinv + AB * IB * AB.t();
  return(matinv);
  
}


// sparse matrix
arma::sp_mat Dfun(int n)
{
  int npair = n*(n-1)/2;
  int index1;
  arma::mat loc = zeros(2,npair*2);


  for(int j = 0; j < n-1; j++)
  {
    index1 = (2*n -j- 1)*j/2;
    loc(0,span(index1,index1 + n - j - 2)) = linspace<rowvec>(index1,index1+n-j-2,n-j-1);
    loc(1,span(index1,index1 + n - j - 2)) = j*ones<rowvec>(n-j-1);

    loc(0,span(index1 + npair,index1 + n - j - 2 + npair)) = linspace<rowvec>(index1,index1+n-j-2,n-j-1);
    loc(1,span(index1 + npair,index1 + n - j - 2 + npair)) = linspace<rowvec>(j+1,n -1,n-j-1);
  }

  arma::umat loc1 =  arma::conv_to<arma::umat>::from(loc);

  arma::vec values = zeros<vec>(npair*2);
  values.head(npair) = ones(npair);
  values.tail(npair) = (-1)*ones(npair);

  arma::sp_mat res(loc1,values,npair,n);

  return(res);
}


arma::mat inverser( arma::vec indexy, arma::mat &x, arma::mat &z, 
                    double nu)
{
  int p = x.n_cols;
  int q = z.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat zt = z.t();
  arma::mat matinv = zeros(n0,n0);
  
  
  arma::mat DB = zeros(p,p);
  arma::mat zAz = zeros(q, q);
  arma::mat Az = zeros(p,q);
  arma::mat AB = zeros(n0, p);
  arma::mat AZ = zeros(n0,q);
  
  arma::mat mati = zeros(p,p);
  
  arma::uvec indexi;
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    indexi = find(indexy == uindexy(i));
    
    mati = inv(xt.cols(indexi) * x.rows(indexi) + nIp);
    DB = DB + mati;
    zAz = zAz + zt.cols(indexi)*z.rows(indexi) -  
      zt.cols(indexi)*x.rows(indexi) *mati * xt.cols(indexi)*z.rows(indexi);
    
    Az = Az + mati * xt.cols(indexi)*z.rows(indexi);
    AB.rows(idp1,idp2) = mati;
    AZ.rows(idp1,idp2) = mati*xt.cols(indexi)* z.rows(indexi);
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  arma::mat zAzinv = zAz.i();
  arma::mat Azt = Az.t();
  arma::mat temp1 = inv(Ip - DB - Az * zAzinv* Azt);
  
  arma::mat A1 = AB * temp1 * AB.t();
  arma::mat A2 = AZ * zAzinv * Az.t();
  arma::mat A3 = A2 * temp1;
  arma::mat A22 = A3 * A2.t();
  arma::mat A12 = A3 * AB.t();
  
  matinv = matinv + AZ * zAzinv * AZ.t() + A1 + A22 + A12 + A12.t();
  return(matinv);
  
}


arma::mat inverserx( arma::vec indexy, arma::mat &x,double nu)
{
  int p = x.n_cols;
  
  arma::vec uindexy = unique(indexy);
  int n = uindexy.size();
  
  int n0 = n*p;
  arma::mat Ip = 1/nu * eye(p,p);
  arma::mat nIp = nu* n * eye(p,p);
  arma::mat xt = x.t();
  arma::mat matinv = zeros(n0,n0);
  
  
  arma::mat DB = zeros(p,p);
  arma::mat AB = zeros(n0, p);
  arma::mat mati = zeros(p,p);
  
  arma::uvec indexi;
  int idp1 = 0;
  int idp2 = p - 1;
  for(int i=0; i < n; i++){
    
    indexi = find(indexy == uindexy(i));
    
    mati = inv(xt.cols(indexi) * x.rows(indexi) + nIp);
    DB = DB + mati;
    AB.rows(idp1,idp2) = mati;
    matinv.submat(idp1,idp1,idp2,idp2) = mati; // output
    
    idp1 = idp1 + p;
    idp2 = idp2 + p;
  }
  
  
  arma::mat IB = inv(Ip - DB);
  
  matinv = matinv + AB * IB * AB.t();
  return(matinv);
  
}
