// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,  
                           arma::rowvec const &mean,  
                           arma::mat const &sigma, 
                           bool const logd = false) { 
  using arma::uword;
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}

arma::vec Mahalanobis(arma::mat const &x, 
                      arma::vec const &center, 
                      arma::mat const &cov) {
  arma::mat x_cen = x.t();
  x_cen.each_col() -= center;
  arma::solve(x_cen, arma::trimatl(chol(cov).t()), x_cen);
  x_cen.for_each( [](arma::mat::elem_type& val) { val = val * val; } );
  return arma::sum(x_cen, 0).t();    
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat const &x, 
                       arma::vec const &mean, 
                       arma::mat const &sigma, 
                       bool const logd = false) { 
  arma::vec const distval = Mahalanobis(x,  mean, sigma);
  double const logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec const logretval = 
    -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
    
    if (logd)
      return logretval;
    return exp(logretval);
}

// [[Rcpp::export]]
arma::mat mvrnorm_arma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}