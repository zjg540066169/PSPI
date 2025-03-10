#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RCPP_H_
#define RCPP_H_
#include <Rcpp.h>
#endif

#ifndef NS_BASIS_H_
#define NS_BASIS_H_
#include "NS_basis.h"
#endif

#ifndef RCPPDIST_H_
#define RCPPDIST_H_
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#endif

#include <algorithm>

using namespace Rcpp;
using namespace arma;



class NS{
public:
  NS(){};
  
  NS(NumericVector X, NumericVector y, long K, double sigma, int order = 3){
    this->X = X;
    this->y = as<arma::vec>(y);
    this->K = K;
    this->n = X.length();
    this->sigma = sigma;
    
    basis = new NS_basis(X, K, order);
    complete_basis = as<arma::mat>(basis->get_basis());
    if(K > 2)
      ns_basis = as<arma::mat>(basis->get_ns_part());
    lr_basis = complete_basis.cols(0, 1);
    arma::vec theta_init = update_beta_mean(complete_basis, this->y);
  
    lr_coefficient = theta_init.subvec(0, 1);
    
    if(K > 2){
      ns_coefficient = theta_init.subvec(2, K - 1);
    }
  };
  
  virtual void update() = 0;
  virtual void update(double sigma)= 0;
  virtual void update(NumericVector sigma)= 0;
  
  void set_Y(NumericVector y){
    this->y = as<arma::vec>(y);
  }
  
  NumericVector predict(NumericVector X_test){
    complete_basis_test = as<arma::mat>(basis->predict(X_test));
    if(K > 2)
      ns_basis_test = complete_basis_test.cols(2, K-1);
    lr_basis_test = complete_basis_test.cols(0, 1);
    arma::vec ns_outcome_test = lr_basis_test * lr_coefficient;
    if(K > 2)  
      ns_outcome_test = ns_outcome_test + ns_basis_test * ns_coefficient;
    return wrap(ns_outcome_test);
  }
  
  NumericMatrix get_basis(){
    return wrap(complete_basis);
  }
  
  NumericMatrix get_basis_test(){
    return wrap(complete_basis_test);
  }
  
  NumericMatrix get_ns_part(){
    return wrap(ns_basis);
  }
  
  NumericMatrix get_lr_part(){
    return wrap(lr_basis);
  }
  
  NumericVector get_boundary_knots(){
    return basis->get_boundary_knots();
  }
  
  NumericVector get_internal_knots(){
    return basis->get_internal_knots();
  }
  
  NumericVector get_knots(){
    return basis->get_knots();
  }
  
  NumericVector get_theta(){
    return wrap(arma::join_cols(lr_coefficient, ns_coefficient));
  }
  
  NumericVector get_ns_outcome(){
    return wrap(ns_outcome);
  }
  
  double get_gamma(){
    return gamma;
  }
  
  arma::mat inv_X_T_X(arma::mat mat_X){
    return(arma::inv_sympd(mat_X.t() * mat_X));
  }
  
  arma::vec update_beta_mean(arma::mat mat_X, arma::vec mat_Y){
    return(inv_X_T_X(mat_X) * mat_X.t() * mat_Y);
  }
  
  arma::vec update_beta(arma::mat mat_X, arma::vec mat_Y, double sigma){
    vec beta_mean = update_beta_mean(mat_X, mat_Y);
    mat beta_var = inv_X_T_X(mat_X) * sigma * sigma;
    return arma::vectorise(rmvnorm(1, beta_mean, beta_var));
  }
  
  arma::vec update_beta(arma::mat mat_X, arma::vec mat_Y, arma::vec sigma) {
    this->Sigma = arma::diagmat(sigma % sigma);
    //Rcout<<this->Sigma<<std::endl;
    // Form the diagonal weight matrix.
    inv_Sigma = arma::inv_sympd(Sigma);
    
    // Compute the weighted cross-product matrix.
    arma::mat XtWX = mat_X.t() * inv_Sigma * mat_X;
    
    // The posterior variance for beta.
    arma::mat beta_var = arma::inv_sympd(XtWX);
    
    // The posterior mean for beta (weighted least squares estimate).
    arma::vec beta_mean = beta_var * (mat_X.t() * inv_Sigma * mat_Y);
    
    // Sample from the multivariate normal with mean beta_mean and covariance beta_var.
    return arma::vectorise(rmvnorm(1, beta_mean, beta_var));
  }
  
  // arma::vec update_beta(arma::mat mat_X, arma::vec mat_Y, arma::mat Sigma) {
  //   this->Sigma = Sigma;
  //   
  //   // Form the diagonal weight matrix.
  //   arma::mat W = arma::inv_sympd(this->Sigma);
  //   
  //   // Compute the weighted cross-product matrix.
  //   arma::mat XtWX = mat_X.t() * W * mat_X;
  //   
  //   // The posterior variance for beta.
  //   arma::mat beta_var = arma::inv(XtWX);
  //   
  //   // The posterior mean for beta (weighted least squares estimate).
  //   arma::vec beta_mean = beta_var * (mat_X.t() * W * mat_Y);
  //   
  //   // Sample from the multivariate normal with mean beta_mean and covariance beta_var.
  //   return arma::vectorise(rmvnorm(1, beta_mean, beta_var));
  // }
  
  
  
  void update_outcome(){
    lr_outcome = lr_basis * lr_coefficient;
    ns_outcome = lr_outcome;
    if(K > 2){
      ns_part_outcome = ns_basis * ns_coefficient;
      ns_outcome = ns_outcome + ns_part_outcome;
    }
    //Rcout << ns_outcome << std::endl;
  }
  
  double rinvgamma(double a, double b){
    double s = R::rgamma(a, 1 / b);
    return 1 / s;
  }
  
  double qinvgamma(double p, double a, double b) {
    // R::qgamma returns the quantile for the Gamma distribution with given shape and scale.
    double q = R::qgamma(1.0 - p, a, 1.0 / b, true, false);
    return 1.0 / q;
  }
  
  
protected:
  NS_basis * basis;
  
  long K;
  long n;
  double sigma;
  double gamma;
  
  arma::mat Sigma;
  arma::mat inv_Sigma;
  
  arma::vec lr_coefficient;
  arma::vec ns_coefficient;
  
  
  NumericVector X;
  arma::vec y;
  
  arma::mat complete_basis;
  arma::mat ns_basis;
  arma::mat lr_basis;
  
  arma::mat complete_basis_test;
  arma::mat ns_basis_test;
  arma::mat lr_basis_test;
  
  arma::vec ns_outcome;
  arma::vec lr_outcome;
  arma::vec ns_part_outcome;
  
};



// 
// // [[Rcpp::export]]
// List test_NS(NumericVector X, NumericVector X_test, NumericVector y, long K){
//   NS * a = new NS(X, y, K, 1);
//   Rcout << a->get_theta() << std::endl;
//   for(int i = 1; i < 10000; ++i)
//     a->update();
//   Rcout << a->get_theta() << std::endl;
//   return List::create(Named("ns_predict") = a->predict(X_test), Named("eta") = a->get_eta(), Named("gamma") = a->get_gamma(), Named("theta") = a->get_theta());
//   //return List::create(Named("ns_outcome") = wrap(a->ns_outcome), Named("lr_outcome") = wrap(a->lr_outcome), Named("ns_part_outcome") = wrap(a->ns_part_outcome), Named("eta") = a->get_eta(), Named("gamma") = a->get_gamma(), Named("theta") = a->get_theta(), Named("boundary_knots") = a->get_boundary_knots(), Named("knots") = a->get_knots(), Named("internal_knots") = a->get_internal_knots(), Named("ns_part") = a->get_ns_part(), Named("lr_part") = a->get_lr_part(), Named("basis") = a->get_basis());
// };
