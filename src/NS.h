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
  
  NS(NumericVector X, NumericVector y, long K, double sigma){
    this->X = X;
    this->y = as<arma::vec>(y);
    this->K = K;
    this->n = X.length();
    this->sigma = sigma;
    basis = new NS_basis(X, K);
    complete_basis = as<arma::mat>(basis->get_basis());
    ns_basis = as<arma::mat>(basis->get_ns_part());
    lr_basis = complete_basis.cols(0, 1);
    arma::vec theta_init = update_beta_mean(complete_basis, this->y);
    
  
    
    lr_coefficient = theta_init.subvec(0, 1);
    ns_coefficient = theta_init.subvec(2, K - 1);
    
    eta = rinvgamma(0.5, 1);
    gamma = sqrt(rinvgamma(0.5, 1 / eta));
    update_outcome();
  };
  
  void set_Y(NumericVector y){
    this->y = as<arma::vec>(y);
  }
  
  void update(){
    lr_coefficient = update_beta(lr_basis, this->y - ns_part_outcome, sigma);
    lr_outcome = lr_basis * lr_coefficient;
    arma::mat theta_var = arma::inv_sympd(1 / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + 1 / pow(sigma, 2) * ns_basis.t() * ns_basis);
    arma::vec theta_mean = arma::inv_sympd(pow(sigma, 2) / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + ns_basis.t() * ns_basis) * ns_basis.t() * (this->y - lr_outcome);
    ns_coefficient = arma::vectorise(rmvnorm(1, theta_mean, theta_var));
    ns_part_outcome = ns_basis * ns_coefficient;
    ns_outcome = lr_outcome + ns_part_outcome;

    gamma = sqrt(rinvgamma((K-1.0) / 2.0, 1 / eta + 0.5 * arma::dot(ns_coefficient, ns_coefficient)));
    eta = rinvgamma(0.5, 1.0 + 1.0 / (gamma * gamma));
  }
  
  void update(double sigma){
    this->sigma = sigma;
    lr_coefficient = update_beta(lr_basis, this->y - ns_part_outcome, this->sigma);
    lr_outcome = lr_basis * lr_coefficient;
    arma::mat theta_var = arma::inv_sympd(1 / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + 1 / pow(this->sigma, 2) * ns_basis.t() * ns_basis);
    arma::vec theta_mean = arma::inv_sympd(pow(this->sigma, 2) / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + ns_basis.t() * ns_basis) * ns_basis.t() * (this->y - lr_outcome);
    ns_coefficient = arma::vectorise(rmvnorm(1, theta_mean, theta_var));
    ns_part_outcome = ns_basis * ns_coefficient;
    ns_outcome = lr_outcome + ns_part_outcome;
    
    gamma = sqrt(rinvgamma((K-1.0) / 2.0, 1 / eta + 0.5 * arma::dot(ns_coefficient, ns_coefficient)));
    eta = rinvgamma(0.5, 1.0 + 1.0 / (gamma * gamma));
  }
  
  NumericVector predict(NumericVector X_test){
    complete_basis_test = as<arma::mat>(basis->predict(X_test));
    ns_basis_test = complete_basis_test.cols(2, K-1);
    lr_basis_test = complete_basis_test.cols(0, 1);
    arma::vec ns_outcome_test = lr_basis_test * lr_coefficient + ns_basis_test * ns_coefficient;
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
  
  double get_eta(){
    return eta;
  }
  
  double get_gamma(){
    return gamma;
  }
  
  NumericVector get_ns_outcome(){
    return wrap(ns_outcome);
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

  void update_outcome(){
    lr_outcome = lr_basis * lr_coefficient;
    ns_part_outcome = ns_basis * ns_coefficient;
    ns_outcome = lr_outcome + ns_part_outcome;
  }
  
  double rinvgamma(double a, double b){
    double s = R::rgamma(a, 1 / b);
    return 1 / s;
  }
  
private:
  NS_basis * basis;
  
  long K;
  long n;
  double sigma;
  
  arma::vec lr_coefficient;
  arma::vec ns_coefficient;
  double gamma;
  double eta;
  
  
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




// // [[Rcpp::export]]
// List test_NS(NumericVector X, NumericVector X_test, NumericVector y, long K){
//   NS * a = new NS(X, X_test, y, K, 1);
//   Rcout << a->get_theta() << std::endl;
//   for(int i = 1; i < 10000; ++i)
//     a->update();
//   Rcout << a->get_theta() << std::endl;
//   return List::create(Named("ns_predict") = a->predict(), Named("ns_outcome") = wrap(a->ns_outcome), Named("eta") = a->get_eta(), Named("gamma") = a->get_gamma(), Named("theta") = a->get_theta());
//   //return List::create(Named("ns_outcome") = wrap(a->ns_outcome), Named("lr_outcome") = wrap(a->lr_outcome), Named("ns_part_outcome") = wrap(a->ns_part_outcome), Named("eta") = a->get_eta(), Named("gamma") = a->get_gamma(), Named("theta") = a->get_theta(), Named("boundary_knots") = a->get_boundary_knots(), Named("knots") = a->get_knots(), Named("internal_knots") = a->get_internal_knots(), Named("ns_part") = a->get_ns_part(), Named("lr_part") = a->get_lr_part(), Named("basis") = a->get_basis());
// };
