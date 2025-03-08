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
#include "NS.h"
#endif

#ifndef RCPPDIST_H_
#define RCPPDIST_H_
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#endif

#include <algorithm>

using namespace Rcpp;
using namespace arma;



class NS_HC: public NS{
public:
  NS_HC(){};
  
  NS_HC(NumericVector X, NumericVector y, long K, double sigma, int order = 3, bool horseshoe = false, double A_lambda = 1.0, double A_gamma = 1.5): NS(X, y, K, sigma, order){
    // this->X = X;
    // this->y = as<arma::vec>(y);
    // this->K = K;
    // this->n = X.length();
    // this->sigma = sigma;
    this->A_lambda = A_lambda;
    this->A_gamma = A_gamma;
    this->horseshoe = horseshoe;
    
    eta = rinvgamma(0.5, 1 / pow(A_gamma , 2));
    gamma = sqrt(rinvgamma(0.5, 1 / eta));
    
    if(horseshoe){
      for(int k = 2; k <= K-1; ++k){
        double v = rinvgamma(0.5, 1 / pow(A_lambda, 2));
        vk.insert_rows(vk.n_rows, 1);
        vk(vk.n_rows - 1) = v;
        
        lambda_k.insert_rows(lambda_k.n_rows, 1);
        lambda_k(lambda_k.n_rows - 1) = sqrt(rinvgamma(0.5, 1 / v));
      }
    }
    
      
    update_outcome();
  };
  
  void update(){
    if(K > 2){
      lr_coefficient = update_beta(lr_basis, this->y - ns_part_outcome, this->sigma);
    }else{
      lr_coefficient = update_beta(lr_basis, this->y, this->sigma);
    }
    lr_outcome = lr_basis * lr_coefficient;
    if(K > 2){
      if(horseshoe){
        arma::mat Lambda = arma::diagmat(pow(lambda_k, 2)) * pow(gamma, 2);
        arma::mat theta_var = arma::inv_sympd(arma::inv_sympd(Lambda) + 1 / pow(sigma, 2) * ns_basis.t() * ns_basis);
        arma::vec theta_mean = arma::inv_sympd(pow(sigma, 2) * arma::inv_sympd(Lambda) + ns_basis.t() * ns_basis) * ns_basis.t() * (this->y - lr_outcome);
        ns_coefficient = arma::vectorise(rmvnorm(1, theta_mean, theta_var));
      }else{
        arma::mat theta_var = arma::inv_sympd(1 / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + 1 / pow(this->sigma, 2) * ns_basis.t() * ns_basis);
        arma::vec theta_mean = arma::inv_sympd(pow(this->sigma, 2) / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + ns_basis.t() * ns_basis) * ns_basis.t() * (this->y - lr_outcome);
        ns_coefficient = arma::vectorise(rmvnorm(1, theta_mean, theta_var));
      }
      
      ns_part_outcome = ns_basis * ns_coefficient;
      ns_outcome = lr_outcome + ns_part_outcome;
      if(horseshoe){
        for(int k = 0; k < K-2; ++k){
          lambda_k[k] = sqrt(rinvgamma(1.0, 1.0 / vk[k] + 0.5 * pow(ns_coefficient[k] / gamma, 2)));
          vk[k] = rinvgamma(0.5, (1.0 / (A_lambda*A_lambda)) + 1.0 / (lambda_k[k] * lambda_k[k]));
        }
        gamma = sqrt(rinvgamma((K-1.0) / 2.0, 1 / eta + 0.5 * arma::dot(ns_coefficient / lambda_k, ns_coefficient / lambda_k)));
      }else{
        gamma = sqrt(rinvgamma((K-1.0) / 2.0, 1 / eta + 0.5 * arma::dot(ns_coefficient, ns_coefficient)));
      }
      eta = rinvgamma(0.5, (1.0/(A_gamma*A_gamma)) + 1.0 / (gamma * gamma));
    }
  }
  
  void update(double sigma){
    this->sigma = sigma;
    if(K > 2){
      lr_coefficient = update_beta(lr_basis, this->y - ns_part_outcome, this->sigma);
    }else{
      lr_coefficient = update_beta(lr_basis, this->y, this->sigma);
    }
    lr_outcome = lr_basis * lr_coefficient;
    if(K > 2){
      if(horseshoe){
        arma::mat Lambda = arma::diagmat(pow(lambda_k, 2)) * pow(gamma, 2);
        arma::mat theta_var = arma::inv_sympd(arma::inv_sympd(Lambda) + 1 / pow(sigma, 2) * ns_basis.t() * ns_basis);
        arma::vec theta_mean = arma::inv_sympd(pow(sigma, 2) * arma::inv_sympd(Lambda) + ns_basis.t() * ns_basis) * ns_basis.t() * (this->y - lr_outcome);
        ns_coefficient = arma::vectorise(rmvnorm(1, theta_mean, theta_var));
      }else{
        arma::mat theta_var = arma::inv_sympd(1 / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + 1 / pow(this->sigma, 2) * ns_basis.t() * ns_basis);
        arma::vec theta_mean = arma::inv_sympd(pow(this->sigma, 2) / pow(gamma, 2) * arma::eye<arma::mat>(K - 2, K - 2) + ns_basis.t() * ns_basis) * ns_basis.t() * (this->y - lr_outcome);
        ns_coefficient = arma::vectorise(rmvnorm(1, theta_mean, theta_var));
      }
      
      
      ns_part_outcome = ns_basis * ns_coefficient;
      ns_outcome = lr_outcome + ns_part_outcome;
      Rcout << lr_coefficient.t() << " " << ns_coefficient.t() << std::endl;
      
      
      if(horseshoe){
        for(int k = 0; k < K-2; ++k){
          lambda_k[k] = sqrt(rinvgamma(1.0, 1.0 / vk[k] + 0.5 * pow(ns_coefficient[k] / gamma, 2)));
          vk[k] = rinvgamma(0.5, (1.0 / (A_lambda*A_lambda)) + 1.0 / (lambda_k[k] * lambda_k[k]));
        }
        gamma = sqrt(rinvgamma((K-1.0) / 2.0, 1 / eta + 0.5 * arma::dot(ns_coefficient / lambda_k, ns_coefficient / lambda_k)));
      }else{
        gamma = sqrt(rinvgamma((K-1.0) / 2.0, 1 / eta + 0.5 * arma::dot(ns_coefficient, ns_coefficient)));
      }
      eta = rinvgamma(0.5, (1.0/(A_gamma*A_gamma)) + 1.0 / (gamma * gamma));
    }else{
      Rcout << lr_coefficient.t() << std::endl;
    }
    
    
  }
  
  double get_eta(){
    return eta;
  }
  
  
  
  
protected:

  double eta;
  
  arma::vec lambda_k;
  arma::vec vk;
  
  double A_lambda;
  double A_gamma;
  
  bool horseshoe;

};


// // 
// // 
// // [[Rcpp::export]]
// List test_NS(NumericVector X, NumericVector X_test, NumericVector y, long K, bool horseshoe = false){
//   NS_HC * a = new NS_HC(X, y, K, 1, 3, horseshoe);
//   for(int i = 1; i < 10000; ++i)
//     a->update();
//   NumericVector y_test = a->predict(X_test);
//   return List::create(Named("Y_test") = y_test, Named("ns_predict") = a->predict(X_test), Named("eta") = a->get_eta(), Named("gamma") = a->get_gamma(), Named("theta") = a->get_theta());
//   //return List::create(Named("ns_outcome") = wrap(a->ns_outcome), Named("lr_outcome") = wrap(a->lr_outcome), Named("ns_part_outcome") = wrap(a->ns_part_outcome), Named("eta") = a->get_eta(), Named("gamma") = a->get_gamma(), Named("theta") = a->get_theta(), Named("boundary_knots") = a->get_boundary_knots(), Named("knots") = a->get_knots(), Named("internal_knots") = a->get_internal_knots(), Named("ns_part") = a->get_ns_part(), Named("lr_part") = a->get_lr_part(), Named("basis") = a->get_basis());
// };
