#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef RCPP_H_
#define RCPP_H_
#include <Rcpp.h>
#endif

#ifndef RCPPDIST_H_
#define RCPPDIST_H_
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#endif

using namespace Rcpp;
using namespace arma;



class regression_model{
public:
  regression_model(){};
  
  regression_model(NumericMatrix x_train, NumericVector y_train, bool has_intecept = true){
    this->has_intecept = has_intecept;
    this->n = y_train.length();
    this->p = x_train.ncol();
    if(!has_intecept){
      this->p++;
      this->X = NumericMatrix(this->n, this->p);
      this->X(_, 0) = this->X(_, 0) + 1;
      for(int i = 0 ; i < x_train.ncol(); ++i){
        this->X(_, i+1) = x_train(_, i);
      }
    }else{
      this->X = NumericMatrix(this->n, this->p);
      this->X(_, 0) = this->X(_, 0) + 1;
      for(int i = 0 ; i < x_train.ncol(); ++i){
        this->X(_, i) = x_train(_, i);
      }
    }
    this->Y = clone(y_train);
    beta = update_beta_mean(this->X, this->Y);
    //sigma = 2.173695;
    sigma = initial_sigma(this->X, this->Y, this->beta);
    //Rcout << beta << std::endl;
    //Rcout << sigma << std::endl;
  };
  
  void update(){
    this->beta = update_beta(this->X, this->Y, this->sigma);
    this->sigma = sqrt(update_sigma2(this->X, this->Y, this->beta));
  };
  
  void update(double sigma){
    this->sigma = sigma;
    this->beta = update_beta(this->X, this->Y, this->sigma);
  };
  
  void set_data(NumericMatrix x_train, NumericVector y_train, bool has_intecept = true){
    this->has_intecept = has_intecept;
    this->n = y_train.length();
    this->p = x_train.ncol();
    if(!has_intecept){
      this->p++;
      this->X = NumericMatrix(this->n, this->p);
      this->X(_, 0) = this->X(_, 0) + 1;
      for(int i = 0 ; i < x_train.ncol(); ++i){
        this->X(_, i+1) = x_train(_, i);
      }
    }else{
      this->X = NumericMatrix(this->n, this->p);
      this->X(_, 0) = this->X(_, 0) + 1;
      for(int i = 0 ; i < x_train.ncol(); ++i){
        this->X(_, i) = x_train(_, i);
      }
    }
    this->Y = clone(y_train);
  };
  
  NumericVector predict(NumericMatrix x_predict, bool has_intecept = true){
    NumericMatrix X_pre;
    if(!has_intecept){
      X_pre = NumericMatrix(x_predict.nrow(), x_predict.ncol() + 1);
      X_pre(_, 0) = X_pre(_, 0) + 1;
      for(int i = 0 ; i < x_predict.ncol(); ++i){
        X_pre(_, i+1) = x_predict(_, i);
      }
    }else{
      X_pre = clone(x_predict);
    }
    vec x_beta = as<mat>(X_pre) * as<vec>(beta);
    NumericVector y_predict(x_predict.nrow());
    for(int i = 0 ; i < x_predict.nrow(); ++i){
      y_predict[i] = R::rnorm(x_beta[i], sigma);
    }
    return y_predict;
  };
  
  NumericMatrix inv_X_T_X(NumericMatrix x){
    mat mat_X = as<arma::mat>(x);
    return(wrap(arma::inv_sympd(mat_X.t() * mat_X)));
  }
  
  NumericVector update_beta_mean(NumericMatrix x, NumericVector y){
    mat mat_X = as<arma::mat>(x);
    vec mat_Y = as<arma::vec>(y);
    return(wrap(arma::inv_sympd(mat_X.t() * mat_X) * mat_X.t() * mat_Y));
  }
  
  NumericVector update_beta(NumericMatrix x, NumericVector y, double sigma){
    mat mat_X = as<arma::mat>(x);
    vec mat_Y = as<arma::vec>(y);
    vec beta_mean = as<arma::vec>(update_beta_mean(x, y));
    mat beta_var = as<arma::mat>(inv_X_T_X(x)) * sigma * sigma;
    //Rcout << beta_var << std::endl;
    return wrap(rmvnorm(1, beta_mean, beta_var));
  }
  
  double update_sigma2(NumericMatrix x, NumericVector y, NumericVector beta_){
    mat mat_X = as<arma::mat>(x);
    vec mat_Y = as<arma::vec>(y);
    vec mat_beta = as<arma::vec>(beta_);
    double a = (y.length()) / 2;
    mat b = ((mat_Y - mat_X * mat_beta).t() * (mat_Y - mat_X * mat_beta)) / 2;
    return rinvgamma(a, b(0, 0));
  }
  
  double initial_sigma(NumericMatrix x, NumericVector y, NumericVector beta_){
    mat mat_X = as<arma::mat>(x);
    vec mat_Y = as<arma::vec>(y);
    vec mat_beta = as<arma::vec>(beta_);
    double a = y.length();
    mat b = ((mat_Y - mat_X * mat_beta).t() * (mat_Y - mat_X * mat_beta)) / 2;
    //Rcout << b << std::endl;
    return sqrt(b(0, 0) / (a - 1));
  }
  
  NumericVector get_beta(){
    return this->beta;
  }
  
  double get_sigma(){
    return this->sigma;
  }
  
  double rinvgamma(double a, double b){
    double s = R::rgamma(a, 1 / b);
    return 1 / s;
  }
 
private:
  bool has_intecept;
  
  long n;
  long p;
  NumericMatrix X;
  NumericVector Y;
  
  NumericVector beta;
  double sigma;
};




// [[Rcpp::export]]
List test_regression(NumericMatrix X, NumericVector Y, long nburn, long npost, NumericMatrix X_test){
  regression_model * a = new regression_model(X, Y, false);
  a->update();
  NumericMatrix beta_n(npost, a->get_beta().length());
  NumericVector sigma_n(npost);
  NumericMatrix y_predict(npost, X_test.nrow());
  for(int i = 0 ; i < nburn + npost; ++i){
    a->update();
    if(i >= nburn){
      beta_n(i - nburn, _) = a->get_beta();
      sigma_n[i - nburn] = a->get_sigma();
      y_predict(i - nburn, _) = a->predict(X_test, false);
    }
    
  }
  return List::create(Named("beta_n") = beta_n, Named("sigma_n") = sigma_n, Named("y_predict") = y_predict);
};
