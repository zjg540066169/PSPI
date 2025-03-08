#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef CBART_H_
#define CBART_H_
#include "BARTforCausal.h"
#endif


#ifndef RCPP_H_
#define RCPP_H_
#include <Rcpp.h>
#endif


using namespace Rcpp;
/*** R
bartModelMatrix=function(X, numcut=0L, usequants=FALSE, type=7,
                         rm.const=FALSE, cont=FALSE, xinfo=NULL) {
  #print(xinfo)
  X.class = class(X)[1]
  
  if(X.class=='factor') {
    X.class='data.frame'
    X=data.frame(X=X)
  }
  grp=NULL
  if(X.class=='data.frame') {
    print(X)
    p=dim(X)[2]
    
    xnm = names(X)
    for(i in 1:p) {
      print(i)
      if(is.factor(X[[i]])) {
        #print(i)
        Xtemp = nnet::class.ind(X[[i]])
        colnames(Xtemp) = paste(xnm[i],1:ncol(Xtemp),sep='')
        X[[i]]=Xtemp
        m=ncol(Xtemp)
        grp=c(grp, rep(m, m))
      } else {
        X[[i]]=cbind(X[[i]])
        colnames(X[[i]])=xnm[i]
        grp=c(grp, 1)
        ##grp=c(grp, i)
      }
    }
    Xtemp=cbind(X[[1]])
    if(p>1) for(i in 2:p) Xtemp=cbind(Xtemp, X[[i]])
    X=Xtemp
  }
  else if(X.class=='numeric' | X.class=='integer') {
    X=cbind(as.numeric(X))
    ##grp=1
  }
  else if(X.class=='NULL') return(X)
  else if(X.class!='matrix')
    stop('Expecting either a factor, a vector, a matrix or a data.frame')
  
  #print(789)
  N <- nrow(X)
  p <- ncol(X)
  
  xinfo. <- matrix(nrow=p, ncol=numcut)
  nc <- numcut
  rm.vars <- c()
  #print(789)
  if(length(xinfo)==0 & N>0 & p>0 & (rm.const | numcut[1]>0)) {
    for(j in 1:p) {
      X.class <- class(X[1, j])[1]
      
      if(X.class=='numeric' | X.class=='integer') {
        xs <- unique(sort(X[ , j]))
        k <- length(xs)
        nc[j] <- numcut
        
        if(k %in% 0:1) { # deal with constant variables
          rm.vars <- c(rm.vars, -j)
          nc[j] <- 1
          if(k==0) xs <- NA
        }
        else if(cont)
          xs <- seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
        else if(k<numcut) {
          xs <- 0.5*(xs[1:(k-1)]+xs[2:k]) #  if k < numcut, use middle point between values to split
          nc[j] <- k-1
        }
        else if(usequants) { 
          xs <- quantile(X[ , j], type=type,
                         probs=(0:(numcut+1))/(numcut+1))[-c(1, numcut+2)]
          names(xs) <- NULL
        }
        else xs <-
          seq(xs[1], xs[k], length.out=numcut+2)[-c(1, numcut+2)]
      }
      else
        stop(paste0('Variables of type ', X.class, ' are not supported'))
      
      xinfo.[j, 1:nc[j] ] <- xs
    }
  }
  
  X <- data.matrix(X)
  #print(xinfo)
  #print(dim(xinfo.))
  if(length(xinfo)>0) {
    if(is.list(xinfo)) for(j in 1:p){
      
      xinfo.[j, 1:length(xinfo[[j]])] <- xinfo[[j]]
    } 
    else if(is.matrix(xinfo)) xinfo. <- xinfo
    else stop('Only a list or a matrix can be provided for xinfo')
    
    for(j in 1:p) nc[j] <- sum(!is.na(xinfo.[j, ]))
  }
  
  xinfo <- xinfo.
  
  if(rm.const && length(rm.vars)>0 &&
     !(length(rm.vars)==p && all((1:p)==(-rm.vars)))) {
    X <- X[ , rm.vars]
    nc <- nc[rm.vars]
    xinfo <- xinfo[rm.vars, ]
    grp <- grp[rm.vars]
  }
  else if(length(rm.vars)==0 || (length(rm.vars)==p && all((1:p)==(-rm.vars))))
    rm.vars <- 1:p
  
  dimnames(xinfo) <- list(dimnames(X)[[2]], NULL)
  
  if(all(numcut==0)) return(X)
  else return(list(X=X, numcut=as.integer(nc), rm.const=rm.vars,
                   xinfo=xinfo, grp=grp))
}
*/


class model_1: public BARTforCausal{
public:
  model_1(NumericMatrix X_, NumericVector Y_, NumericVector Z_, NumericVector pi_, NumericMatrix X_test_, bool binary, long ntrees_s) : BARTforCausal(X_, Y_, Z_, pi_, X_test_, binary, ntrees_s){
    Z_1 = (Z == 1.0);
    //main_bart = new bart_model(sliceRows(cbind(X, pi), !Z_1), Y[!Z_1], 100L, false, false, false, 200);
    main_bart = new bart_model(cbind(X, pi), Y, 100L, false, false, false, 200);
    //main_bart = new bart_model(cbind(X, pi), Y, 100L, false, false, false, 200);
    main_bart->update(100, 100, 1, false, 10L);
    if(!this->binary)
      sigma = main_bart->get_sigma();
    else
      sigma = 1;
    bart_pre = colMeans(main_bart->predict(cbind(this->X, this->pi)));
    //main_bart_mean = mean(bart_pre);
    //bart_pre = bart_pre - main_bart_mean;
    
    //Z_1 = (Z == 1.0);
    X_Z = sliceRows(cbind(X, pi), Z_1);
    Y_Z = Y[Z_1] - bart_pre[Z_1];
    cbart = new bart_model(X_Z, Y_Z, 100L, false, false, false, ntrees_s);
    cbart->update(sigma, 50, 50, 1, false, 10L);
    Z_cbart = NumericVector(Y.length());
    this->update_Z_cbart();
  };
  
  void update_Z_cbart(){
    cbart_pre = colMeans(cbart->predict(X_Z));
    Z_cbart[Z_1] = cbart_pre;
  }
  
  void update(bool verbose = false) override{
    main_bart->set_data(cbind(X, pi), Y - Z_cbart);
    main_bart->update(sigma, 1, 1, 1, false, 10L);
    
    bart_pre = colMeans(main_bart->predict(cbind(X, pi)));
    Y_Z = Y[Z_1] - bart_pre[Z_1];
    cbart->set_data(X_Z, Y_Z);
    cbart->update(sigma, 1, 1, 1, false, 10L);
    this->update_Z_cbart();
    if(!this->binary){
      double rss = sum(pow(Y - Z_cbart - bart_pre, 2));
      sigma = main_bart->get_invchi(n, rss);
      if (verbose)
        Rcout << rss << "  " << sigma << std::endl;
    }else{
      for(int i = 0; i < n; ++i){
        if(Y[i] < 0){
          NumericVector mean_y = rtruncnorm(1, bart_pre[i] + Z_cbart[i], sigma, R_NegInf, 0);
          Y[i] = mean_y[0];
        }else{
          NumericVector mean_y = rtruncnorm(1, bart_pre[i] + Z_cbart[i], sigma, 0, R_PosInf);
          Y[i] = mean_y[0];
        }
      }
    }
  };
  
  List predict(NumericVector pi_test) override{
    long N = X_test.nrow();
    NumericVector Z_1 (N, 1);
    NumericVector Z_0 (N, 0);
    NumericVector outcome_0 = colMeans(main_bart->predict(cbind(X_test, pi_test)));// - main_bart_mean;
    NumericVector outcome_1 = outcome_0 + colMeans(cbart->predict(cbind(X_test, pi_test)));
    NumericVector outcome_0_hidden(N);
    NumericVector outcome_1_hidden(N);
    if(this->binary){
      for(int i = 0; i < N; ++i){
        outcome_1_hidden[i] = outcome_1[i];
        outcome_0_hidden[i] = outcome_0[i];
        
        outcome_1[i] = R::rbinom(1, R::pnorm(outcome_1_hidden[i], 0, 1, true, false));
        outcome_0[i] = R::rbinom(1, R::pnorm(outcome_0_hidden[i], 0, 1, true, false));
      }
    }else{
      for(int i = 0; i < N; ++i){
        outcome_1[i] = outcome_1[i] + R::rnorm(0, sigma);
        outcome_0[i] = outcome_0[i] + R::rnorm(0, sigma);
      }
    }
    return List::create(Named("outcome_1") = outcome_1, Named("outcome_0") = outcome_0, Named("outcome_1_hidden") = outcome_1_hidden, Named("outcome_0_hidden") = outcome_0_hidden);
    
  };
  
  List get_posterior() override{
    return List::create(
      Named("sigma") = sigma,
      Named("bart_pre") = bart_pre,
      Named("cbart_pre") = cbart_pre,
      Named("Z_cbart") = Z_cbart,
      Named("Y_hat") = Y
    );
  };
  
private:
  double sigma;
  NumericVector bart_pre;
  bart_model * cbart;
  LogicalVector Z_1;
  NumericMatrix X_Z;
  NumericVector Y_Z;
  NumericVector Z_cbart;
  NumericVector cbart_pre;
  
  //double main_bart_mean;
};