#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef DIST_H_
#define DIST_H_
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
#endif


#ifndef BART_H_
#define BART_H_
#include "bart_model.h"
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


class BARTforCausal{
public:
  BARTforCausal(NumericMatrix X_, NumericVector Y_, NumericVector Z_, NumericVector pi_, bool binary, long ntrees_s = 200){
    X = clone(X_);
    Y = Y_;
    Z = clone(Z_);
    pi = clone(pi_);
    main_bart = NULL;
    this->binary = binary;
    n = Y.length();
    this->ntrees_s = ntrees_s;
    //Rcout << "father initial";
  }
  
  BARTforCausal(NumericMatrix X_, LogicalVector Y_, NumericVector Z_, NumericVector pi_, bool binary, long ntrees_s = 200){
    X = clone(X_);
    Y = as<NumericVector>(clone(Y_)) * 2 - 1;
    Z = clone(Z_);
    pi = clone(pi_);
    main_bart = NULL;
    this->binary = binary;
    n = Y.length();
    this->ntrees_s = ntrees_s;
    //Rcout << "father initial";
  }
  
  
  virtual void update(bool verbose = false) = 0;
  virtual List predict(NumericMatrix X_test, NumericVector pi_test) = 0;
  virtual List get_posterior() = 0;
  
  NumericMatrix get_X(){
    return X;
  };
  NumericVector get_Y(){
    return Y;
  };
  NumericVector get_pi(){
    return pi;
  };
  NumericVector get_Z(){
    return Z;
  };
  
  void update_X(NumericMatrix X_new){
    X = X_new;
  }
  
  void update_Y(NumericVector Y_new){
    Y = Y_new;
  }
  
  void update_pi(NumericVector pi_new){
    pi = pi_new;
  }
  
  void update_pi(LogicalVector Z_new){
    Z = Z_new;
  }
  
  NumericMatrix cbind(NumericMatrix a, NumericVector b){
    int row = a.nrow();
    int col = a.ncol();
    
    NumericMatrix result(row, col+1);
    for(int i = 0 ; i < col; ++i){
      result(_, i) = a(_, i);
    }
    result(_, col) = b;
    return result;
  }
  
  NumericMatrix cbind(NumericMatrix a, NumericVector b, NumericVector c){
    int row = a.nrow();
    int col = a.ncol();
    
    NumericMatrix result(row, col+2);
    for(int i = 0 ; i < col; ++i){
      result(_, i) = a(_, i);
    }
    result(_, col) = b;
    result(_, col + 1) = c;
    return result;
  }
  

  NumericMatrix sliceRows(NumericMatrix mat, LogicalVector vec) {
      int n = mat.nrow();
      int m = mat.ncol();
      int count = sum(vec);
      
      // Create a new NumericMatrix to store the result
      NumericMatrix result(count, m);
      
      // Fill the result matrix with the selected rows
      int rowIndex = 0;
      for (int i = 0; i < n; i++) {
        if (vec[i]) {
          for (int j = 0; j < m; j++) {
            result(rowIndex, j) = mat(i, j);
          }
          rowIndex++;
        }
      }
      return result;
  }
  
protected:
  NumericMatrix X;
  NumericVector Y;
  NumericVector Z;
  NumericVector pi;
  bart_model * main_bart;
  bool binary;
  long n;
  long ntrees_s;
};