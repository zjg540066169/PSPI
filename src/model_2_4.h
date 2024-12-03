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


class model_2_4: public BARTforCausal{
public:
  model_2_4(NumericMatrix X_, NumericVector Y_, NumericVector Z_, NumericVector pi_, bool binary) : BARTforCausal(X_, Y_, Z_, pi_, binary){
    main_bart = new bart_model(cbind(X, pi), Y, 100L, false, false, false, 200);
    main_bart->update(50, 50, 1, false, 10L);
    sigma = main_bart->get_sigma();
    bart_pre = colMeans(main_bart->predict(cbind(this->X, this->pi)));
    
    Z_1 = (Z == 1.0);
    X_Z = sliceRows(X, Z_1);
    Y_Z = Y[Z_1] - bart_pre[Z_1];
    pi_Z = NumericMatrix(Y_Z.length(), 1, as<NumericVector>(pi[Z_1]).begin());
    cbart_pi = new bart_model(pi_Z, Y_Z, 100L, false, false, false, 100);
    cbart_pi->update(50, 50, 1, false, 10L);
    cbart_pi_pre = colMeans(cbart_pi->predict(pi_Z));
    //Rcout << cbart_pi_pre << std::endl;
    cbart = new bart_model(X_Z, Y_Z - cbart_pi_pre, 100L, false, false, false, 100);
    cbart->update(50, 50, 1, false, 10L);
    cbart_pre = colMeans(cbart->predict(X_Z));
    cbart_pre_mean = mean(cbart_pre);
    cbart_pre = cbart_pre - cbart_pre_mean;
    
    Z_cbart = NumericVector(Y.length());
    this->update_Z_cbart();
  };
  
  void update_Z_cbart(){
    //Rcout << sum(cbart_pre) << std::endl;
    NumericVector Z_cbart_Z_1 = cbart_pre + cbart_pi_pre;
    //Rcout << Z_cbart_Z_1 << std::endl;
    Z_cbart[Z_1] = Z_cbart_Z_1;
    //Rcout << Z_cbart << std::endl;
  }
  
  void update(bool verbose = false) override{
    main_bart->set_data(cbind(X, pi), Y - Z_cbart);
    main_bart->update(1, 1, 1, verbose, 10L);
    sigma = main_bart->get_sigma();
    bart_pre = colMeans(main_bart->predict(cbind(X, pi)));
    
    Y_Z = Y[Z_1] - bart_pre[Z_1];
    NumericVector Y_Cb = Y_Z - cbart_pi_pre;

    cbart->set_data(X_Z, Y_Z - cbart_pi_pre);
    cbart->update(1, 1, 1, false, 10L);
    cbart_pre = colMeans(cbart->predict(X_Z));
    cbart_pre_mean = mean(cbart_pre);
    cbart_pre = cbart_pre - cbart_pre_mean;
    // 
    // 
    NumericVector y_te = Y_Z - cbart_pre;
    cbart_pi->set_data(pi_Z, Y_Z - cbart_pre);
    cbart_pi->update(1, 1, 1, false, 10L);
    cbart_pi_pre = colMeans(cbart_pi->predict(pi_Z));
    //Rcout << cbart_pi_pre << std::endl;
    this->update_Z_cbart();
  };
  
  List predict(NumericMatrix X_test, NumericVector pi_test) override{
    long n = X_test.nrow();
    NumericMatrix pi_test_Z = NumericMatrix(n, 1, pi_test.begin());
    NumericVector outcome_0 = colMeans(main_bart->predict(cbind(X_test, pi_test)));
    NumericVector outcome_1 = outcome_0 + colMeans(cbart->predict(X_test)) - cbart_pre_mean + colMeans(cbart_pi->predict(pi_test_Z));
    return List::create(Named("outcome_1") = outcome_1, Named("outcome_0") = outcome_0);
  };
  
  List get_posterior() override{
    return List::create(
      Named("sigma") = sigma,
      Named("bart_pre") = bart_pre,
      Named("cbart_pre") = cbart_pre,
      Named("Z_cbart") = Z_cbart,
      Named("cbart_pi_pre") = cbart_pi_pre,
      Named("cbart_pre_mean") = cbart_pre_mean
    );
  };
  
private:
  double sigma;
  NumericVector bart_pre;
  bart_model * cbart;
  bart_model * cbart_pi;
  
  
  LogicalVector Z_1;
  NumericMatrix X_Z;
  NumericVector Y_Z;
  NumericMatrix pi_Z;
  NumericVector Z_cbart;
  NumericVector cbart_pre;
  double cbart_pre_mean;
  NumericVector cbart_pi_pre;
};