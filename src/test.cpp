#ifndef ARMADILLO_H_
#define ARMADILLO_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#endif

#ifndef SPLINE_ARMADILLO_H_
#define SPLINE_ARMADILLO_H_
#include <splines2Armadillo.h>
// [[Rcpp::depends(splines2)]]
#endif

using namespace Rcpp;
using namespace splines2;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



// [[Rcpp::export]]
NumericMatrix timesTwo(NumericVector x, NumericVector y) {
  BSpline * bsp_obj = new BSpline(x, 10);
  arma::vec boundary_knots = bsp_obj->get_knot_sequence();
  bsp_obj->set_x(y);
  //bsp_obj->set_knot_sequence(boundary_knots);
  return(wrap(bsp_obj->basis(false)));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(1:1000, 1000:2000)
*/
