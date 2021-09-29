// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "settings.h"
#include <iostream>

#include <stdio.h>

#include "parameters.h"		/* Read parameters interactively */
#include "mclasses/kernel_functions.h"
#include "mclasses/matrix_h2.h"
#include "mclasses/matrix_h.h"


/*
 * 
 */
//[[Rcpp::export]]
void cholesky_hmatrix(Rcpp::XPtr<matrix_h> gp){
    gp->cholesky(); 
}


//[[Rcpp::export]]
NumericMatrix solve_hmatrix(Rcpp::XPtr<matrix_h> gp,
                            arma::vec x){
    return wrap( gp->triang_invmul_vec(x));  
}
//[[Rcpp::export]]
double log_sqrt_determinant_hmatrix(Rcpp::XPtr<matrix_h> gp){
   return(gp->log_determinant());
}

//[[Rcpp::export]]
void save_h_matrix(string file, Rcpp::XPtr<matrix_h> gp){
  gp->save_h_netcdf(file);
}

//[[Rcpp::export]]
Rcpp::XPtr<matrix_h> read_h_matrix(string file){
  matrix_h *hmat  = new matrix_h(file); 
  Rcpp:XPtr<matrix_h> rv(hmat,true); 
  return rv; 
}

//[[Rcpp::export]]
arma::vec multiply_hmat_vec(Rcpp::XPtr<matrix_h> h, arma::vec v){
  return h->multiply_hmat_vec(v);
}