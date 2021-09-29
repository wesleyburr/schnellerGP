// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "settings.h"
#include <iostream>
#include <string>
#include <stdio.h>

#include "parameters.h"		/* Read parameters interactively */
#include "mclasses/kernel_functions.h"
#include "mclasses/matrix_h2.h"
#include "mclasses/matrix_h.h"


/*greate_gp_h2lib
 *  @X        - Matrix of X coordinates to make the GP kernel
 *  @gp_parms - the parameters for the gaussian process relative to the
 *              current chosen kernel
 *  @kernel_type - the type of kernel current values are in
 *                  the enum type kernel_types
 *
 *  Purpose: Return a pointer to a new instance of a h2_lib object.
 *
 */
//[[Rcpp::export]]
Rcpp::XPtr<matrix_h2>  create_h2_kernel(arma::mat X,
                                        arma::mat hp,
                                        NumericVector kt){
  
  matrix_h2 *gp;
  kernel_type ktype = (kernel_type)kt[0];
  
  switch(ktype){
    case kernel_type::gaussian:
      
      gp = new matrix_h2(X, kernel_sqr_exp, hp,75,2,1e-8);
     // gp->draw_matrix(); 
      break;
    case kernel_type::gaussian_truncate:   
      gp = new matrix_h2(X, kernel_sqr_exp_truncate, hp,50,2,1e-8);
      break; 
    case kernel_type::matern:
      gp = new matrix_h2(X, kernel_matern_exp, hp,50,2,1e-8);
      //stop("The matern kernel has not yet been implemented.\n");
      break;
    case kernel_type::diagonal:
      gp = new matrix_h2(X,kernel_diagonal, hp);
      break;
    default:
      stop("The requested kernel is unknown.\n");
    }
  Rcpp::XPtr<matrix_h2> rv(gp,true);
  return rv;
}

//[[Rcpp::export]]
NumericMatrix  solve_h2_vector(Rcpp::XPtr<matrix_h2> gp,
                               arma::vec d_precond,
                               arma::mat X,
                               double tolerance){
  diagonal_precond pc;
  pc.n = d_precond.n_rows;
  pc.d = new double[pc.n];
  for (uint i = 0; i < pc.n; i++){
      pc.d[i] = d_precond[i];
  }
  arma::mat rv = gp->solve_precond_vect(X,&pc,tolerance,5000);
  delete (pc.d);
  return wrap(rv);
}

//[[Rcpp::export]]
NumericMatrix  solve_h2_vector_start(Rcpp::XPtr<matrix_h2> gp,
                               arma::vec d_precond,
                               arma::mat X,
                               arma::mat start, 
                               double tolerance){
  diagonal_precond pc;
  pc.n = d_precond.n_rows;
  pc.d = new double[pc.n];
  for (uint i = 0; i < pc.n; i++){
    pc.d[i] = d_precond[i];
  }
  arma::mat rv = gp->solve_precond_vect_wstart(X,&pc,tolerance,5000,start);
  delete (pc.d);
  return wrap(rv);
}

//[[Rcpp::export]]
void add_diag(Rcpp::XPtr<matrix_h2> gp,
              arma::vec dv){
  //FIX ME: change to avector not arma:vec!
  pavector pvec = new_avector (dv.n_elem); 	
  for (uint i = 0; i < pvec->dim; i++ ){
    pvec->v[i] = dv[i]; 
  }
  gp->add_diag_avector(pvec);
  del_avector(pvec); 
  
}

//[[Rcpp::export]]
NumericMatrix  multiply_h2_vector(Rcpp::XPtr<matrix_h2> gp,
                                  arma::vec dv){
  
  return wrap(gp->multiply_vector(dv));
   
}

//[[Rcpp::export]]
NumericMatrix  multiply_h2_stdnorm(Rcpp::XPtr<matrix_h2> gp){
  //FIX ME: SET SEED
  arma::vec dv = arma::randn<arma::vec>(gp->getdim()); 
  return wrap(gp->multiply_vector(dv));
  
}

//
//
//[[Rcpp::export]]
Rcpp::XPtr<matrix_h> convert_h2_to_h_matrix(Rcpp::XPtr<matrix_h2> gp){
    phmatrix h = gp->convert_to_hmatrix();
    matrix_h *hclass = new matrix_h(h);
    Rcpp:XPtr<matrix_h>  rv(hclass,true);
    del_hmatrix(h);
    return rv;
}



//[[Rcpp::export]]
void save_h2_matrix(string file,Rcpp::XPtr<matrix_h2> gp){
  gp->save_h2_netcdf(file); 
}

//[[Rcpp::export]]
arma::umat get_h2_idx(Rcpp::XPtr<matrix_h2> gp){
  return gp->getIndex()+1;   
}

