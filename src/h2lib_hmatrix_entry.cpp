// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <settings.h>
#include <iostream>
#include <stdio.h>
#include <parameters.h>		/* Read parameters interactively */
#include <mclasses/kernel_functions.h>
#include <mclasses/matrix_h2.h>
#include <mclasses/matrix_h.h>

/*
 *  Set of Rcpp import wrappers for functions from
 *  H2lib - some are called directly from R functions,
 *  but most are setup to be called from other C++ functions
 *  custom-written as part of the source.
 *
 *  Pulls from three header files in /src/include/mclasses/:
 *  * kernel_functions
 *  * matrix_h2
 *  * matrix_h
 *
 *  Creates exports for:
 *  * cholesky_hmatrix
 *  * solve_hmatrix
 *  * log_sqrt_determinant
 *  * save_h_matrix
 *  * read_h_matrix
 *  * multiply_hmat_vec
 */


//******************************************************************************

//' cholesky_hmatrix
//'
//' Convert a matrix into its Cholesky form.
//'
//' @param gp An H-matrix, likely created by create_h2_kernel or equivalent.
//' @return Cholesky decomposition of input.
//' @export
//' @examples
//'  A <- matrix(data = rnorm(1000),
//'              nrow = 500,
//'              ncol = 2)
//'  gp_h2 <- create_h2_kernel(X = as.matrix(A),
//'                            hp = as.matrix(c(16, 1), ncol = 1),
//'                            kt = 1)
//'  gp_h  = convert_h2_to_h_matrix(gp = gp_h2)
//'  cholesky_hmatrix(gp = gp_h)
//'
//[[Rcpp::export]]
void cholesky_hmatrix(Rcpp::XPtr<matrix_h> gp){
    gp->cholesky();
}

//******************************************************************************

//' solve_hmatrix
//'
//' Solve Ax=b, for A an H-matrix
//'
//' Solve a matrix equation Ax=b for x, where A is approximated in H-matrix form.
//'
//' @param gp An H-matrix, likely created by create_h2_kernel() or equivalent.
//' @param b A vector of constants.
//' @return Solution of matrix operation, x.
//' @export
//' @examples
//'  A <- matrix(data = rnorm(1000),
//'              nrow = 500,
//'              ncol = 2)
//'  gp_h2 <- create_h2_kernel(X = as.matrix(A),
//'                            hp = as.matrix(c(16, 1), ncol = 1),
//'                            kt = 1)
//'  gp_h  = convert_h2_to_h_matrix(gp_h2)
//'  cholesky_hmatrix(gp_h)
//'  b <- matrix(data = rnorm(500), nrow = 500, ncol = 1)
//'  x <- solve_hmatrix(gp = gp_h, b = b)
//'
//[[Rcpp::export]]
NumericMatrix solve_hmatrix(Rcpp::XPtr<matrix_h> gp,
                            arma::vec b){
    return wrap( gp -> triang_invmul_vec(b));
}

//******************************************************************************

//' log_sqrt_determinant_hmatrix
//'
//' Compute the log determinant of an H-matrix
//'
//' Compute the log-determinant of the input H-matrix. Note: no square root
//' is taken at any point in this process, so name is confusing.
//'
//' @param gp An H-matrix, likely created by create_h2_kernel() or equivalent.
//' @export
//' @return Log determinant of the input H-matrix.
//'
//[[Rcpp::export]]
double log_sqrt_determinant_hmatrix(Rcpp::XPtr<matrix_h> gp){
   return(gp -> log_determinant());
}

//******************************************************************************

//' save_h_matrix
//'
//' Save an H-matrix to disk in NetCDF binary format.
//'
//' @param file System file reference or file name.
//' @param gp An H-matrix, likely created by create_h2_kernel or equivalent.
//' @export
//' @return Nothing.
//'
//[[Rcpp::export]]
void save_h_matrix(string file, Rcpp::XPtr<matrix_h> gp){
  gp->save_h_netcdf(file);
}

//******************************************************************************

//' read_h_matrix
//'
//' Read an H-matrix from disk in NetCDF binary format.
//'
//' Read an H-matrix from disk in NetCDF binary format. Partner function
//' to save_h_matrix().
//'
//' @param file System file reference or file name.
//' @export
//' @return rv An H-matrix.
//' @examples
//'
//'
//'
//[[Rcpp::export]]
Rcpp::XPtr<matrix_h> read_h_matrix(string file){
  matrix_h *hmat  = new matrix_h(file);
  Rcpp:XPtr<matrix_h> rv(hmat, true);
  return rv;
}

//******************************************************************************

//' multiply_hmat_vec
//'
//' Multiply an H-matrix by a vector.
//'
//' @param h An H-matrix.
//' @param v A vector.
//' @export
//' @return The matrix product of h*v
//'
//[[Rcpp::export]]
arma::vec multiply_hmat_vec(Rcpp::XPtr<matrix_h> h, arma::vec v){
  return h -> multiply_hmat_vec(v);
}
