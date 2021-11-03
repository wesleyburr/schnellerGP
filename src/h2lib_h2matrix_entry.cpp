// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <string>

#include "settings.h"
#include "stdio.h"
#include "parameters.h"		/* Read parameters interactively */
#include "mclasses/kernel_functions.h"
#include "mclasses/matrix_h2.h"
#include "mclasses/matrix_h.h"

using namespace Rcpp;

/*
 *  Set of Rcpp import wrappers for functions from
 *  H2lib - some are called directly from R functions,
 *  but most are setup to be called from other C++ functions
 *  custom-written as part of the source. Focuses on H2 functions.
 *
 *  Pulls from three header files in /src/include/mclasses/:
 *  * kernel_functions
 *  * matrix_h2
 *  * matrix_h
 *
 *  Creates exports for:
 *  * create_h2_kernel
 *  * solve_h2_vector
 *  * solve_h2_vector_start
 *  * add_diag
 *  * get_h2_idx
 *  * multiply_h2_vector
 *  * multiply_h2_stdnorm
 *  * convert_h2_to_h_matrix
 *  * save_h2_matrix
 *  * read_h2_matrix
 *
 */

//******************************************************************************

//' create_h2_kernel
//'
//' Return a pointer to a new instance of a h2_lib object.
//'
//' @param X  Matrix of X coordinates to use to form the GP kernel.
//' @param hp Matrix of parameters for the Gaussian Process relative to the
//'             current chosen kernel.
//' @param kt Kernel Type: the type of kernel to use as a number from the list of
//' current supported kernels: gaussian, gaussian_truncate, matern and diagonal.
//' @export
//' @return An H2 matrix with the specified kernel.
//' @examples
//'  # Note the use
//'  A <- matrix(data = rnorm(1000),
//'              nrow = 500,
//'              ncol = 2)
//'  gp_h2 <- create_h2_kernel(X = as.matrix(A),
//'                            hp = as.matrix(c(16, 1), ncol = 1),
//'                            kt = 1)
//[[Rcpp::export]]
Rcpp::XPtr<matrix_h2> create_h2_kernel(arma::mat X,
                                       arma::mat hp,
                                       NumericVector kt){

  matrix_h2 *gp;
  kernel_type ktype = (kernel_type)kt[0];

  switch(ktype){
    case kernel_type::gaussian:
      gp = new matrix_h2(X, kernel_sqr_exp, hp, 75, 2, 1e-8);
     // gp->draw_matrix();
      break;
    case kernel_type::gaussian_truncate:
      gp = new matrix_h2(X, kernel_sqr_exp_truncate, hp, 50, 2, 1e-8);
      break;
    case kernel_type::matern:
      gp = new matrix_h2(X, kernel_matern_exp, hp, 50, 2, 1e-8);
      break;
    case kernel_type::diagonal:
      gp = new matrix_h2(X, kernel_diagonal, hp);
      break;
    default:
      stop("The requested kernel is unknown.\n");
    }
  Rcpp::XPtr<matrix_h2> rv(gp, true);
  return rv;
}

//******************************************************************************

//' solve_h2_vector
//'
//' Solve Ax=b using H2 matrices.
//'
//' Solve a linear system $Ax=b$ with the preconditioned generalized minimal
//' residual method from H2lib. Note that the underlying H2 class function defaults to
//' iter = 2000 and eps = 1e-14 (eps is overwritten by argument tolerance).
//'
//' @param gp  An H2-matrix
//' @param d_precond   A vector of pre-conditions.
//' @param b  A vector b, for equation $gp*x = b$, with dimension matching that of gp.
//' @param tolerance  Tolerance for solver, passed on to solve_precond_vect
//' @param iter  Maximum number of iterations for solver, passed on to solve_precond_vect
//' @export
//' @return
//'  Given inputs gp and b, solves for x in $Ax=b$, and returns it.
//' @examples
//'
//'
//'
//[[Rcpp::export]]
NumericMatrix solve_h2_vector(Rcpp::XPtr<matrix_h2> gp,
                              arma::vec d_precond,
                              arma::vec b,
                              double tolerance,
                              uint iter){
  // given d_precond as a vector of diagonal pre-conditions,
  // reform into struct as given in matrix_h2.h
  diagonal_precond pc;
  pc.n = d_precond.n_rows;
  pc.d = new double[pc.n];
  for (uint i = 0; i < pc.n; i++){
      pc.d[i] = d_precond[i];
  }

  // pass b, &pc, tolerance, and iter
  // to arguments named armas::vec tb, diagonal_precond *p, double eps and uint iter
  arma::vec rv = gp -> solve_precond_vect(b, &pc, tolerance, iter);
  delete (&pc);
  return wrap(rv);
}

//' solve_h2_vector_start
//'
//' Solve Ax=b using H2 matrices with initial conditions.
//'
//' Solve a linear system $Ax=b$ with the preconditioned generalized minimal
//' residual method from H2lib, with provided initial condition.
//' Note that the underlying H2 class function defaults to
//' iter = 2000 and eps = 1e-14 (eps is overwritten by argument tolerance).
//'
//' @param gp  An H2-matrix
//' @param d_precond   A vector of pre-conditions.
//' @param b  A vector b, for equation $gp*x = b$, with dimension matching that of gp.
//' @param start  Initial conditions for solver.
//' @param tolerance  Tolerance for solver, passed on to solve_precond_vect
//' @param iter  Maximum number of iterations for solver, passed on to solve_precond_vect
//' @export
//' @return
//'  Given inputs gp and b, solves for x in $Ax=b$, and returns it.
//' @examples
//'
//'
//'
//[[Rcpp::export]]
NumericMatrix solve_h2_vector_start(Rcpp::XPtr<matrix_h2> gp,
                                    arma::vec d_precond,
                                    arma::vec b,
                                    arma::vec start,
                                    double tolerance,
                                    uint iter){
  diagonal_precond pc;
  pc.n = d_precond.n_rows;
  pc.d = new double[pc.n];
  for (uint i = 0; i < pc.n; i++){
    pc.d[i] = d_precond[i];
  }
  arma::vec rv = gp -> solve_precond_vect_wstart(b, &pc, start, tolerance, iter);
  delete (&pc);
  return wrap(rv);
}

//******************************************************************************

//' add_diag
//'
//' Add provided vector 'dv' to diagonal of provided H2 matrix gp.
//'
//' @param gp  An H2-matrix.
//' @param dv  A vector to be added to the diagonal.
//' @export
//' @return
//'   Takes gp, and adds dv to the diagonal in place.
//' @examples
//'
//[[Rcpp::export]]
void add_diag(Rcpp::XPtr<matrix_h2> gp,
              arma::vec dv) {

  pavector pvec = new_pointer_avector((double*)dv.memptr(),dv.n_rows);
  gp -> add_diag_avector(pvec);
  del_avector(pvec);
}

//******************************************************************************

//' multiply_h2_vector
//'
//' Left-multiply provided vector 'dv' by provided H2 matrix gp.
//'
//' @param gp  An H2-matrix.
//' @param dv  A vector to by multiplied by gp, that is, gp * dv.
//' @export
//' @return
//'   Takes gp and dv and computes the multiplication gp * dv.
//' @examples
//'
//[[Rcpp::export]]
NumericMatrix multiply_h2_vector(Rcpp::XPtr<matrix_h2> gp,
                                 arma::vec dv){
  return wrap(gp -> multiply_vector(dv));
}

//******************************************************************************

//' multiply_h2_stdnorm
//'
//' Multiply an H2 matrix by Standard Normals
//'
//' Given an H2 matrix 'gp', right-multiply it by a vector 'rv' of N(0,1) random numbers.
//'
//' @param gp  An H2-matrix.
//' @export
//' @return
//'   The multiplication gp * rv. Note that the seed in the global environment
//'   is used, as the random normal generation internally to this function
//'   is just rnorm().
//' @examples
//'
//[[Rcpp::export]]
NumericMatrix multiply_h2_stdnorm(Rcpp::XPtr<matrix_h2> gp) {
  // Rcpp export automatically wraps this function in a way which fulfills the
  // logic described by the package authors:
  //
  //    ... you must call GetRNGState prior to using them and then PutRNGState
  //    afterwards. These functions (respectively) read .Random.seed and then
  //    write it out after use.
  //  (https://gallery.rcpp.org/articles/random-number-generation/)
  //
  // Basically: random numbers which are drawn from R functions work
  // as you'd expect, using the RNG stack that R maintains.

  // getdim() is an internal for H2 that returns the number of rows
  arma::vec dv = rnorm(gp -> getdim());
  return wrap(gp -> multiply_vector(dv));
}


//******************************************************************************

//' convert_h2_to_h_matrix
//'
//' Given an H2 matrix 'gp', convert it to a H matrix.
//'
//' @param gp  An H2-matrix.
//' @export
//' @return
//'   The converted form of the matrix.
//' @examples
//'
//[[Rcpp::export]]
Rcpp::XPtr<matrix_h> convert_h2_to_h_matrix(Rcpp::XPtr<matrix_h2> gp){
    phmatrix h = gp -> convert_to_hmatrix();
    matrix_h *hclass = new matrix_h(h);
    Rcpp:XPtr<matrix_h>  rv(hclass, true);
    del_hmatrix(h);
    return rv;
}

//******************************************************************************

//' save_h2_matrix
//'
//' Save an H2-matrix to disk in NetCDF binary format.
//'
//' @param file System file reference or file name.
//' @param gp An H-matrix, likely created by create_h2_kernel or equivalent.
//' @export
//' @return Nothing.
//'
//[[Rcpp::export]]
void save_h2_matrix(string file,Rcpp::XPtr<matrix_h2> gp){
  gp->save_h2_netcdf(file);
}

//******************************************************************************

//' read_h2_matrix
//'
//' Read an H2-matrix from disk in NetCDF binary format.
//'
//' Read an H2-matrix from disk in NetCDF binary format. Partner function
//' to save_h2_matrix().
//'
//' @param file System file reference or file name.
//' @export
//' @return rv An H2-matrix.
//'
//[[Rcpp::export]]
Rcpp::XPtr<matrix_h2> read_h2_matrix(string file){
  matrix_h2 *hmat  = new matrix_h2(file);
  Rcpp:XPtr<matrix_h2> rv(hmat, true);
  return rv;
}

//******************************************************************************

//' get_h2_idx
//'
//' Given an H2 matrix gp, extract the indexes.
//'
//' @param gp  An H2 matrix.
//' @export
//' @return  The indexes of the input.
//'
//[[Rcpp::export]]
arma::umat get_h2_idx(Rcpp::XPtr<matrix_h2> gp){
  return gp -> getIndex()+1;
}

