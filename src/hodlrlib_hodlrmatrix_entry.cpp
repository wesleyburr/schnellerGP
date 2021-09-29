#include <RcppEigen.h>
#include <Rcpp.h>
#include <mclasses/squaredeMat.hpp> //squared exponential kernel
#include <mclasses/maternMat.hpp> //squared exponential kernel
#include "mclasses/HODLR_Matrix.hpp"
#include "mclasses/HODLR_Tree.hpp"
#include "mclasses/squaredMat_ptau.hpp"
#include "mclasses/maternMat_ptau.hpp"
#include "mclasses/maternMat_pX.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
int get_n_levels(int N, int M) {
  int n_levels = log(N / M) / log(2);
  return n_levels;
}

// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_MaternP1(Mat X, double sigma, double rho,
                                                         double tol, int M) {
  /*
   Rcpp setup the HODLR matrix for use in a compressed Gaussian Process with sqr exp kernel.
   X      -matrix of locations
   sigma  -function std deviation
   rho    -length-scale
   tol    -specified tolerance for accuracy of calculations
   M      -Max submatrix size
   */
  bool is_sym = true;
  bool is_pd = true;
  int N =  X.rows();
  Matern_p1_Kernel* K  = new Matern_p1_Kernel(X, N, sigma, rho);
  int n_levels      = log(N / M) / log(2);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;
}


// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_sqrExp(Mat X, double sigma, double rho,
                                                       double tol, int M) {
  /*
   Rcpp setup the HODLR matrix for use in a compressed Gaussian Process with sqr exp kernel.
   X      -matrix of locations
   sigma  -function std deviation
   rho    -length-scale
   tol    -specified tolerance for accuracy of calculations
   M      -Max submatrix size
   */
  bool is_sym = true;
  bool is_pd = true;
  int N =  X.rows();
  SQRExponential_Kernel* K  = new SQRExponential_Kernel(X, N, sigma, rho);
  int n_levels      = log(N / M) / log(2);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_sqrExpP1(Mat X, double sigma, double rho,
                                                       double tol, int M) {
  /*
   Rcpp setup the HODLR matrix for use in a compressed Gaussian Process with sqr exp kernel.
   X      -matrix of locations
   sigma  -function std deviation
   rho    -length-scale
   tol    -specified tolerance for accuracy of calculations
   M      -Max submatrix size
   */
  bool is_sym = true;
  bool is_pd = true;
  int N =  X.rows();
  SQRExponential_p1_Kernel* K  = new SQRExponential_p1_Kernel(X, N, sigma, rho);
  int n_levels      = log(N / M) / log(2);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;
}



// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_Matern(Mat X, double sigma, double rho,
                                                       double tol, int M) {
  /*
   Rcpp setup the HODLR matrix for use in a compressed Gaussian Process with matern kernel.
   X      -matrix of locations
   sigma  -function std deviation
   rho    -length-scale
   tol    -specified tolerance for accuracy of calculations
   M      -Max submatrix size
   */
  bool is_sym = true;
  bool is_pd = true;
  int N =  X.rows();
  Matern_Kernel* K         = new Matern_Kernel(X, N, sigma, rho);
  int n_levels      = log(N / M) / log(2);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;
}

// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_Matern_tP(Mat X, Mat tP, double sigma, double rho,
                                                       double tol, int M) {
  /*
   Rcpp setup the HODLR matrix for use in a compressed Gaussian Process with matern kernel.
   X      -matrix of locations
   sigma  -function std deviation
   rho    -length-scale
   tol    -specified tolerance for accuracy of calculations
   M      -Max submatrix size
   */
  bool is_sym = true;
  bool is_pd = true;
  int N =  X.rows();
 // Matern_Kernel* K         = new Matern_Kernel(X, N, sigma, rho);
  Matern_pX_Kernel* K   = new Matern_pX_Kernel(X, N, sigma, rho, tP);
  int n_levels      = log(N / M) / log(2);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;
}

// [[Rcpp::export]]
NumericVector simulate_compressedMatrixGP_TP(int N, Rcpp::XPtr<HODLR_Tree> GPobj,NumericVector tP) {
  /*
   Rcpp simulate a draw from a compressed Gaussian Process.
   N     - number of obs
   GPobj - H matrix, already assembled and factorized
   Return simulated draw from a compressed GP.
   */

  NumericVector tW = rnorm(N, 0, 1);
  tW = tW*sqrt(tP);
  Eigen::Map<Eigen::MatrixXd> ttW(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(tW));
  Eigen::MatrixXd W = ttW;

  Eigen::MatrixXd returnV = GPobj->symmetricFactorProduct(W);

  return wrap(returnV);
}


// [[Rcpp::export]]
NumericVector simulate_compressedMatrixGP(int N, Rcpp::XPtr<HODLR_Tree> GPobj) {
  /*
   Rcpp simulate a draw from a compressed Gaussian Process.
   N     - number of obs
   GPobj - H matrix, already assembled and factorized
   Return simulated draw from a compressed GP.
   */

  NumericVector tW = rnorm(N, 0, 1);
  Eigen::Map<Eigen::MatrixXd> ttW(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(tW));
  Eigen::MatrixXd W = ttW;

  Eigen::MatrixXd returnV = GPobj->symmetricFactorProduct(W);

  return wrap(returnV);
}

// [[Rcpp::export]]
Eigen::MatrixXd getSymmetricFactor(Rcpp::XPtr<HODLR_Tree> GPobj) {
  /*
   Explicitly builds and returns the symmetric factor W.
   GPobj -H matrix, already assembled and factorized
   Return W.
   */

  Eigen::MatrixXd returnW = GPobj->getSymmetricFactor();

  return returnW;
}

// [[Rcpp::export]]
Eigen::MatrixXd matmatProduct(Rcpp::XPtr<HODLR_Tree> GPobj, Eigen::MatrixXd M) {
  /*
   Obtain the matrix-matrix / matrix-vector product of the given matrix / vector M,
   with the matrix that is abstracted by the instance of Kernel.
   GPobj - H matrix, already assembled and factorized
   M - matrix or vector.
   Return GPobj M.
   */

  Eigen::MatrixXd returnMat = GPobj->matmatProduct(M);

  return returnMat;
}

// [[Rcpp::export]]
NumericVector solve_HODLR(Rcpp::XPtr<HODLR_Tree> GPobj, Mat b) {
  /*
   Applies the inverse of the matrix(abstracted by the Kernel object) on the given vector b.
   GPobj - H matrix, already assembled and factorized
   b - vector.
   Return GPobj^-1 b
   */

  Eigen::MatrixXd returnx = GPobj->solve(b);

  return wrap(returnx);
}
