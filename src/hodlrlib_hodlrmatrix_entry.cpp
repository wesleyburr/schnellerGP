#include <RcppEigen.h>
#include <Rcpp.h>
#include <squaredeMat.hpp> //squared exponential kernel
#include <maternMat.hpp> //squared exponential kernel
#include <HODLR_Matrix.hpp>
#include <HODLR_Tree.hpp>
#include <squaredMat_ptau.hpp>
#include <maternMat_ptau.hpp>
#include <maternMat_pX.hpp>
#include <iostream>
using namespace std; 

using namespace Rcpp;

//' get_n_levels
//'
//' Given integers N and M, computes n levels.
//'
//' @param N  Integer.
//' @param M  Integer.
//' @export
//' @return Returns log(N / M) / log(2).
//'
// [[Rcpp::export]]
int get_n_levels(int N, int M) {

  int n_levels = log(N / M) / log(2);
  return n_levels;

}

//******************************************************************************//

//' setup_compressedMatrixGP_Matern
//'
//' Set Up Compressed Gaussian Process with Matern Kernel
//'
//' Given a matrix of locations X, standard deviation sigma, and length
//' scale rho, sets up a compressed matrix for a Gaussian process using
//' a Matern kernel, with maximum submatrix size M. Default is the case $p=2$.
//'
//' @param X  Matrix of locations.
//' @param sigma  Function standard deviation.
//' @param rho  Length-Scale.
//' @param tol  Tolerance for accuracy of calculations.
//' @param M  Maximum sub-matrix size.
//' @export
//' @return After setting up compressed memory storage, returns a pointer
//'  to a HODLR_Tree object. Specifically, in this case, the Matern kernel
//'  (default) is the case of $p=2$, or $nu = 5/2$.
//'
//'    K(r) = sigma^2 * [1 + sqrt(5) * r / rho + 5/3 * (r / rho)^2] * exp[-sqrt(5) * r / rho]
//'
//'  For this implementation, the diagonal elements for the nugget are set to 1e-8.
// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_Matern(Mat X, double sigma, double rho,
                                                       double tol, int M) {
  bool is_sym = true;
  bool is_pd = true;
  int N = X.rows();
  Matern_Kernel* K = new Matern_Kernel(X, N, sigma, rho);
  int n_levels = get_n_levels(N, M);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T,true);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  delete K; 
  return ptr;

}

//' setup_compressedMatrixGP_MaternP1
//'
//' Set Up Compressed Gaussian Process with Matern (P1) Kernel
//'
//' Given a matrix of locations X, standard deviation sigma, and length
//' scale rho, sets up a compressed matrix for a Gaussian process using
//' a Matern kernel, with maximum submatrix size M. This is the $p=1$ case.
//'
//' @param X  Matrix of locations.
//' @param sigma  Function standard deviation.
//' @param rho  Length-Scale.
//' @param tol  Tolerance for accuracy of calculations.
//' @param M  Maximum sub-matrix size.
//' @export
//' @return After setting up compressed memory storage, returns a pointer
//'  to a HODLR_Tree object. Specifically, in this case, the Matern_p1_Kernel
//'  refers to a kernel with $p=1$ and $nu=3/2$:
//'
//'    K(r) = sigma^2 * [1 + sqrt(3) * r / rho] * exp[-sqrt(3) * r / rho]
//'
//'  For this implementation, the diagonal elements for the nugget are set to 1.
// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_MaternP1(Mat X,
                                                         double sigma,
                                                         double rho,
                                                         double tol,
                                                         int M) {
  bool is_sym = true;
  bool is_pd = true;
  int N = X.rows();
  Matern_p1_Kernel* K  = new Matern_p1_Kernel(X, N, sigma, rho);
  int n_levels = get_n_levels(N, M);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T,true);
  T -> assembleTree(is_sym, is_pd);
  T -> factorize();
  return ptr;

}

//' setup_compressedMatrixGP_Matern_tP
//'
//' Set Up Compressed Gaussian Process with Matern (tP) Kernel
//'
//' Given a matrix of locations X, standard deviation sigma, and length
//' scale rho, sets up a compressed matrix for a Gaussian process using
//' a Matern kernel, with maximum submatrix size M. This is $p=2$, with
//' user-provided reciprocal diagonal components for the nugget.
//'
//' @param X  Matrix of locations.
//' @param tP  Specified kernel.
//' @param sigma  Function standard deviation.
//' @param rho  Length-Scale.
//' @param tol  Tolerance for accuracy of calculations.
//' @param M  Maximum sub-matrix size.
//' @export
//' @return After setting up compressed memory storage, returns a pointer
//'  to a HODLR_Tree object. Specifically, in this case, the Matern_pX_Kernel
//'  refers to a kernel such that $p=2$, i.e.,
//'
//'    K(r) = sigma^2 * (1 + sqrt(5) * r / rho + 5/3 * (r / rho)^2) * exp(-sqrt(5) * r / rho)
//'
//'  but with the nugget diagonals being set to 1/P, for P a N*1 matrix.
// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_Matern_tP(Mat X,
                                                          Mat tP,
                                                          double sigma,
                                                          double rho,
                                                          double tol,
                                                          int M) {
  bool is_sym = true;
  bool is_pd = true;
  int N = X.rows();
  Matern_pX_Kernel* K = new Matern_pX_Kernel(X, N, sigma, rho, tP);
  int n_levels = get_n_levels(N, M);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T,true);
  delete K; 
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;

}

//' setup_compressedMatrixGP_sqrExpP1
//'
//' Set Up Compressed Gaussian Process with Square Exponential (P1) Kernel
//'
//' Given a matrix of locations X, standard deviation sigma, and length
//' scale rho, sets up a compressed matrix for a Gaussian process using
//' the Squared Exponential kernel.
//'
//' @param X  Matrix of locations.
//' @param tP   Specified kernel.
//' @param sigma  Function standard deviation.
//' @param rho  Length-Scale.
//' @param tol  Tolerance for accuracy of calculations.
//' @param M  Maximum sub-matrix size.
//' @export
//' @return After setting up compressed memory storage, returns a pointer
//'  to a HODLR_Tree object. Specifically, in this case, kernel is the
//'  Squared Exponential Kernel:
//'
//'    K(r) = sigma^2 * exp[ -rho * ||r||^2 ]
//'
//'  with the nugget diagonals set to 1.0, as with the Matern p1 case.
// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_sqrExp_tP(Mat X,
                                                          Mat tP,
                                                          double sigma,
                                                          double rho,
                                                          double tol,
                                                          int M) {
  bool is_sym = true;
  bool is_pd = true;
  int N = X.rows();
  SQRExponential_pX_Kernel* K  = new SQRExponential_pX_Kernel(X, N, sigma, rho, tP);
  int n_levels = get_n_levels(N, M);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T,true);
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  delete K; 
  return ptr;
  
}

//******************************************************************************//

//' setup_compressedMatrixGP_sqrExp
//'
//' Set Up Compressed Gaussian Process with Square Exponential Kernel
//'
//' Given a matrix of locations X, standard deviation sigma, and length
//' scale rho, sets up a compressed matrix for a Gaussian process using
//' the Squared Exponential kernel.
//'
//' @param X  Matrix of locations.
//' @param sigma  Function standard deviation.
//' @param rho  Length-Scale.
//' @param tol  Tolerance for accuracy of calculations.
//' @param M  Maximum sub-matrix size.
//' @export
//' @returns After setting up compressed memory storage, returns a pointer
//'  to a HODLR_Tree object. Specifically, in this case, kernel is the
//'  Squared Exponential Kernel:
//'
//'    K(r) = sigma^2 * exp[ -rho * ||r||^2 ]
//'
//'  with the nugget diagonals set to 1e-8, as with the Matern default.
// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_sqrExp(Mat X,
                                                       double sigma,
                                                       double rho,
                                                       double tol,
                                                       int M) {
  bool is_sym = true;
  bool is_pd = true;
  int N = X.rows();
  SQRExponential_Kernel* K = new SQRExponential_Kernel(X, N, sigma, rho);
  int n_levels = get_n_levels(N, M);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T,true);
  delete K; 
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;

}

//' setup_compressedMatrixGP_sqrExpP1
//'
//' Set Up Compressed Gaussian Process with Square Exponential (P1) Kernel
//'
//' Given a matrix of locations X, standard deviation sigma, and length
//' scale rho, sets up a compressed matrix for a Gaussian process using
//' the Squared Exponential kernel.
//'
//' @param X  Matrix of locations.
//' @param sigma  Function standard deviation.
//' @param rho  Length-Scale.
//' @param tol  Tolerance for accuracy of calculations.
//' @param M  Maximum sub-matrix size.
//' @export
//' @return After setting up compressed memory storage, returns a pointer
//'  to a HODLR_Tree object. Specifically, in this case, kernel is the
//'  Squared Exponential Kernel:
//'
//'    K(r) = sigma^2 * exp[ -rho * ||r||^2 ]
//'
//'  with the nugget diagonals set to 1.0, as with the Matern p1 case.
// [[Rcpp::export]]
Rcpp::XPtr<HODLR_Tree> setup_compressedMatrixGP_sqrExpP1(Mat X, double sigma, double rho,
                                                       double tol, int M) {
  bool is_sym = true;
  bool is_pd = true;
  int N = X.rows();
  SQRExponential_p1_Kernel* K  = new SQRExponential_p1_Kernel(X, N, sigma, rho);
  int n_levels = get_n_levels(N, M);
  HODLR_Tree* T = new HODLR_Tree(n_levels, tol, K);
  Rcpp::XPtr<HODLR_Tree> ptr(T,true);
  delete K; 
  T->assembleTree(is_sym, is_pd);
  T->factorize();
  return ptr;

}

//******************************************************************************//

//' simulate_compressedMatrixGP
//'
//' Simulate a draw from a compressed GP.
//'
//' @param N Number of observations to draw.
//' @param GPobj  The H matrix, assembled and factorized for use.
//' @export
//' @return  Given N and GPobj, generate N random N(0, 1)s, then
//'    use HODLR code to do a quick symmetric factor product. Returns
//'    the product, which are the draws from the process using compression.
//'
// [[Rcpp::export]]
NumericVector simulate_compressedMatrixGP(int N, Rcpp::XPtr<HODLR_Tree> GPobj) {

  NumericVector tW = rnorm(N, 0, 1);
  Eigen::Map<Eigen::MatrixXd> ttW(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(tW));
  Eigen::MatrixXd W = ttW;
  Eigen::MatrixXd returnV = GPobj -> symmetricFactorProduct(W);
  return wrap(returnV);

}

//' simulate_compressedMatrixGP_TP
//'
//' Simulate a draw from a compressed Gaussian Process with **something scaled**.
//'
//' @param N Number of observations to draw.
//' @param GPobj  The H matrix, assembled and factorized for use.
//' @param tP  A scale factor (or scale vector??) **
//' @export
//' @return  Given N and GPobj, generate N random N(0, 1)s, scale using tP, then
//'    use HODLR code to do a quick symmetric factor product. Returns
//'    the product, which are the draws from the process using compression.
//'
// [[Rcpp::export]]
NumericVector simulate_compressedMatrixGP_TP(int N,
                                             Rcpp::XPtr<HODLR_Tree> GPobj,
                                             NumericVector tP) {
  NumericVector tW = rnorm(N, 0, 1);
  tW = tW*sqrt(tP);
  Eigen::Map<Eigen::MatrixXd> ttW(Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(tW));
  Eigen::MatrixXd W = ttW;
  Eigen::MatrixXd returnV = GPobj -> symmetricFactorProduct(W);
  return wrap(returnV);
}

//' getSymmetricFactor
//'
//' Extract symmetric factor from an H matrix
//'
//' Build and return the symmetric factor W from an H matrix, already
//' assembled and factorized.
//'
//' @param GPobj  The H matrix, assembled and factorized for use.
//' @export
//' @return  Extracts the symmetric factor from the Gaussian Process
//' compressed object.
//'
// [[Rcpp::export]]
Eigen::MatrixXd getSymmetricFactor(Rcpp::XPtr<HODLR_Tree> GPobj) {

  Eigen::MatrixXd returnW = GPobj->getSymmetricFactor();
  return returnW;

}

//' matmatProduct
//'
//' Compute a matrix-matrix product
//'
//' Compute a matrix-matrix or matrix-vector product of the abstracted
//' kernel matrix and the provided matrix M. Overloads the matmatProduct()
//' function from the Eigen package.
//'
//' @param GPobj  The H matrix, assembled and factorized for use.
//' @param M  A matrix of doubles, as dynamically accessed via Eigen.
//' @export
//' @return  Computes the matrix-matrix product of GPobj and M, which works
//' even if M is a N*1 matrix.
//'
// [[Rcpp::export]]
Eigen::MatrixXd matmatProduct(Rcpp::XPtr<HODLR_Tree> GPobj, Eigen::MatrixXd M) {

  Eigen::MatrixXd returnMat = GPobj -> matmatProduct(M);
  return returnMat;

}

//' solve_HODLR
//'
//' Solve Ax=b using HODLR
//'
//' Solves the matrix equation Ax=b by using the kernel object's abstraction
//' inverse, i.e., A^{-1} * b.
//'
//' @param GPobj  The H matrix, assembled and factorized for use.
//' @param b  A matrix (vector), the right-hand side of the Ax=b equation.
//' @export
//' @return  Solves the matrix equation Ax=b for an A matrix which is already
//' assembled and factorized by HODLR.
//'
// [[Rcpp::export]]
NumericVector solve_HODLR(Rcpp::XPtr<HODLR_Tree> GPobj, Mat b) {

  Eigen::MatrixXd returnx = GPobj->solve(b);
  return wrap(returnx);
}


//' log_det_HODLR
//'
//' Returns the log-determinant of a  HODLR matrix
//'
//' Givne a HODLR matrix A returns \eqn{\log | A |}
//' 
//'
//' @param GPobj  The H matrix, assembled and factorized for use.
//' 
//' @export
//' @return \eqn{\log |A|} where \eqn{|A|} is the determinant of the matrix.,   
//' 
//'
// [[Rcpp::export]]
NumericVector log_det_HODLR(Rcpp::XPtr<HODLR_Tree> GPobj) {

  double  returnx = GPobj->logDeterminant(); 
  return wrap(returnx);
}

