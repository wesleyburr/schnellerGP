#include <mclasses/matrix_h.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <getdiag.h>

using namespace std;

// Create new H-matrix from file reference, NetCDF
matrix_h::matrix_h(string file){
  hmatrix  = read_cdfcomplete_hmatrix(file.c_str());
}

// Create new H-matrix from phmatrix source; uses clone_hmatrix()
matrix_h::matrix_h(phmatrix ph){

  hmatrix = clone_hmatrix(ph);
  is_cholesky = false;
};

// Destructor call
matrix_h::~matrix_h(){
  del_hmatrix(hmatrix);
};

// Compute the triangular inverse multiplication by vector, for solve_
arma::vec matrix_h::triang_invmul_vec(arma::vec xp){
  pavector pv = new_avector(xp.n_rows);
  arma::vec rVal(xp.n_rows);

  for (uint i = 0; i < pv->dim; i++){
    pv->v[i] = xp[i];
  }

  if (is_cholesky){
    triangularsolve_hmatrix_avector	(	true,
                                      false,
                                      false,
                                      hmatrix,
                                      pv);
  } else {
    Rcpp::warning("Warning: Input H-matrix not Cholesky. Making a copy - this will add computation time.");
    // clone hmatrix, Cholesky decompose it, then
    // solve the equation using this copy
    phmatrix hmatrix2 = clone_hmatrix(hmatrix);
    ptruncmode tm = new_abseucl_truncmode();
    choldecomp_hmatrix(hmatrix2, tm, 1e-16);
    triangularsolve_hmatrix_avector	(	true,
                                      false,
                                      false,
                                      hmatrix2,
                                      pv);
    del_hmatrix(hmatrix2);
  }
  for (uint i = 0; i < pv->dim; i++){
    rVal[i]  = pv->v[i];
  }
  del_avector(pv);
  return rVal;
}


arma::vec matrix_h::multiply_hmat_vec(arma::vec v){
  pavector pv = new_avector(v.n_rows);
  pavector Y  = new_zero_avector(v.n_rows);

  for (uint i = 0; i < pv->dim; i++){
    pv->v[i] = v[i];
  }

  if (is_cholesky){
    // the hmatrix has been factored and only the lower triangular
    // part matters the upper part is crap.
    triangularmul_hmatrix_avector(true,
                                  false,
                                  false,
                                  hmatrix,
                                  pv);
    for (uint i = 0; i < pv->dim; i++){
      v[i] = pv->v[i];
    }
  }else{
    mvm_hmatrix_avector(1.0,false,hmatrix,pv,Y);
    for (uint i = 0; i < Y->dim; i++){
      v[i] = Y->v[i];
    }
  }

  del_avector(pv);
  del_avector(Y);

  return v;
};

arma::mat matrix_h::get_matrix(){
  uint rows = getrows_hmatrix(hmatrix);
  uint cols = getcols_hmatrix(hmatrix);
  pamatrix C = new_zero_amatrix(rows,cols);
  pamatrix B = new_identity_amatrix(rows,cols);
  addmul_hmatrix_amatrix_amatrix (1.0, false, hmatrix, false, B, false, C);

  arma::mat rv((int)rows,(int)cols);

  for (int i = 0; i < rv.n_rows; i++){
    for (int j = 0; j < rv.n_cols; j++){
      rv(i,j) = C->a[i+j*rv.n_cols];
    }
  }

  del_amatrix(C);
  del_amatrix(B);
  return rv;
}

/*
 *  @cholesky() - Cholesky decomposition of the current hmatrix
 *  in place.  That is the current hmatrix is overwritten.
 *
 *
 */
void matrix_h::cholesky(){
  ptruncmode tm = new_abseucl_truncmode();
  choldecomp_hmatrix(hmatrix,tm,1e-16);
  //choldecomp2_hmatrix(hmatrix,tm,1e-16);
  is_cholesky = true; //set this true so all subsequent operations
  // are done on the lower triangular part
  del_truncmode(tm);
}


double matrix_h::log_determinant(){
  //
  if (!is_cholesky){
    this->cholesky();
  }
  pavector diag = new_avector(getrows_hmatrix(hmatrix));
  getdiag_hmatrix(hmatrix, diag);
  double log_prod = 0.0;
  for (uint i = 0; i < diag->dim ; i ++){
    log_prod += log(diag->v[i]);
  }
  return log_prod;
}


void matrix_h::save_h_netcdf(string file){

  if (hmatrix){
      write_cdfcomplete_hmatrix(hmatrix,file.c_str());
  }
  return;
}
