#include <RcppArmadillo.h>
#include <string>


#ifndef CLASS_MATRIX_H_LIB_H
#define CLASS_MATRIX_H_LIB_H

// files needed from H2lib

#include <iostream>

#include "hmatrix.h"
#include "harith.h"
#include "harith2.h"

#include "matrixnorms.h"
#include "parameters.h"
#include "krylovsolvers.h"
using namespace std; 

class matrix_h{
public:
     
     matrix_h(phmatrix ph);
     matrix_h(string file);
     ~matrix_h();

     void cholesky();
     double log_determinant(); 
     arma::vec multiply_hmat_vec(arma::vec v);
     arma::mat get_matrix();
     void save_h_netcdf(string file);
     
     arma::vec triang_invmul_vec(arma::vec xp); 
	
private:
     phmatrix hmatrix;
     bool is_cholesky;
};

#endif
