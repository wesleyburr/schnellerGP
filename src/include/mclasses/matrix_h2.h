// GP.h
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifndef _MATRIX_H2_H
#define _MATRIX_H2_H
#include <iostream>
#include <string>
#include <omp.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

#include "settings.h"
#include "h2matrix.h"
#include "cluster.h"
#include "clustergeometry.h"
#include "truncation.h"
#include "settings.h"
#include "cov_matrix.h"
#include "cov_matrix_build.h"
#include "basic.h"

// files needed from H2lib
#include "hmatrix.h"
#include "harith.h"
#include "harith2.h"

#include "h2compression.h"
#include "h2matrix.h"
#include "h2arith.h"
#include "matrixnorms.h"
#include "parameters.h"
#include "krylovsolvers.h"

using namespace Rcpp;
using namespace std;
pclustergeometry
    creategeometry_amatrix(pcamatrix a);

struct diagonal_precond{
        uint n; // number of entries
        double *d;// diagonal entries
};

pclustergeometry  creategeometry_amatrix(pcamatrix a);
void add_diag_avector_h2lib(ph2matrix G, pavector a);

class matrix_h2{
public:

     matrix_h2(ph2matrix ph2);
     matrix_h2(string file);
     matrix_h2(const arma::mat imat,
                const double compression = 1e-10,
                const uint leafsize = 50,
                double eta = 1.0);//, const double compression=1e-10,
                                     // uint leafsize = 50, double eta = 1);

     matrix_h2(arma::mat tX, field  (*kernel_fun)(const double *xx, const double *yy, void *data),
                  arma::mat hp, uint tleafsize = 50, double eta = 2, double tolerance = 1e-10);

     ~matrix_h2();
    ////////////////////////////////////////////////////////////
     uint getdim(){
         return getrows_h2matrix(gp_h2);
     }
     void  draw_matrix(){
       cairo_t  *cr;
       printf("Draw hmatrix to \"hm_p1.pdf\"\n");
       cr = new_cairopng( 2048.0, 2048.0);
      // cr = new_cairopdf("./hm_p1.pdf", 1024.0, 1024.0);
       draw_cairo_h2matrix(cr, gp_h2, false, 0);
       write_cairopng(cr,"../hm_p1.png");
       cairo_destroy(cr);
     }
     arma::mat getMatrix(){
          pamatrix test_temp  = convert_h2matrix_amatrix (false,gp_h2);
          arma::mat rV(test_temp->a, test_temp->rows, test_temp->cols);
          del_amatrix(test_temp);
          return rV;
     }
     arma::umat getIndex(){


         arma::umat rVI(rows,1);
         for (unsigned int i = 0; i < rows; i++){
             rVI(i,0) = idx[i];
         }
//         cout <<"Always here" << endl;
         return rVI;
     }

     arma::vec multiply_vector(arma::vec v);
     void add_diag_avector( pavector a);
     arma::mat orderByIndex(arma::mat Y);

     ////////////////////////////////////////////////////////
     arma::vec solve_precond_vect(arma::vec tb,
                                 diagonal_precond *p,
                                 double eps=1e-14,
                                 uint iter=2000);


     arma::vec solve_precond_vect_wstart(arma::vec tb,
                                         diagonal_precond *p,
                                         arma::vec start,
                                         double eps=1e-14,
                                         uint iter=2000);

     phmatrix convert_to_hmatrix();

     void save_h2_netcdf(string file);
private:
     pclusterbasis rb, cb;
     arma::mat X;          // matrix points as well as cluster geometry object
     pclustergeometry cg;  // cluster geometry object
     pcovmatrix km;        // kernel matrix object
     ph2matrix gp_h2;      // h2matrix to use
     double eta;           // set approximation tolerence, defaults to 0.1;
     pcluster root;
     pblock broot;
     uint *idx;
     uint leafsize;
     uint rows,cols;
};

#endif

