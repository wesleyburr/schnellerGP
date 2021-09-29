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

#include "amatrix.h"
#include "h2compression.h"
#include "h2matrix.h"
#include "h2arith.h"
#include "matrixnorms.h"
#include "parameters.h"
#include "krylovsolvers.h"

class h2matrix{

     public:


     private:
          pclusterbasis rb, cb;
          arma::mat X; // co-ordinates of the GP
          pcovmatrix km;     // kernel matrix object
          ph2matrix  h2;      // h2matrix to use
          pcluster  root;
          pblock    broot;
          uint      *idx;
          uint      leafsize;

};
