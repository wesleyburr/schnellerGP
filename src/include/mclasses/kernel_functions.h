
#include <iostream>
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

#ifndef _KERNEL_FUNCTIONS_H
#define _KERNEL_FUNCTIONS_H

enum kernel_type {gaussian = 1, matern = 2, diagonal = 3, gaussian_truncate =4};

HEADER_PREFIX  field
kernel_sqr_exp(const double *xx, const double *yy, void *data);


HEADER_PREFIX  field
kernel_diagonal(const double *xx, const double *yy, void *data);

HEADER_PREFIX field
kernel_sqr_exp_truncate(const double *xx, const double *yy, void *data);

HEADER_PREFIX field
kernel_matern_exp(const double *xx, const double *yy, void *data);

#endif