
#ifndef _COV_MATRIX_BUILD_H_
#define _COV_MATRIX_BUILD_H_

/* ------------------------------------------------------------
 * Active blocks,
 * i.e., blocks that have to be approximated by the current cluster
 * basis
 * ------------------------------------------------------------ */
#include "cov_matrix.h"
#include "settings.h"
#include "h2matrix.h"
#include "cluster.h"
#include "clustergeometry.h"

/* ------------------------------------------------------------
 * Active blocks,
 * i.e., blocks that have to be approximated by the current cluster
 * basis
 * ------------------------------------------------------------ */

typedef struct _compcovactive compcovactive;
typedef compcovactive *pcompcovactive;

struct _compcovactive {
  pccovmatrix G;			/* Source matrix */
  pcblock  b;			    /* Block */

  amatrix   A;			    /* Permuted or compressed submatrix */
  real      weight;		  /* Weight factor */
  pcompcovactive next;	/* Next block in list */
};

/* ------------------------------------------------------------
 * Passive blocks,
 * i.e., blocks that do not have to be approximated by the current
 * cluster basis, but may have to be approximated by its descendants
 * ------------------------------------------------------------ */

typedef struct _compcovpassive compcovpassive;
typedef compcovpassive *pcompcovpassive;

struct _compcovpassive {
  pccovmatrix G;
  pcblock  b;
  pcompcovpassive next;
};


HEADER_PREFIX pclusterbasis
  buildrowbasis_covmatrix(pccovmatrix G, pcblock b, pctruncmode tm, real eps);

HEADER_PREFIX pclusterbasis
  buildcolbasis_covmatrix(pccovmatrix G, pcblock b, pctruncmode tm, real eps);

HEADER_PREFIX void
  collectdense_h2matrix_covmatrix(pccovmatrix a, pcclusterbasis rb, pcclusterbasis cb,
                                  pamatrix s);

HEADER_PREFIX ph2matrix
build_projected_covmatrix_h2matrix(pccovmatrix G, pcblock b,
                                     pclusterbasis rb, pclusterbasis cb);

HEADER_PREFIX void
fill_kernelmatrix(const uint *ridx, const uint *cidx, pccovmatrix km,
                    pamatrix N);

#endif
