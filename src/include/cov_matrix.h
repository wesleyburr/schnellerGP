
#include "cov_matrix_build.h"

#ifndef COVMATRIX_H
#define COVMATRIX_H
/** @defgroup kernelmatrix kernelmatrix
 *  @brief Approximation of kernel matrices.
 *  @{ */

/** @brief Data required to approximate a kernel matrix. */
typedef struct _covmatrix covmatrix;

/** @brief Pointer to a @ref kernelmatrix object. */
typedef covmatrix * pcovmatrix;

/** @brief Pointer to a constant @ref kernelmatrix object. */
typedef const covmatrix * pccovmatrix;

#include "settings.h"
#include "h2matrix.h"
#include "cluster.h"
#include "clustergeometry.h"

#define MAX_HPS 10

/** @brief Representation of a kernel matrix an its approximation.
 *
 *  A kernel matrix is a matrix with entries of the form
 *  @f$g_{ij} = k(x_i,x_j)@f$, where @f$k@f$ is a kernel functions
 *  and @f$(x_i)_{i\in\Idx}@f$ are points in a suitable space.
 *
 *  If the kernel functions is locally smooth, it can be approximated
 *  by interpolation, and this gives rise to blockwise low-rank
 *  approximations. */
struct _covmatrix {

        /** @brief Leading dimension, i.e., increment used to switch from one column to the next.  */
        uint ld;
        /** @brief Number of rows. */
        uint rows;
        /** @brief Number of columns.  */
        uint cols;
          /** @brief Spatial dimension. */
        uint dim;
          /** @brief Kernel function. */
        field (*kernel)(const real *xx, const real *yy, void *data);
        /** @brief Data for the kernel function. */
        double hparm[MAX_HPS];
          /** @brief Number of points. */
        uint points;
        /** @brief Coordinates of points. */
        real **x;
        /** @brief Interpolation order (i.e., number of interpolation points). */
        uint m;
        /** @brief Interpolation points for the reference interval @f$[-1,1@f$. */
        real *xi_ref;
};

struct gp_data{
        uint i; // row i dim of matrix
        uint j; // col j dim of matrix
        uint d; //dimensions of data
        double hparm[MAX_HPS]; // assuming no more than 10 hyper parameters
};

HEADER_PREFIX void
print_cluster(pcluster c);


HEADER_PREFIX void
print_clusterB(pclusterbasis cb);
/** @brief Create an empty @ref kernelmatrix object.
 *
 *  @param dim Spatial dimension.
 *  @param points Number of points.
 *  @param m Interpolation order.
 *  @returns New object. */
HEADER_PREFIX pcovmatrix
        new_covmatrix(uint dim, uint points, uint m, const double *hp);

/** @brief Delete a @ref covmatrix object.
 *
 *  @param km Object to be deleted. */
HEADER_PREFIX void
          del_covmatrix(pcovmatrix km);

/** @brief Create a @ref clustergeometry object for a @ref covmatrix
 *  object.
 *
 *  @param km Description of the cov matrix, particularly the points.
 *  @returns @ref clustergeometry object for the given points and dimension. */
HEADER_PREFIX pclustergeometry
          creategeometry_covmatrix(pccovmatrix km);

/** @brief Fill a @ref clusterbasis using interpolation.
 *
 *  @param km Description of the cov matrix.
 *  @param cb Cluster basis to be filled. */
HEADER_PREFIX void
          fill_clusterbasis_covmatrix(pccovmatrix km, pclusterbasis cb);

/** @brief Fill a @ref h2matrix using interpolation.
 *
 *  @param km Description of the cov matrix.
 *  @param G Matrix to be filled. */
HEADER_PREFIX void
          fill_h2matrix_covmatrix(pccovmatrix km, ph2matrix G);

/** @} */
#endif
