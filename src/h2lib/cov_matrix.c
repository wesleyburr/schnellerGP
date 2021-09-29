/* ------------------------------------------------------------
 * This is the file "kernelmatrix.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2018
 * ------------------------------------------------------------ */

#include "cov_matrix.h"
#include "cov_matrix_build.h"
#include "h2matrix.h"
#include "basic.h"
#include "stdio.h"
#include "cov_matrix_build.h"

/* ------------------------------------------------------------
 * Create an empty @ref kernelmatrix object
 * ------------------------------------------------------------ */

pcovmatrix
new_covmatrix(uint dim, uint points, uint m, const double *hp)
{
  pcovmatrix km;
  real *x0;
  uint i;

  /* Initialize basic structure */
  km = (pcovmatrix) allocmem(sizeof(covmatrix));
  km->dim = dim;
  km->points = points;
  km->m = m;

  /* Empty kernel callback function */
  km->kernel = NULL;

  /* Initialize arrays for point coordinates */
  km->x = (real **) allocmem(sizeof(real *) * points);
  km->x[0] = x0 = allocreal(points * dim);
  for(i=1; i<points; i++) {
    x0 += dim;
    km->x[i] = x0;
  }
  //copy over hyper parameters
  for (i = 0; i < MAX_HPS; i++){
    km->hparm[i] = hp[i];
  }
  /* Initialize Chebyshev points */
  km->xi_ref = allocreal(m);
  for(i=0; i<m; i++)
    km->xi_ref[i] = REAL_COS(M_PI * (m - i - 0.5) / m);

  return km;
}

/* ------------------------------------------------------------
 * Delete a @ref kernelmatrix object
 * ------------------------------------------------------------ */

void
del_covmatrix(pcovmatrix km)
{
  if (km != NULL){
            freemem(km->x[0]);
            freemem(km->x);
            freemem(km);
  }
}

/* ------------------------------------------------------------
 * Create a clustergeometry object
 * ------------------------------------------------------------ */

pclustergeometry
creategeometry_covmatrix(pccovmatrix km)
{
  uint dim = km->dim;
  uint points = km->points;
  const real **x = (const real **) km->x;
  pclustergeometry cg;
  uint i, j;

  cg = new_clustergeometry(dim, points);

  for(i=0; i < points; i++)
    for(j=0; j <dim; j++)
      cg->x[i][j] = cg->smin[i][j] = cg->smax[i][j] = x[i][j];

  return cg;
}

/* ------------------------------------------------------------
 * Create transformed interpolation points for a bounding box
 * ------------------------------------------------------------ */

static real **
transform_points(uint dim, real *bmin, real *bmax,
                 uint m, const real *xi_ref)
{
  real **xi;
  real *xi0;
  real mid, rad;
  uint i, j;

  xi = (real **) allocmem(sizeof(real *) * dim);

  xi[0] = xi0 = (real *) allocmem(sizeof(real) * m * dim);
  for(i=1; i<dim; i++) {
    xi0 += m;
    xi[i] = xi0;
  }

  for(i=0; i<dim; i++) {
    mid = 0.5 * (bmax[i] + bmin[i]);
    rad = 0.5 * (bmax[i] - bmin[i]);

    for(j=0; j<m; j++){
      xi[i][j] = mid + rad * xi_ref[j];
    }
  }

  return xi;
}

/* ------------------------------------------------------------
 * Evaluate a Lagrange polynomial
 * ------------------------------------------------------------ */
static real
eval_lagrange(uint m, const real *xi, uint i, real t)
{
    real d, n;
    uint j;

    d = n = 1.0;

    for(j=0; j<i; j++) {
      d *= t - xi[j];
      n *= xi[i] - xi[j];
    }

    for(j=i+1; j<m; j++) {
      d *= t - xi[j];
      n *= xi[i] - xi[j];
    }

    return d / n;
}

/* ------------------------------------------------------------
 * Fill a leaf matrix
 * ------------------------------------------------------------ */
void print_cluster(pcluster c){

  if (c->sons){
    for (int i = 0; i < c->sons; i++){
      print_cluster(c->son[i]);
    }
  }else{
    printf("[%1.3f,%1.3f]\n",c->bmin[0],c->bmax[0]);
  }
}

void print_clusterB(pclusterbasis cb){

  if (cb->sons){
    for (int i = 0; i < cb->sons; i++){
      print_clusterB(cb->son[i]);
    }

  }else{
    print_matlab_amatrix(&(cb->V));
  }

}


/* ------------------------------------------------------------
 * Fill a dense matrix leaf node
 * This function gets each element in the kernel and directly fills it into
 * the tree
 * ------------------------------------------------------------ */
void
fillDense_covmatrix(const uint *ridx, const uint *cidx, pccovmatrix km,
                      pamatrix N)
  {
    const real **x = (const real **) km->x;
    uint rows = N->rows;
    uint cols = N->cols;
    pfield Na = N->a;
    longindex ldN = N->ld;
    uint i, j, ii, jj;
    /*TO DO: fix the bk to have all
     * of the COV kernel data right
     * now it just has the (x,y)
     * coordinates in the matrix.
     */
    uint bk[3];
    bk[2] = km->dim;
    if(ridx) {
      if(cidx) {
        for(j=0; j<cols; j++) {
          jj =    cidx[j];

          for(i=0; i<rows; i++) {
            ii = ridx[i];
            bk[0] = ii; bk[1] = jj;
            Na[i+j*ldN] = km->kernel(x[ii], x[jj], bk); /*km->data);*/
          }
        }
      }
      else {
        assert(cols <= km->points);

        for(j=0; j<cols; j++) {
          for(i=0; i<rows; i++) {
            ii = ridx[i];
            bk[0] = ii; bk[1] = j;
            Na[i+j*ldN] = km->kernel(x[ii], x[j], bk); /*km->data);*/
          }
        }
      }
    }else {
      assert(rows <= km->points);

      if(cidx) {
        for(j=0; j<cols; j++) {
          jj = cidx[j];

          for(i=0; i<rows; i++){
            bk[0] = i; bk[1] = jj;
          }
          Na[i+j*ldN] = km->kernel(x[i], x[jj], bk); /*km->data);*/
        }
      } else {
        assert(cols <= km->points);

        for(j=0; j<cols; j++){
          for(i=0; i<rows; i++){
            bk[0] = i; bk[1] = j;
          }
        }
        Na[i+j*ldN] = km->kernel(x[i], x[j], bk); /*km->data);*/
      }
    }
}



/*actually fill the damn H2 Matrix thing*/
void
fill_h2matrix_covmatrix(pccovmatrix km, ph2matrix G)
{
    uint rsons, csons;
    uint i, j;

    if(G->son) {
      rsons = G->rsons;
      csons = G->csons;

      for(j=0; j<csons; j++){
        for(i=0; i<rsons; i++){
          fill_h2matrix_covmatrix(km, G->son[i+j*rsons]);
        }
      }
    }
    else if(G->u){ /*admissable leaf compressed*/
      //fillS_kernelmatrix(G->rb->t, G->cb->t, km, &G->u->S);
      //printf("[%u]\n", G->rb->sons);
      //printf("[%u]\n", G->cb->sons);
      /*for (int i = 0; i < G->rb->t->size; i++){
        for (int j = 0; j < G->cb->t->size; j++){
          printf("[%u, %u]\n", G->rb->t->idx[i],G->cb->t->idx[j]);
        }
      }*/
      collectdense_h2matrix_covmatrix(km, G->rb,G->cb, &G->u->S);

      //collectdense_h2matrix(G,)
    } else { /*inadmissable leaf full rank*/
      assert(G->f);
      fillDense_covmatrix(G->rb->t->idx, G->cb->t->idx, km, G->f);
    }

}



