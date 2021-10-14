// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "mclasses/matrix_h2.h"
#include "settings.h"
#include <iostream>
#include <string>
#include "stdio.h"
#include "parameters.h"		/* Read parameters interactively */


pclustergeometry
creategeometry_amatrix(pcamatrix a){

    uint points = a->rows;
    pclustergeometry cg;
    uint i;

    cg = new_clustergeometry(1, points);

    for(i=0; i < points; i++){
      cg->x[i][0] = cg->smin[i][0] = cg->smax[i][0] = double(i)/double(points);
    }

    return cg;
}

static void
diag_pre(void *pdata, pavector r)
{
    diagonal_precond  *A = (diagonal_precond*) pdata;
    for (uint i = 0; i < A->n; i++){
      r->v[i] = r->v[i]/A->d[i];
    }

}


matrix_h2::matrix_h2(string file){
  gp_h2  = read_cdfcomplete_h2matrix(file.c_str());
  rb = clone_clusterbasis(gp_h2->rb);
  cb = clone_clusterbasis(gp_h2->cb);
  cg = NULL;
  idx = NULL;
  root  = NULL;
  broot = NULL;
}


matrix_h2::matrix_h2(ph2matrix ph2){
    rb = clone_clusterbasis(ph2->rb);
    cb = clone_clusterbasis(ph2->cb);
    gp_h2  = clone_h2matrix(ph2, rb, cb);
    cg = NULL;
    idx = NULL;
    root  = NULL;
    broot = NULL;
}


matrix_h2::matrix_h2(const arma::mat imat,
                     const double compression,
                     const uint leafsize,
                     double eta){

  pamatrix matrix  =  new_amatrix(imat.n_rows,imat.n_cols);


  for (int i = 0; i < imat.n_rows; i++){
    for(int j =0; j < imat.n_cols; j++){
      setentry_amatrix(matrix,i,j, imat(i,j));
    }
  }

  rows = imat.n_rows; cols = imat.n_cols;
  idx = (uint *) allocmem(sizeof(uint) * imat.n_rows);

  for(int i=0; i < imat.n_rows; i++){
    idx[i] = i;
  }

  cg = creategeometry_amatrix(matrix);

  root  = build_adaptive_cluster (cg, imat.n_rows , idx, leafsize);
  broot = build_strict_block(root, root, &eta, admissible_max_cluster);

  cb = NULL; // don't need to mess with this
  rb = NULL;
  ptruncmode trunc_mode = new_abseucl_truncmode();
  gp_h2 = compress_amatrix_h2matrix(matrix,broot, trunc_mode, compression);
  del_truncmode(trunc_mode);
}

// constructor to build the x variable
matrix_h2::matrix_h2(arma::mat tX,
                     field (*kernel_fun)(const double *xx, const double *yy, void *data),
                     arma::mat hp,
                     uint tleafsize,
                     double eta,
                     double tolerance){
    pcovmatrix km;     // kernel matrix object
    leafsize = tleafsize;
    X = tX;

    double *hyp_p = new double[MAX_HPS];

    for (int i = 0; i < MAX_HPS; i++){

      if (i >= hp.n_rows){
        hyp_p[i] = 0;
      }else{
        hyp_p[i] = hp(i,0);
      }
    }

    km = new_covmatrix(X.n_cols , X.n_rows, 3, hyp_p); // size of the matrix default interpolation order = 4
    km->kernel = kernel_fun;

    delete(hyp_p);

    for (int i = 0; i < X.n_rows; i++){
      for (int j = 0; j < X.n_cols; j++){
        km->x[i][j] = X(i,j);
      }
    }

    cg = creategeometry_covmatrix(km);

    idx = (uint *) allocmem(sizeof(uint) * X.n_rows);
    rows =  cols = X.n_rows;
    for(int i=0; i < X.n_rows; i++){
      idx[i] = i;
    }

    root = build_pca_cluster(cg, X.n_rows, idx, leafsize);
    broot = build_strict_block(root, root, &eta, admissible_max_cluster );
    ptruncmode trunc_mode = 	new_abseucl_truncmode();
    rb = buildrowbasis_covmatrix (km, broot, trunc_mode, tolerance);
    cb = buildcolbasis_covmatrix (km, broot, trunc_mode, tolerance);
    gp_h2 = build_projected_covmatrix_h2matrix(km,broot, rb, cb);
    del_truncmode(trunc_mode);
    del_covmatrix(km);

}

matrix_h2::~matrix_h2(){
    // remove the allocated h2lib objects
    // which include the h2lib tree and all corresponding
    // objects.
    if (idx !=NULL){
      delete idx;
      idx = NULL;
    }

    if (cg != NULL){
      del_clustergeometry(cg);
      cg = NULL;
    }
    if (root != NULL){
      del_cluster(root);
      root = NULL;
    }
    if (broot !=NULL){
      del_block(broot);
      broot = NULL;
    }
    if (gp_h2 != NULL){
      del_h2matrix(gp_h2);
      gp_h2 = NULL;
    }
}

////////////////////////////
//
//  Defaults to iter = 2000 and eps = 1e-14 if not otherwise set, via
//  the defaults in solve_pgmres_h2matrix_avector() and in the .h file.
//  Overwritten by user-passed eps and iter in this.
//
//  Call from solve_h2_vector():
//      arma::vec rv = gp -> solve_precond_vect(b, &pc, tolerance, iter);
arma::vec matrix_h2::solve_precond_vect(arma::vec tb,
                                        diagonal_precond *p,
                                        double eps,
                                        uint iter){

  // tb is the b vector for the equation Ax=b; A is the parent 'owner'
  // of this method (as it is hooked off matrix_h2::)
  pavector b = new_pointer_avector((double*)tb.memptr(), tb.n_rows);
  // initialize x as a vector of 0s
  pavector x = new_zero_avector(tb.n_rows);


  // Parameters
  // A	      System matrix, should be invertible.
  // prcd	    Callback function for preconditioner $N$.
  // pdata  	Data for prcd callback function.
  // b	      Right-hand side vector.
  // x	      Initial guess, will be overwritten by approximate solution.
  // eps	    Relative accuracy $\epsilon$, the method stops if $\|N(Ax-b)\|_2 \leq \epsilon \|N b\|_2$.
  // maxiter	Maximal number of iterations. maxiter=0 means that the number of iterations is not bounded.
  // kmax	    Maximal dimension of Krylov subspace.
  uint ik = solve_pgmres_h2matrix_avector(gp_h2,
                                          diag_pre,
                                          p,
                                          b,  // RHS
                                          x,  // initial guess of 0
                                          eps,
                                          iter,
                                          200);

  if (ik == iter){ warning("Solver failed to converge!."); };
  pstopwatch sw;	/* Stopwatch for time measuring */
  double time;		/* Variable for time measurement */

  arma::mat retVal(x -> v, x -> dim, 1);
  del_avector(b);
  del_avector(x);
  return retVal;
}

////////////////////////////
//
//  Defaults to iter = 2000 and eps = 1e-14 if not otherwise set, via
//  the defaults in solve_pgmres_h2matrix_avector(). Overwritten by user-passed
//  eps and iter in this.
//
//  Call from solve_h2_vector_start():
//      arma::vec rv = gp -> solve_precond_vect_wstart(b, &pc, tolerance, iter, start);
arma::vec matrix_h2::solve_precond_vect_wstart(arma::vec tb,
                                               diagonal_precond *p,
                                               arma::vec start,
                                               double eps,
                                               uint iter){

  pavector b = new_pointer_avector((double*)tb.memptr(),tb.n_rows);
  // starting point for x, passed below
  pavector x = new_pointer_avector((double*)start.memptr(), start.n_rows);

  uint ik = solve_pgmres_h2matrix_avector(gp_h2,
                                          diag_pre,
                                          p,
                                          b,  // RHS
                                          x,  // starting point via start
                                          eps,
                                          iter,
                                          200);
  if (ik == iter){warning("Solver failed to converge!.");};
  pstopwatch sw;	/* Stopwatch for time measuring */
  double time;		/* Variable for time measurement */

  arma::vec retVal(x -> v, x -> dim, 1);

  del_avector(b);
  del_avector(x);

  return retVal;
}

// compute gp_h2 * v, for provided parent H2 matrix gp_h2 - in
// the internal notation, gp_h2 * pv = y, multiply to obtain y
arma::vec matrix_h2::multiply_vector(arma::vec v){
  arma::vec rv(v.n_rows);

  pavector pv = new_pointer_avector((double*)v.memptr(), v.n_rows);
  pavector y =  new_zero_avector(v.n_rows);

  // addeval_h2matrix_avector(field alpha, pch2matrix h2, pcavector x, pavector y)
  addeval_h2matrix_avector( 1.0, gp_h2, pv, y);

  for (int i = 0; i < v.n_rows; i++ ){
    rv[i] = y->v[i];
  }

  del_avector(pv); del_avector(y);
  return rv;
}


void add_diag_avector_h2lib(ph2matrix G, pavector a){
  const uint *ridx, *cidx;
  uint rsize, csize;
  uint rsons, csons;

  if (G->son) {
    for (uint j = 0; j < G->csons; j++) {
      for (uint i = 0; i < G->rsons; i++) {
        add_diag_avector_h2lib(G->son[i + j * G->rsons], a);
      }
    }
  }else if (G->u){

  }else if (G->f) {
    //only worry about full rank leaves
    ridx = G->rb->t->idx;
    cidx = G->cb->t->idx;
    rsize = G->rb->t->size;
    csize = G->cb->t->size;

    //print_matlab_amatrix(G->f);
    for (uint j = 0; j < csize; j++){
      for (uint i = 0; i < rsize; i++){
        // we are at the diagonal
        if (ridx[i]==cidx[j]){
          G->f->a[i + j * G->f->ld] = G->f->a[i + j * G->f->ld] + a->v[ridx[i]];

        }
      }
    }
    //update_h2matrix(G);
  }
}

void matrix_h2::add_diag_avector( pavector a){
    add_diag_avector_h2lib(gp_h2, a);
}

phmatrix matrix_h2::convert_to_hmatrix(){
  phmatrix pm = convert_h2matrix_hmatrix (gp_h2);
  return pm;
};

void matrix_h2::save_h2_netcdf(string file){
  char save_file[file.size() + 1];
  file.copy(save_file,file.size()+1);
  if ( gp_h2){
     write_cdfcomplete_h2matrix(gp_h2,file.c_str());
  }
  return;
}
