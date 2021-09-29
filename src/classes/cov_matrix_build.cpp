#include "stdio.h"
#include "cov_matrix.h"
#include "block.h"
#include "amatrix.h"
#include "cov_matrix_build.h"

/* ------------------------------------------------------------
 * Destructors
 * ------------------------------------------------------------ */

static void
del_compcovactive(pcompcovactive ca)
{
  pcompcovactive next;

  while (ca) {
    next = ca->next;

    uninit_amatrix(&ca->A);
    freemem(ca);

    ca = next;
  }
}

static void
del_compcovpassive(pcompcovpassive cp)
{
  pcompcovpassive next;

  while (cp) {
    next = cp->next;

    freemem(cp);

    cp = next;
  }
}
/* ***************************
 *
 *
 *
 * ****************************/
static void
addcol_comp(pccluster cc, pccovmatrix G, pcblock b, pctruncmode tm,
	    pcompcovactive * active, pcompcovpassive * passive)
{
  pcompcovactive ca;
  pcompcovpassive cp;
  pamatrix  Bhat;
  const uint *ridx, *cidx;
  uint      rsize, csize;
  real      norm, weight;
  uint      rsons, csons;
  uint      i, j;
  const real **x = (const real **) G->x;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    /* Check whether there is a son matching the cluster cc */
    j = 0;
    while (j < csons && b->son[j * rsons]->cc != cc){
      j++;
    }
    /* If there is, check the sons recursively */
    if (j < csons) {
      for (i = 0; i < rsons; i++){
          addcol_comp(cc, G, b->son[i + j * rsons], tm, active, passive);
      }
    } else {			/* Otherwise, this matrix block is passive */
      assert(b->cc == cc);

      cp = (pcompcovpassive) allocmem(sizeof(compcovpassive));
      cp->G = G;
      cp->b = b;
      cp->next = *passive;
      *passive = cp;
    }
  }
  else if (b->a) {
    assert(b->cc == cc);

    /* Create an active block */
    rsize = b->rc->size;
    csize = b->cc->size;
    ca = (pcompcovactive) allocmem(sizeof(compcovactive));
    ca->G = G;
    ca->b = b;
    Bhat = init_amatrix(&ca->A, csize, rsize);
    ca->next = *active;
    *active = ca;

    /* Copy entries from original matrix G */
    ridx = b->rc->idx;
    cidx = b->cc->idx;

    gp_data gpD;
    gpD.d = G->dim;
    for (int z = 0; z < MAX_HPS; z++){ // there are no more than MAX_HP
         gpD.hparm[z] = G->hparm[z];   // for the kernel
    }

    for (j = 0; j < rsize; j++){
      for (i = 0; i < csize; i++){
           gpD.i = ridx[j]; gpD.j = cidx[i];
           Bhat->a[i + j * Bhat->ld] = G->kernel(x[gpD.i],x[gpD.j],&gpD);
       /*  uint exD[3]; exD[0] =  ridx[j];
          exD[1] = cidx[i]; exD[2] = G->dim;
          Bhat->a[i + j * Bhat->ld] = G->kernel(x[exD[0]],x[exD[1]],&exD);*/
          //Bhat->a[i + j * Bhat->ld] = G->a[ridx[j] + cidx[i] * G->ld];
      }
    }
    /* Compute weight factor if necessary */
    weight = 1.0;
    if (tm && tm->blocks) {
      norm = (tm->frobenius ? normfrob_amatrix(Bhat) : norm2_amatrix(Bhat));
      if (norm > 0.0){
          weight = 1.0 / norm;
      }
    }
    ca->weight = weight;
  }
}
/*
 *
 *
 *
 *
 *
 */
 /* ------------------------------------------------------------
 * Add a row block to the lists for a given cluster.
 * ------------------------------------------------------------ */
static void
addrow_comp(pccluster rc, pccovmatrix G, pcblock b, pctruncmode tm,
	    pcompcovactive * active, pcompcovpassive * passive)
{
  pcompcovactive ca;
  pcompcovpassive cp;
  pamatrix  Ahat;
  const uint *ridx, *cidx;
  size_t    ldA, ldG;
  uint      rsize, csize;
  real      norm, weight;
  uint      rsons, csons;
  uint      i, j;
  const real **x = (const real **) G->x;

  if (b->son) {
    rsons = b->rsons;
    csons = b->csons;

    /* Check whether there is a son matching the cluster rc */
    i = 0;
    while (i < rsons && b->son[i]->rc != rc){
      i++;
    }
    /* If there is, check the sons recursively */
    if (i < rsons) {
      for (j = 0; j < csons; j++){
     	addrow_comp(rc, G, b->son[i + j * rsons], tm, active, passive);
      }
    }
    else {			/* Otherwise, this matrix block is passive */
      assert(b->rc == rc);

      cp = (pcompcovpassive) allocmem(sizeof(compcovpassive));
      cp->G = G;
      cp->b = b;
      cp->next = *passive;
      *passive = cp;
    }
  }else if (b->a) {
    assert(b->rc == rc);

    /* Create an active block */
    rsize = b->rc->size;
    csize = b->cc->size;
    ca = (pcompcovactive) allocmem(sizeof(compcovactive));
    ca->G = G;
    ca->b = b;

    Ahat = init_amatrix(&ca->A, rsize, csize);
    ca->next = *active;
    *active = ca;

    /* Copy entries from original matrix G */
    ldG = G->ld;
    ldA = Ahat->ld;
    ridx = b->rc->idx;
    cidx = b->cc->idx;

    gp_data gpD;
    gpD.d = G->dim;
    for (int z = 0; z < MAX_HPS; z++){ // there are no more than MAX_HP
         gpD.hparm[z] = G->hparm[z];   // for the kernel
    }

    for (j = 0; j < csize; j++){
      for (i = 0; i < rsize; i++){
           // FIX ME - include other information lenghtscale etc.
          gpD.i = ridx[i]; gpD.j = cidx[j];
          Ahat->a[i + j * ldA] = G->kernel(x[gpD.i],x[gpD.j],&gpD);
      }
     }
    /* Compute weight factor if necessary */
    weight = 1.0;
    if (tm && tm->blocks) {
      norm = (tm->frobenius ? normfrob_amatrix(Ahat) : norm2_amatrix(Ahat));
      if (norm > 0.0){
          weight = 1.0 / norm;
      }
    }
    ca->weight = weight;
  }
}



/* ***************************\
 *
 *
 *
 *
 *
 *
 * **************************/
 static    pclusterbasis
buildbasis_comp(pccluster t, bool colbasis,
		     pcompcovactive active, pcompcovpassive passive, pctruncmode tm,
		     real eps)
{
       pclusterbasis cb, cb1;
       amatrix   tmp1, tmp2, tmp3, tmp4;
       realavector tmp5;
       pamatrix  Ahat, Ahat0, Ahat1;
       pamatrix  Q, Q1;
       prealavector sigma;
       pcompcovactive active1, ca, ca1;
       pcompcovpassive passive1, cp;
       real      zeta_age, zeta_level;
       uint      i, off, m, n, k;

       zeta_age = (tm ? tm->zeta_age : 1.0);
       zeta_level = (tm ? tm->zeta_level : 1.0);

       cb = new_clusterbasis(t);

       if (cb->sons > 0) {
         assert(cb->sons == t->sons);

         off = 0;
         m = 0;
         for (i = 0; i < t->sons; i++) {
           active1 = 0;
           passive1 = 0;

           /* Check for passive blocks that become active in son */
           if (colbasis) {
               for (cp = passive; cp; cp = cp->next){
                 addcol_comp(t->son[i], cp->G, cp->b, tm, &active1, &passive1);
               }
            }else {
               for (cp = passive; cp; cp = cp->next){
                 addrow_comp(t->son[i], cp->G, cp->b, tm, &active1, &passive1);
               }
            }

                /* Add submatrices for already active blocks to list */
            for (ca = active; ca; ca = ca->next) {
                    ca1 = (pcompcovactive) allocmem(sizeof(compcovactive));
                    ca1->G = ca->G;
                    ca1->b = ca->b;
                    init_sub_amatrix(&ca1->A, &ca->A, t->son[i]->size, off, ca->A.cols,
                               0);
                    ca1->weight = ca->weight * zeta_age;
                    ca1->next = active1;
                    active1 = ca1;
            }

           /* Create cluster basis for son */
           cb1 = buildbasis_comp(t->son[i], colbasis, active1, passive1, tm,
                                eps * zeta_level);
           ref_clusterbasis(cb->son + i, cb1);

           /* Clean up block lists */
           del_compcovactive(active1);
           del_compcovpassive(passive1);

           off += t->son[i]->size;
           m += cb1->k;
         }

         assert(off == t->size);

         for (ca = active; ca; ca = ca->next) {
           Ahat = init_amatrix(&tmp1, m, ca->A.cols);

           /* Combine projections of son's submatrices */
           off = 0;
           m = 0;
           for (i = 0; i < t->sons; i++) {
               /* Matrix prepared by recursive call */
               Ahat1 = init_sub_amatrix(&tmp2, &ca->A, cb->son[i]->k, off, Ahat->cols, 0);

               /* Appropriate submatrix in Ahat */
               Ahat0 = init_sub_amatrix(&tmp3, Ahat, cb->son[i]->k, m, Ahat->cols, 0);

               /* Copy submatrix to its place */
               copy_amatrix(false, Ahat1, Ahat0);
               uninit_amatrix(Ahat1);
               uninit_amatrix(Ahat0);

               off += t->son[i]->size;
               m += cb->son[i]->k;
           }
           assert(off == t->size);
           assert(m == Ahat->rows);

           copy_sub_amatrix(false, Ahat, &ca->A);
           uninit_amatrix(Ahat);
         }
       } else{
         m = t->size;
       }
       /* Determine number of columns of SVD matrix */
       n = 0;
       for (ca = active; ca; ca = ca->next){
         n += ca->A.cols;
       }

       /* Quick exit if zero columns */
       if (n == 0) {
         resize_clusterbasis(cb, 0);
         return cb;
       }

       /* Combine submatrices */
       Ahat = init_amatrix(&tmp1, m, n);
       n = 0;
       for (ca = active; ca; ca = ca->next) {
         Ahat0 = init_sub_amatrix(&tmp2, Ahat, m, 0, ca->A.cols, n);

         copy_sub_amatrix(false, &ca->A, Ahat0);
         if (ca->weight != 1.0){
           scale_amatrix(ca->weight, Ahat0);
         }
         uninit_amatrix(Ahat0);

         n += ca->A.cols;
       }
       assert(n == Ahat->cols);

       /* Compute singular value decomposition */
       k = UINT_MIN(m, n);
       Q = init_amatrix(&tmp2, m, n);
       sigma = init_realavector(&tmp5, k);
       svd_amatrix(Ahat, sigma, Q, 0);
       uninit_amatrix(Ahat);

       /* Find appropriate rank */
       k = findrank_truncmode(tm, eps, sigma);
       uninit_realavector(sigma);

       /* Set rank of new cluster basis */
       resize_clusterbasis(cb, k);

       if (cb->sons == 0) {
         assert(cb->V.rows == Q->rows);
         assert(cb->V.cols == k);

         /* Copy cluster basis matrix to new cluster basis */
         Q1 = init_sub_amatrix(&tmp1, Q, Q->rows, 0, k, 0);
         copy_amatrix(false, Q1, &cb->V);
         uninit_amatrix(Q1);
       }
       else {
         off = 0;
         for (i = 0; i < cb->sons; i++) {
           assert(cb->son[i]->E.rows == cb->son[i]->k);
           assert(cb->son[i]->E.cols == k);

           /* Copy transfer matrix to new cluster basis */
           Q1 = init_sub_amatrix(&tmp1, Q, cb->son[i]->k, off, k, 0);
           copy_amatrix(false, Q1, &cb->son[i]->E);
           uninit_amatrix(Q1);

           off += cb->son[i]->k;
         }
         assert(off == m);
       }

       /* Prepare projection matrices for father */
       Q1 = init_sub_amatrix(&tmp3, Q, m, 0, k, 0);
       for (ca = active; ca; ca = ca->next) {
         /* Copy original matrix */
         Ahat = init_amatrix(&tmp1, m, ca->A.cols);
         copy_sub_amatrix(false, &ca->A, Ahat);

         /* Pick submatrix for result */
         Ahat0 = init_sub_amatrix(&tmp4, &ca->A, k, 0, ca->A.cols, 0);
         clear_amatrix(Ahat0);

         /* Compute projection */
         addmul_amatrix(1.0, true, Q1, false, Ahat, Ahat0);
         uninit_amatrix(Ahat0);
         uninit_amatrix(Ahat);
       }
       uninit_amatrix(Q1);
       uninit_amatrix(Q);

       /* Now we're done with this subtree */
       update_clusterbasis(cb);

       return cb;
}

pclusterbasis
buildrowbasis_covmatrix(pccovmatrix G, pcblock b, pctruncmode tm, real eps)
{
  pclusterbasis rb;
  pcompcovactive active;
  pcompcovpassive passive;

  active = 0;
  passive = 0;
  addrow_comp(b->rc, G, b, tm, &active, &passive);

  rb = buildbasis_comp(b->rc, false, active, passive, tm, eps);

  del_compcovactive(active);
  del_compcovpassive(passive);

  return rb;
}

pclusterbasis
buildcolbasis_covmatrix(pccovmatrix G, pcblock b, pctruncmode tm, real eps)
{
  pclusterbasis cb;
  pcompcovactive active;
  pcompcovpassive passive;

  active = 0;
  passive = 0;
  addcol_comp(b->cc, G, b, tm, &active, &passive);

  cb = buildbasis_comp(b->cc, true, active, passive, tm, eps);

  del_compcovactive(active);
  del_compcovpassive(passive);

  return cb;
}

/* ------------------------------------------------------------
 * Orthogonal projection
 * ------------------------------------------------------------ */
void
collectdense_h2matrix_covmatrix(pccovmatrix a, pcclusterbasis rb, pcclusterbasis cb,
                                pamatrix s)
{
    const real **x = (const real **) a->x;
    amatrix   tmp1, tmp2;
    pamatrix  s1, s2;
    pccluster row, col;
    uint      i, j;

    assert(s->rows == rb->k);
    assert(s->cols == cb->k);

    clear_amatrix(s);

{

    if (rb->sons > 0) {
        if (cb->sons > 0) {		/* rb has sons, cb has sons */

            for (j = 0; j < cb->sons; j++) {
              s1 = init_amatrix(&tmp1, rb->k, cb->son[j]->k);
              clear_amatrix(s1);

                for (i = 0; i < rb->sons; i++) {
                  s2 = init_amatrix(&tmp2, rb->son[i]->k, cb->son[j]->k);
                  collectdense_h2matrix_covmatrix(a, rb->son[i], cb->son[j], s2);
                  addmul_amatrix(1.0, true, &rb->son[i]->E, false, s2, s1);
                  uninit_amatrix(s2);
                }
              addmul_amatrix(1.0, false, s1, false, &cb->son[j]->E, s);
              uninit_amatrix(s1);
            }
        } else {			/* rb has sons, cb is a leaf */
            for (i = 0; i < rb->sons; i++) {
              s2 = init_amatrix(&tmp2, rb->son[i]->k, cb->k);
              collectdense_h2matrix_covmatrix(a, rb->son[i], cb, s2);
              addmul_amatrix(1.0, true, &rb->son[i]->E, false, s2, s);
              uninit_amatrix(s2);
            }
        }
    } else {

      if (cb->sons > 0) {		/* rb is a leaf, cb has sons */
       //   #pragma omp for
          for (j = 0; j < cb->sons; j++) {
            s1 = init_amatrix(&tmp1, rb->k, cb->son[j]->k);
            collectdense_h2matrix_covmatrix(a, rb, cb->son[j], s1);
            addmul_amatrix(1.0, false, s1, false, &cb->son[j]->E, s);
            uninit_amatrix(s1);
          }
      }else {			/* rb is a leaf, cb is a leaf */
          row = rb->t;
          col = cb->t;
          s1 = init_amatrix(&tmp1, rb->k, col->size);
          clear_amatrix(s1);
          s2 = init_amatrix(&tmp2, row->size, col->size);
          gp_data gpD;
          gpD.d = a->dim;
          for (int z = 0; z < MAX_HPS; z++){ // there are no more than MAX_HP
               gpD.hparm[z] = a->hparm[z];   // for the kernel
          }


          for (j = 0; j < col->size; j++){
            for (i = 0; i < row->size; i++){
              gpD.i = row->idx[i];
              gpD.j = col->idx[j];

               //FIX ME
              s2->a[i + j * s2->ld] = a->kernel(x[gpD.i],x[gpD.j],&gpD);
              // a[row->idx[i] + col->idx[j] * a->ld];
            }
          }
          addmul_amatrix(1.0, true, &rb->V, false, s2, s1);
          addmul_amatrix(1.0, false, s1, false, &cb->V, s);

          uninit_amatrix(s2);
          uninit_amatrix(s1);
      }
    }
}
}

/* ------------------------------------------------------------
 * Orthogonal projection
 * ------------------------------------------------------------ */

ph2matrix
build_projected_covmatrix_h2matrix(pccovmatrix G, pcblock b,
                                   pclusterbasis rb, pclusterbasis cb)
{
          const real **x = (const real **) G->x;
          ph2matrix h2, h21;
          pamatrix  f;
          pclusterbasis rb1, cb1;
          pccluster rc, cc;
          const uint *ridx, *cidx;
          uint      rsize, csize;
          uint      rsons, csons;
          uint      i, j;

          if (b->son) {
                    rsons = b->rsons;
                    csons = b->csons;
                    h2 = new_super_h2matrix(rb, cb, rsons, csons);

                    for (j = 0; j < csons; j++) {
                              cb1 = cb;
                              if (b->son[0]->cc != b->cc) {
                                        assert(j < cb->sons);
                                        cb1 = cb->son[j];
                              }

                              for (i = 0; i < rsons; i++) {
                                        rb1 = rb;
                                        if (b->son[0]->rc != b->rc) {
                                                  assert(i < rb->sons);
                                                  rb1 = rb->son[i];
                                        }

                                        h21 = build_projected_covmatrix_h2matrix(G, b->son[i + j * rsons], rb1,
                                                                                 cb1);
                                        ref_h2matrix(h2->son + i + j * rsons, h21);
                              }
                    }
          }else if (b->a) {
                    h2 = new_uniform_h2matrix(rb, cb);
                    collectdense_h2matrix_covmatrix(G, rb, cb, &h2->u->S);
          }else {
                    h2 = new_full_h2matrix(rb, cb);
                    f = h2->f;

                    rc = rb->t;
                    cc = cb->t;
                    ridx = rc->idx;
                    cidx = cc->idx;
                    rsize = rc->size;
                    csize = cc->size;

                    assert(h2->f->rows == rc->size);
                    assert(h2->f->cols == cc->size);

                    gp_data gpD;
                    gpD.d = G->dim;
                    for (int z = 0; z < MAX_HPS; z++){ // there are no more than MAX_HP
                         gpD.hparm[z] = G->hparm[z];   // for the kernel
                    }

                    for (j = 0; j < csize; j++){
                         for (i = 0; i < rsize; i++){
                              gpD.i = ridx[i]; gpD.j = cidx[j];
                              //FIX ME

                              f->a[i + j * f->ld] = G->kernel(x[gpD.i],x[gpD.j],&gpD);
                              //G->a[ridx[i] + cidx[j] * G->ld];
                         }
                    }
          }

          update_h2matrix(h2);
          return h2;
}


