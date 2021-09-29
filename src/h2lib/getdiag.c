
#include "avector.h"
#include "hmatrix.h"

void
getdiag_hmatrix(pchmatrix hm, pavector diag)
{
  pfield fa;
  longindex ldf;
  uint cols;
  uint i, ii;

  assert(hm->rc == hm->cc);
  
  if(hm->f) {
    fa = hm->f->a;
    ldf = hm->f->ld;
    cols = hm->f->cols;
    
    for(i=0; i<cols; i++) {
      ii = hm->rc->idx[i];
      
      assert(ii < diag->dim);
      
      diag->v[ii] = fa[i+i*ldf];
    }
  }
  else {
    assert(hm->son);
    assert(hm->rsons == hm->csons);

    for(i=0; i<hm->csons; i++)
      getdiag_hmatrix(hm->son[i+i*hm->rsons], diag);
  }
}

static void
getdiag(pchmatrix hm, pavector diagp, uint off)
{
  pfield fa;
  longindex ldf;
  uint cols;
  uint i, off1;

  assert(hm->rc == hm->cc);

  if(hm->f) {
    assert(hm->f->cols + off <= diagp->dim);
    
    fa = hm->f->a;
    ldf = hm->f->ld;
    cols = hm->f->cols;

    for(i=0; i<cols; i++)
      diagp->v[i+off] = fa[i+i*ldf];
  }
  else {
    assert(hm->son);
    assert(hm->rsons == hm->csons);
    
    off1 = off;
    for(i=0; i<hm->csons; i++) {
      getdiag(hm->son[i+i*hm->rsons], diagp, off1);

      off1 += hm->son[i+i*hm->rsons]->rc->size;
    }
    assert(off1 == off + hm->rc->size);
  }
}

void
getdiag_nopermute_hmatrix(pchmatrix hm, pavector diagp)
{
  getdiag(hm, diagp, 0);
}
