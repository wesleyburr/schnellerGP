
#ifndef _GET_DIAG_H
#define _GET_DIAG_H
#include "avector.h"
#include "hmatrix.h"

HEADER_PREFIX void
getdiag_hmatrix(pchmatrix hm, pavector diag);


HEADER_PREFIX void
getdiag_nopermute_hmatrix(pchmatrix hm, pavector diagp);

#endif 