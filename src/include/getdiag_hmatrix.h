/** @brief Copy the diagonal of an @ref hmatrix, i.e., fill a vector
 *    @f$d_i \gets g_{ii}@f$.
 *
 *  This function uses the original index order, i.e., the vector
 *  has to be large enough to accomodate indices from <tt>hm->rc->idx</tt>
 *  and <tt>hm->cc->idx</tt>.
 *
 *  @param hm Source matrix with <tt>hm->rc==hm->cc</tt>.
 *  @param diag Target vector. */
HEADER_PREFIX void
getdiag_hmatrix(pchmatrix hm, pavector diag);

/** @brief Copy the diagonal of an @ref hmatrix, i.e., fill a vector
 *    @f$d_i \gets g_{ii}@f$.
 *
 *  This function uses the internal index order defined by the
 *  cluster tree, i.e., it can be called for submatrices and
 *  the vector only has to have a dimension of <tt>hm->rc->size</tt>.
 *
 *  @param hm Source matrix with <tt>hm->rc==hm->cc</tt>.
 *  @param diagp Target vector. */
HEADER_PREFIX void
getdiag_nopermute_hmatrix(pchmatrix hm, pavector diagp);