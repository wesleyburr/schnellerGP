#include "mclasses/kernel_functions.h"

//////////////////////////////////////////////////////
//Function: kernel_sqr_exp()
//Purpose : Produce the (i,j) entry of a covariance matrix
//          using the squared exponential kernel.
//          Nugget is 1e-8. Length-Scale is parameter 0 in gpd->hparm
 field
kernel_sqr_exp(const double *xx, const double *yy, void *data)
{
  double norm2 = 0.0;

  gp_data* gpd =  (gp_data*) data;

  if (gpd->i == gpd->j){
    return gpd->hparm[0] + 1e-7; // diagonal add 1e-8 as the nugget
  }

  for(int i = 0; i < gpd->d; i++){
    norm2 += (xx[i] - yy[i])*(xx[i] - yy[i]);
  }
  // for the squared exponential term we have hparm[0] as the length scale
  // hparm[1] is the variance component
  return gpd->hparm[1]*REAL_EXP(-gpd->hparm[0]*norm2); // return exp(-c||X_i-X_j||^2_2)

}

//Function: kernel_sqr_exp_truncate()
//Purpose : Produce the (i,j) entry of a covariance matrix
//          using the squared exponential kernel.
//          Nugget is 1e-8. Length-Scale is parameter 0 in gpd->hparm
//          Truncates values less than 1e-7.
field
kernel_sqr_exp_truncate(const double *xx, const double *yy, void *data)
  {
    double norm2 = 0.0;

    gp_data* gpd =  (gp_data*) data;

    if (gpd->i == gpd->j){
      return gpd->hparm[0] + 1e-7; // diagonal add 1e-8 as the nugget
    }

    for(int i = 0; i < gpd->d; i++){
      norm2 += (xx[i] - yy[i])*(xx[i] - yy[i]);
    }
    // for the squared exponential term we have hparm[0] as the length scale
    // hparm[1] is the variance component
    field temp = gpd->hparm[1]*REAL_EXP(-gpd->hparm[0]*norm2);
    return fabs(temp) < 1e-7? 0.0: temp; // return exp(-c||X_i-X_j||^2_2) if > 1e-7, 0 otherwise

  }

//////////////////////////////////////////////////////
//Function: kernel_diag()
//Purpose : Produce the (i,j) entry of a covariance matrix
//          using the diagional kernel. Here every non-diagonal
//          entry is 0
field
kernel_diagonal(const double *xx, const double *yy, void *data)
{
  double norm2 = 0.0;

  gp_data* gpd =  (gp_data*) data;

  if (gpd->i == gpd->j){
    return xx[0] ; // X represents the diagonal elements
  }else{
    return 0;
  }

}

//////////////////////////////////////////////////////
//Function: kernel_matern_exp()
//Purpose : Produce the (i,j) entry of a covariance matrix
//          using the matern exponential kernel.
field
kernel_matern_exp(const double *xx, const double *yy, void *data)
  {
    double norm2 = 0.0;

    gp_data* gpd =  (gp_data*) data;

    if (gpd->i == gpd->j){
      return gpd->hparm[1] + 1e-7; // diagonal add 1e-8 as the nugget
    }

    for(int i = 0; i < gpd->d; i++){
      norm2 += (xx[i] - yy[i])*(xx[i] - yy[i]);
    }
    norm2 = sqrt(norm2)/gpd->hparm[0];
    // for the squared exponential term we have hparm[0] as the length scale
    // hparm[1] is the variance component
    // return gpd->hparm[0] *

    return gpd->hparm[1]*(1 + sqrt(5) * norm2 + (5/3) * (norm2 * norm2)) * exp(-sqrt(5) * norm2); // return exp(-c||X_i-X_j||^2_2)

  }
