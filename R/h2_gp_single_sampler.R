########################################################
#' sample_gp
#'
#' @param h2 H^2 approximation matrix to the covariance matrix.
#' @param h Cholesky decomposition of the covariance matrix
#' @param tX Regression coefficients or weights
#' @param sigma_p Precision of the variance-covariance matrix
#' @param tau Data precision.
#' @importFrom stats rnorm
#'
#'
########################################################
sample_gp <- function(h2,
                      h,
                      y,
                      tX,
                      sigma_p,
                      tau)
{
  tW = t(tX) %*% tX * tau;
  rw = diag(tW)
  idx = get_h2_idx(h2)
  n1 = matrix(data = rnorm(n = length(idx)))
  n2 = matrix(data = rnorm(n = length(idx),
                           mean = 0,
                           sd = sqrt((1/sigma_p) * rw)
                           )
              )

  diag(tW) = sigma_p / rw
  n1[idx]  = sqrt(1/sigma_p) * multiply_hmat_vec(h, n1)
  n2       = sqrt(1/sigma_p) * multiply_h2_vector(h2, as.matrix(n2))

  add_diag(h2, sigma_p / rw)
  n3 = solve_h2_vector(h2,
                       matrix(data = sigma_p/rw),
                       n1 + n2,
                       1e-7)
  n4 = solve_h2_vector(h2,
                       matrix(data = 1 / rw),
                       matrix(data = tW %*% t(tX * tau) %*% y),
                       1e-7)
  add_diag(h2, -sigma_p / rw)
  n3 = sigma_p / rw * n3
  n4 = (1 / sigma_p) * multiply_h2_vector(h2, as.matrix(n4))
  betas = n4 + n3
  return(betas)
}
