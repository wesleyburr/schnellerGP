########################################################
#' Given a specified kernel (currently supporting matern and squared Exponential),
#' and a matrix of locations X, set up a compressed matrix Gaussian process
#' pointer via h2lib.
#'
#' @param X Matrix of locations to use.
#' @param kernel_type String specifying which kernel to use. Defaults to matern; currently supported are 'matern' and 'squared_exponential'.
#' @param kernel_params Named list of parameters for the chosen kernel. For supported kernels at present, must contain both 'sigma' and 'rho'.
#' @param tol Specified tolerance for accuracy of calculations.
#' @param M Maximum submatrix size in the H^2 decomposition.
#'
#' @return ptr Pointer to the compressed approximation matrix Gaussian Process.
#'
########################################################

setup_compressedMatrixGP = function(X,
                                    kernel_type = "matern",
                                    kernel_params,
                                    tol,
                                    M)
{
  stopifnot(kernel_type %in% c("matern", "squared_exponential"))
  stopifnot(is.list(kernel_params))
  if(kernel_type %in% c("matern", "squared_exponential")) {
    if( !all( c("sigma", "rho") %in% names(kernel_params) ) ) {
      stop("The kernel_params must be a named list containing definitions for 'sigma' and 'rho'.")
    }
  }

  if( kernel_type == "matern" ) {

    ptr = .setup_compressedMatrixGP_Matern(X,
                                           kernel_params$sigma,
                                           kernel_params$rho,
                                           tol,
                                           M);

  } else if( kernel_type == "squared_exponential" ) {

    ptr = .setup_compressedMatrixGP_sqrExp(X,
                                           kernel_params$sigma,
                                           kernel_params$rho,
                                           tol,
                                           M);

  } else {  # this case cannot be reached because of the stopifnot() above.
    stop("The specified kernel_type is not supported.")
  }

  return(ptr)
}
