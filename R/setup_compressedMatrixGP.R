########################################################
# setup_compressedMatrixGP
# @X - matrix of locations
# @kernel_type  - string specifying which kernel to use
# @kernel_params - list of parameters for the chosen kernel
# @tol - specified tolerance for accuracy of calculations
# @M  - Max submatrix size
########################################################
setup_compressedMatrixGP = function(X, kernel_type, kernel_params, tol, M)
{
  names_kernParams = names(kernel_params)
  if( kernel_type == "squaredExponential" ){ # If user chooses squared exponential kernel.
    if( ("sigma" %in% names_kernParams) & ("rho" %in% names_kernParams) ){
      ptr = .setup_compressedMatrixGP_sqrExp(X, kernel_params$sigma, kernel_params$rho, tol, M);
    } else{
      print('For squared exponential kernel, kernel_params must be a list defining "sigma" and "rho".')
      return(NULL)
    }
  } else if( kernel_type == "matern" ){ # If user chooses Matern kernel.
    if( ("sigma" %in% names_kernParams) & ("rho" %in% names_kernParams) ){
      ptr = .setup_compressedMatrixGP_Matern(X, kernel_params$sigma, kernel_params$rho, tol, M);
    } else{
      print('For squared exponential kernel, kernel_params must be a list defining "sigma" and "rho".')
      return(NULL)
    }
  } else{
    print("Choose one of kernel_type entries currently supported: 'squaredExponential' or 'matern'.")
    return(NULL)
  }
  return(ptr)
}