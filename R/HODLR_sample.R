#' HODLR_TP_sample
#'
#' Assuming a tensor product f*g, this function samples f given g 
#' using a Gibbs sampler defined in Moran and Wheeler. 

#'
#' @param Y The observed data. 
#'
#' @param TP This parameter represents g. In practice, this value can represent
#' multiple dimensions. 
#' @param x  Values of the covariate associated with f. 
#' @param l_scale The length scale parameter 
#' @param idx   This optional parameter is the index of the x sorted in increasing order. 
#'              If it is not specified, the function will sort x for the Gibbs sampler. 
#'              However, since most samplers only require x to be sorted once, 
#'              this is not efficient.   We recommend passing the index (i.e. order(x)). 
#' @param kernel  The kernel function to use. By default the "matern" kernel is used but the 
#'                squared exponential kernel "sqr-exp" is also avaialble. 
#' @export
#'
#' @examples
#' add_numbers(1, 2) ## returns 3
#'
HODLR_TP_sample <- function(Y,TP,x,cur_tau,l_scale,idx = NULL, kernel = "matern"){
       
        if (is.null(idx)){
                idx = order(x)
        }
        
        tY       <- as.matrix(Y[idx])
        W        <- matrix(TP[idx]^2*cur_tau)
        if (kernel == "matern"){
                Q        <- setup_compressedMatrixGP_Matern(as.matrix(x[idx]),1,l_scale,1e-14,50)
                QP1      <- setup_compressedMatrixGP_Matern_tP(as.matrix(x[idx]),W,1,l_scale,1e-14,50)
        }else if (kernel == "sqr-exp"){
                Q        <- setup_compressedMatrixGP_sqrExp(as.matrix(x[idx]),1,l_scale,1e-14,50)
                QP1      <- setup_compressedMatrixGP_sqrExp_P1(as.matrix(x[idx]),W,1,l_scale,1e-14,50)
        }else{
                stop("Undefined kernel specified.")            
        }
        
        A1       <- simulate_compressedMatrixGP(length(tY),Q)
        A2       <- matmatProduct(Q,matrix(rnorm(length(tY))*sqrt(W),ncol=1))
        R        <- matrix(v1*TP[idx]*cur_tau)
        mTY      <- as.matrix(tY*R)

        T1       <- solve_HODLR(QP1,(1/W)*mTY)
        retval = x 
        retval[idx]<- matmatProduct(Q,T1) +  (1/W)*solve_HODLR(QP1,A1+A2)

        return(retval)

}
