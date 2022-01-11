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
#' @cur_tau  Current precision estimate of Y, i.e. Y = f\times g + \epsilon.
#'           where \epsilon \sim N(0,\tau^-1).
#' @param l_scale The length scale parameter 
#' @param idx   This optional parameter is the index of the x sorted in increasing order. 
#'              If it is not specified, the function will sort x for the Gibbs sampler. 
#'              However, since most samplers only require x to be sorted once, 
#'              this is not efficient.   We recommend passing the index (i.e. order(x)). 
#' @param kernel  The kernel function to use. By default the "matern" kernel is used but the 
#'                squared exponential kernel "sqr-exp" is also available. 
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
        R        <- matrix(TP[idx]*cur_tau)
        mTY      <- as.matrix(tY*R)

        T1       <- solve_HODLR(QP1,(1/W)*mTY)
        retval = x 
        retval[idx]<- matmatProduct(Q,T1) +  (1/W)*solve_HODLR(QP1,A1+A2)

        return(retval)

}

#' HODLR_TP_sample_ls
#'
#' Assuming a tensor product f*g, this function samples the length-scale
#' for f, given the marginal distribution of Y given g. It performs
#' a simple Metropolis step where the proposal distribution is 
#' a uniform(-0.5,0.5) and the prior distribution on the length-scale
#' is a uniform distribution specified by the bounds argument. 
#'
#' @param ls The current value of the length scale being used for f. 
#' @param Y The observed data or the residual given other components in
#' the mean. 
#' @param TP This parameter represents g. In practice, this value can represent
#' multiple dimensions. 
#' @param x  Values of the covariate associated with f.
#' @cur_tau  Current precision estimate of Y, i.e. Y = f\times g + \epsilon.
#'           where \epsilon \sim N(0,\tau^-1).  
#' @param idx   This optional parameter is the index of the x sorted in increasing order. 
#'              If it is not specified, the function will sort x for the Gibbs sampler. 
#'              However, since most samplers only require x to be sorted once, 
#'              this is not efficient.   We recommend passing the index (i.e. order(x)). 
#' @param kernel  The kernel function to use. By default the "matern" kernel is used but the 
#'                squared exponential kernel "sqr-exp" is also available. 
#' @export
#'
#' @examples
#' add_numbers(1, 2) ## returns 3
#'
HODLR_TP_sample_ls <- function(ls,Y,TP,x,cur_tau,bounds,idx = NULL,kernel="matern"){
        
        nls <- runif(1,-0.5,0.5)+ls
        if ( (nls < bounds[1]) || ( nls > bounds[2]) ){
                return(ls)
        }
        if (is.null(idx)){
                idx = order(x)
        }
        
        ny    <- as.matrix(t(Y/(TP)))[idx]
        W     <- as.matrix(1/(TP)^2*(1/cur_tau))[idx]
        if (kernel == "matern"){
                QP1   <- setup_compressedMatrixGP_Matern_tP(as.matrix(x[idx]),as.matrix(1/W),1,ls,1e-14,50)
                QP2   <- setup_compressedMatrixGP_Matern_tP(as.matrix(x[idx]),as.matrix(1/W),1,nls,1e-14,50)
        }else if (kernel == "sqr-exp"){
                function(X, sigma, rho, tol, M)
                QP1   <- setup_compressedMatrixGP_sqrExp_tP(as.matrix(x[idx]),as.matrix(1/W),1,ls,1e-14,50)
                QP2   <- setup_compressedMatrixGP_sqrExp_tP(as.matrix(x[idx]),as.matrix(1/W),1,nls,1e-14,50)
                
        }else{
                stop("Undefined kernel specified.") 
        }
        
        denom_like    <- sum(solve_HODLR(QP1,as.matrix(ny))*ny)
        denom_det     <- log_det_HODLR(QP1)
        num_like      <- sum(solve_HODLR(QP2,as.matrix(ny))*ny)
        num_det       <- log_det_HODLR(QP2)
        test = 0.5*(-num_like - num_det) - 0.5*( -denom_like - denom_det)
        if (exp(test) > runif(1))
        {
                return(nls)
        }else{
                return(ls)
        }
        
}

#' HODLR_TP_sample_cale
#'
#' Assuming a tensor product of k 1-D functions, i.e. \gamma f_1 \cdot f_2\cdots f_k,
#'  this function samples \gamma, which is the scale parameter. 
#'  Here \gamma \sim N(\mu,\alpha^-1), and by default \mu = 0 and \alpha =1.
#'
#' @param Y The observed n \times 1 data vector or the residual given \
#' other components in the mean. 
#' @param TP A matrix of dimension n \times k. Here, each row corresponds
#' to an observation in Y and each column corresponds to one of the k dimensions 
#' in the tensor Product. 
#' @cur_tau  Current precision estimate of Y, i.e. Y = \gamma f_1 \cdot f_2\cdots f_k + \epsilon.
#'           where \epsilon \sim N(0,\tau^-1).  
#' @export
#'
#' @examples
#' add_numbers(1, 2) ## returns 3
HODLR_TP_sample_scale <- function(Y,TP,cur_tau,mean=0,prec=1){
        if (length(Y)!=nrow(TP)){
                stop("The data and the TP dimensions do not match.")
        }
        
        tY = as.matrix(Y)
        tX = as.matrix(apply(TP,1,prod))
        tV = 1/(t(tX)%*%tX * cur_tau + prec)
        tM = tV*(t(tX)%*%(cur_tau*Y) + mean*prec)
        return(rnorm(1,tM,sqrt(tV)))

}
