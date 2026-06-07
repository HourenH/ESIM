## last revised: 30.Sep.2025
## estimation function for ESIM
## code for SIQR is adapted from Wu et al.(2010); 
## code for FSIM is adapted from Bhattacharjee and M\{''}uller (2023). AoS. See also https://github.com/functionaldata/tFrechet
## require packages: 
## - 
library(KernSmooth)
library(quantreg)
library(frechet)

###########################################################
# The implementation of ESIM with LS
###########################################################
#' Normalizing a vector to the unit vector
#'
#' @param x A numeric vector of length p.
#'
#' @returns A p-dimensional vector of unit norm.
#' @export
#'
#' @examples
SpheNormalize <- function(x){
    x / l2norm(x)
}

#' Compute extrinsic local linear regression using the least squares.
#'
#' @param txdat A numeric vector of length n, the index value \eqn{theta^\top X}.
#' @param tydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param exdat A numeric scaler, the evaluation points.
#' @param bw A positive numeric number, the bandwidth for the local linear regression.
#'
#' @returns A numeric vector of length 2d. 
#' The leading d components are the estimate of the local linear regression at exdat, and the last d components the estimate of its first derivative.
#' 
#' @export
#'
#' @examples
ll_LS <- function(txdat, tydat, exdat, bw){
    K = dnorm((txdat - exdat)/bw)/bw 
    s0 = mean(K)
    s1 = mean(K * (txdat - exdat))
    s2 = mean(K * (txdat - exdat)^2)
    
    mu_weight = ((s2 - s1*(txdat - exdat))/(s2*s0 - s1^2)) * K
    mudev_weight = ((s1 - s0*(txdat - exdat))/(s1^2-s2*s0)) * K
    
    mu = colMeans(tydat * matrix(rep(mu_weight, ncol(tydat)), ncol = ncol(tydat)))
    mudev = colMeans(tydat * matrix(rep(mudev_weight, ncol(tydat)), ncol = ncol(tydat)))
    return(c(mu,mudev))
}

#' Compute weighted extrinsic local linear regression using least squares.
#'
#' @param txdat A numeric vector of length n, the index value \eqn{theta^\top X}.
#' @param tydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param exdat A numeric scaler, the evaluation point.
#' @param bw A positive numeric number, the bandwidth for the local linear regression.
#' @param weight A numeric vector of length n, the observation weights.
#'
#' @returns A numeric vector of length 2d.
#' The leading d components are the estimate of the local linear regression at exdat, and the last d components the estimate of its first derivative.
#'
#' @export
#'
#' @examples
ll_LS_weight <- function(txdat, tydat, exdat, bw, weight){
    K = dnorm((txdat - exdat)/bw) * weight / bw 
    s0 = mean(K)
    s1 = mean(K * (txdat - exdat))
    s2 = mean(K * (txdat - exdat)^2)
    
    mu_weight = ((s2 - s1*(txdat - exdat))/(s2*s0 - s1^2)) * K
    mudev_weight = ((s1 - s0*(txdat - exdat))/(s1^2-s2*s0)) * K
    
    mu = colMeans(tydat * matrix(rep(mu_weight, ncol(tydat)), ncol = ncol(tydat)))
    mudev = colMeans(tydat * matrix(rep(mudev_weight, ncol(tydat)), ncol = ncol(tydat)))
    return(c(mu,mudev))
}

#' Cost function of the extrinsic single index model using the least squares.
#'
#' @param param A numeric vector of length (p+1), the first component is the bandwidth and the remaining components are the parameter coefficients.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' 
#' @returns Scalar, the loss value.
cost_ESIM <- function(param, xdat, ydat){
    MaxPanelty = .Machine$double.xmax
    K = dnorm # kernel function
    dx = ncol(xdat); dy = ncol(ydat)
    n = nrow(xdat)
    
    if(dx != (length(param)-1)){
        stop("xdat has different dimensions with param.")
    }
    
    bw = param[1] # bandwidth
    
    theta = param[-1]
    theta = sign(theta[1]) * SpheNormalize(theta) # parametric coefficients
    
    index_ = xdat %*% theta
    
    bw_range = SetBwRange(xin = index_, xout = index_, kernel_type = "gauss")# avoid small or negative bandwidth
    if(bw>bw_range$max || bw < bw_range$min){return(.Machine$double.xmax)}
    
    # leave one out local linear of mu and first derivative
    llmu = sapply(1:n, function(r_){
        ll_LS(txdat = index_[-r_], tydat = ydat[-r_,], exdat = index_[r_], bw)
    })
    llmu = t(llmu)[,1:dy]
    # objective value
    if (any(is.nan(llmu))) {return(.Machine$double.xmax)}
    sum((ydat - llmu)^2)
    
}

#' Estimation of the extrinsic single-index model using the least squares.
#'
#' This function estimates the parameter coefficients \eqn{\theta} and the bandwidth 
#' for an extrinsic single-index model using a leave-one-out local linear regression 
#' with the least squares loss.
#' 
#' @param param A numeric of length (p+1). The first component is the initial value of the bandwidth and the remaining components are the initial value of the parametric coefficients.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#'
#' @returns A list containing:
#' \describe{
#'   \item{theta}{Estimated parameter coefficients, normalized to have unit length and a positive leading component.}
#'   \item{bw}{Estimated bandwidth for local linear regression.}
#'   \item{mu}{An n by d matrix of leave-one-out local linear regression estimates of the response.}
#'   \item{mudev}{An n by d matrix of leave-one-out estimates of the derivative of the regression function.}
#' }
#' 
#' @export
#'
#' @examples
ls_est <- function(param, xdat, ydat){
    dy = ncol(ydat)
    est = optim(param, cost_ESIM, xdat = xdat, ydat = ydat)
    
    
    bw = est$par[1]
    theta = est$par[-1]
    theta = sign(theta[1]) * SpheNormalize(theta)
    
    index = xdat %*% theta
    llmu = t(sapply(1:nrow(xdat), function(r_){
        ll_LS(txdat = index[-r_], tydat = ydat[-r_,], exdat = index[r_], bw)}))
    mu = llmu[,1:dy]
    mudev = llmu[,(dy+1):(2*dy)]
    
    return(list(theta = theta, bw = bw, mu = mu, mudev = mudev))
}

#' Weighted cost function of the extrinsic single-index model using least squares.
#'
#' @param param A numeric vector of length (p+1), the first component is the bandwidth and the remaining components are the parameter coefficients.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param weight A numeric vector of length n, the observation weights for the weighted least squares objective.
#'
#' @returns Scalar, the weighted loss value.
cost_ESIM_weight<- function(param, xdat, ydat, weight){
    MaxPanelty = .Machine$double.xmax
    K = dnorm # kernel function
    dx = ncol(xdat); dy = ncol(ydat)
    n = nrow(xdat)

    if(dx != (length(param)-1)){
        stop("xdat has different dimensions with param.")
    }

    bw = param[1] # bandwidth

    theta = param[-1]
    theta = sign(theta[1]) * SpheNormalize(theta) # parametric coefficients

    index_ = xdat %*% theta

    bw_range = SetBwRange(xin = index_, xout = index_, kernel_type = "gauss")# avoid small or negative bandwidth
    if(bw>bw_range$max || bw < bw_range$min){return(.Machine$double.xmax)}

    # leave one out local linear of mu and first derivative
    llmu = sapply(1:n, function(r_){
        ll_LS_weight(txdat = index_[-r_], tydat = ydat[-r_,], exdat = index_[r_], bw, weight = weight[-r_])
    })
    llmu = t(llmu)[,1:dy]
    # objective value
    if (any(is.nan(llmu))) {return(.Machine$double.xmax)}
    sum(weight * (ydat - llmu)^2)

}

#' Weighted estimation of the extrinsic single-index model using least squares.
#'
#' This function estimates the parameter coefficients \eqn{\theta} and the bandwidth
#' for an extrinsic single-index model using a weighted leave-one-out local linear regression
#' with the least squares loss.
#'
#' @param param A numeric of length (p+1). The first component is the initial value of the bandwidth and the remaining components are the initial value of the parametric coefficients.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param weight A numeric vector of length n, the observation weights for the weighted least squares objective.
#'
#' @returns A list containing:
#' \describe{
#'   \item{theta}{Estimated parameter coefficients, normalized to have unit length and a positive leading component.}
#'   \item{bw}{Estimated bandwidth for local linear regression.}
#'   \item{mu}{An n by d matrix of leave-one-out weighted local linear regression estimates of the response.}
#'   \item{mudev}{An n by d matrix of leave-one-out estimates of the derivative of the regression function.}
#' }
#'
#' @export
#'
#' @examples
ls_est_weight <- function(param, xdat, ydat, weight){
    dy = ncol(ydat)
    est = optim(param, cost_ESIM_weight, xdat = xdat, ydat = ydat, weight = weight)


    bw = est$par[1]
    theta = est$par[-1]
    theta = sign(theta[1]) * SpheNormalize(theta)

    index = xdat %*% theta
    llmu = t(sapply(1:nrow(xdat), function(r_){
        ll_LS_weight(txdat = index[-r_], tydat = ydat[-r_,], exdat = index[r_], bw, weight = weight[-r_])}))
    mu = llmu[,1:dy]
    mudev = llmu[,(dy+1):(2*dy)]

    return(list(theta = theta, bw = bw, mu = mu, mudev = mudev))
}

###########################################################
# The implementation of ESIM with the exponential squared loss
###########################################################
#' Compute extrinsic local linear regression using the exponential squared loss. Only used when interatively estimating theta in function `cost_index_bw'
#'
#' @param txdat A numeric vector of length n, the index value \eqn{theta^\top X}.
#' @param tydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param exdat A numeric scaler, the evaluation points.
#' @param bw A positive numeric number, the bandwidth for the local linear regression.
#' @param weight A numeric vector of length n, the weights for the local linear regresson.
#'
#' @returns A numeric vector of length 2d. 
#' The leading d components are the estimate of the local linear regression at exdat, and the last d components the estimate of its first derivative.
#' 
#' @export
#'
#' @examples
ll_ESL <- function(txdat, tydat, exdat, bw, weight){
    K = dnorm((txdat - exdat)/bw)/bw 
    s0 = mean(weight * K)
    s1 = mean(weight * K * (txdat - exdat))
    s2 = mean(weight * K * (txdat - exdat)^2)
    
    mu_weight = ((s2 - s1*(txdat - exdat))/(s2*s0 - s1^2)) * weight * K
    mudev_weight = ((s1 - s0*(txdat - exdat))/(s1^2-s2*s0)) * weight * K
    
    mu = colMeans(tydat * matrix(rep(mu_weight, ncol(tydat)), ncol = ncol(tydat)))
    mudev = colMeans(tydat * matrix(rep(mudev_weight, ncol(tydat)), ncol = ncol(tydat)))
    return(c(mu,mudev))
}

#' Local linear regression with exponential squared loss (prediction)
#'
#' @param txdat Numeric vector of length n. Training index values \eqn{\hat{\theta}^\top x} 
#' @param tydat An n by d matrix of unit responses.
#' @param exdat Numeric scalar. Target index value at which prediction is made.
#' @param bw Positive numeric scalar. Bandwidth for kernel weighting.
#' @param lambda Positive numeric scalar. Tuning parameter for the exponential squared loss.
#' @param init_param Numeric vector of length \eqn{2d}. Initial parameter values 
#' for optimization, consisting of \eqn{d} entries for the conditional mean and 
#' \eqn{d} entries for its derivative.
#' @param w Numeric vector of length n. Weights for the loss function. Equal weights as default
#'
#' @returns A list with components:
#' \item{mu}{Estimated conditional mean vector of length \eqn{d}.}
#' \item{mu_deriv}{Estimated derivative vector of length \eqn{d}.}

ll_ESL_predict <- function(txdat, tydat, exdat, bw, lambda, init_param, w = NA){
    U = txdat - exdat
    K = dnorm((U)/bw)/bw 
    dy = ncol(tydat)
    n = nrow(tydat)
    
    if (anyNA(w)){w = rep(1, n)}
    
    objFctn <- function(param){
        a = param[1:dy]; b = param[(dy+1):(2*dy)]
        diff = tydat - matrix(a, n, dy, byrow = T) - U %*% t(b)
        mean(w * (1 - exp(- rowSums(diff^2)/lambda)) * K)
    }
    
    res = optim(init_param, objFctn)
    
    list(mu = res$par[1:dy], mu_deriv = res$par[(dy+1):(2*dy)])
}

#' The cost function for estimating the tuning parameter \eqn{\lambda} when using the exponential squared loss. 
#'
#' @param lambda Scalar, the value of tuning parameter.
#' @param delta Scalar, the prespecified value controling robustness. See Section 5 of the paper for details. 
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param mu_fit An n by d matrix, each row a d-dimensional unit vector obtained from the leave-one-out local linear regression using the ESL.
#'
#' @returns Scalar, the loss value.
cost_lambda_median <- function(lambda, delta, ydat, mu_fit){
    e = ydat - mu_fit
    abs(delta - 1 + median(exp(-rowSums(e^2)/lambda)))
}

# Update all leave-one-out IRLS weight columns for fixed local fits.
update_esl_psi <- function(index, ydat, mu, mudev, lambda, w){
    n = nrow(ydat); d = ncol(ydat)

    sapply(1:n, function(r_){
        U = index[-r_] - index[r_]
        fitted = matrix(mu[r_,], nrow = n-1, ncol = d, byrow = TRUE) +
            U %o% mudev[r_,]
        e = ydat[-r_,,drop = FALSE] - fitted
        exp(-rowSums(e^2)/lambda) * w[-r_]
    })
}

#' The weighted cost function for estimating \eqn{\theta} and the bandwidth, using the exponential squared loss.
#'
#' @param param A numeric vector of length (p+1), the initial value for the bandwidth and \eqn{\theta} (p dimensional vector).
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param lambda A positive scalar, the tuning parameter that controls the robustness.
#' @param weight An (n-1) by n matrix, weights for ESL computation. 
#' @param w A numeric vector of length n, weighted loss. Equal weights as default
#'
#' @returns A scalar, the loss value.
cost_index_bw <- function(param, xdat, ydat, lambda, weight, w=NA){
    theta = SpheNormalize(param[-1])
    theta = sign(theta[1]) * theta
    bw = param[1]
    
    if(anyNA(w)){w = rep(1, nrow(ydat))}
    
    index_ = xdat %*% theta
    
    bw_range = SetBwRange(xin = index_, xout = index_, kernel_type = "gauss")
    if(bw>bw_range$max || bw < bw_range$min){return(.Machine$double.xmax)}
    
    # leave one out local linear of mu and first derivative
    llmu = sapply(1:nrow(ydat), function(r_){
        ll_ESL(txdat = index_[-r_], tydat = ydat[-r_,], exdat = index_[r_], bw, weight = weight[, r_])
    })
    llmu = t(llmu)
    # objective value
    if (any(!is.finite(llmu))) {return(.Machine$double.xmax)}
    mean(w* (1 - exp(- rowSums((ydat - llmu[,1:ncol(ydat)])^2)/lambda)))
}

#' Iterative estimation of the extrinsic single-index model using exponential squared loss
#'
#' This function estimates the coefficient vector (index parameter) and the bandwidth 
#' of the extrinsic single-index model (ESIM) under the exponential squared loss (ESL). 
#' The tuning parameter \code{lambda} is treated as fixed.
#' 
#' @param theta_init Numeric vector of length \eqn{p}. Initial value of the coefficient 
#' vector (index parameter).
#' @param bw_init A positive scalar, the initial value of the bandwidth.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param lambda Positive numeric scalar. Fixed tuning parameter in the exponential squared loss.
#' @param w Numerical vector of length n, weighted loss with equal weights by default.
#' @param abstol A small positive scalar, the stopping criterion for the iteration based on parameter changes. Default is 1e-4.
#' @param maxiter An integer, the maximum number of iterations allowed. Default is 100.
#'
#' @returns A list containing:
#' \item{theta}{Estimated coefficient vector (unit-norm, with positive first component).}
#' \item{bw}{Estimated bandwidth parameter.}
#' \item{mu}{Matrix of estimated conditional means \eqn{\hat{\mu}(\theta^\top x_i)} for each observation.}
#' \item{mudev}{Matrix of estimated first-order derivatives of the regression function.}
#' \item{psi}{Matrix of estimated weights from the exponential squared loss.}
#' \item{iter}{Number of iterations used until convergence (or \code{maxiter} if not converged).}

esl_index_bw <- function(theta_init, bw_init, xdat, ydat, lambda, w = NA, abstol=1e-4, maxiter=200){
    dist_ = 1e5
    iter = 0
    n = nrow(xdat); d = ncol(ydat)
    
    if (anyNA(w)) {w = rep(1, n)}
    psi = sapply(1:n, function(r_){w[-r_]})
    
    theta_update = sign(theta_init[1]) * SpheNormalize(theta_init) # normalize
    bw_update = bw_init
    
    while (dist_ >= abstol && iter < maxiter) {
        theta_old = theta_update 
        
        # given theta and bw, updates mu and mudev
        index_update = xdat %*% theta_update # n * 1
        llmu_mudev = t(sapply(1:n, function(r_){
            ll_ESL(txdat = index_update[-r_], tydat = ydat[-r_,], exdat = index_update[r_],
                   bw = bw_update, weight = psi[,r_])})) # n * 2d
        mu_update = llmu_mudev[,1:d]
        mudev_update = llmu_mudev[,(d+1):(2*d)]
        
        # given lambda, mu_update, update weight
        psi = update_esl_psi(index = index_update, ydat = ydat, mu = mu_update,
                             mudev = mudev_update, lambda = lambda, w = w)
        
        # given psi, update theta and bw
        para_update = optim(c(bw_update, theta_update), cost_index_bw, 
                            xdat=xdat, ydat=ydat, lambda=lambda, weight= psi, w = w)
        theta_update = para_update$par[-1]
        theta_update = sign(theta_update[1]) * SpheNormalize(theta_update)
        bw_update = para_update$par[1]
        
        dist_ = SpheGeoDist(theta_update, theta_old)
        iter = iter + 1 # update counter
    }
    if (iter == maxiter) message("ESL fail to converge.")
    
    # given theta and bw, find mu and mudev numerically
    index_update = xdat %*% theta_update # n * 1
    
    mu_mudev = t(sapply(1:n, function(r_){
        fit = ll_ESL_predict(txdat = index_update[-r_], tydat = ydat[-r_,], exdat = index_update[r_],
                                      bw = bw_update, lambda = lambda, init_param = llmu_mudev[r_,], w = w[-r_])
        c(fit$mu, fit$mu_deriv)
    }))
    
    mu_update = mu_mudev[,1:d]
    mudev_update = mu_mudev[,(d+1):(2*d)]
    psi = update_esl_psi(index = index_update, ydat = ydat, mu = mu_update,
                         mudev = mudev_update, lambda = lambda, w = w)
    if (any(is.na(mu_update))) message("NA detected in mu.")
    return(list(theta = theta_update, bw = bw_update, 
                mu = mu_update, mudev = mudev_update, psi = psi, iter = iter))
}


#' Estimation of extrinsic single-index model using the exponential squared loss
#' 
#' The function implements Algorithm 1 by iteratively updating the tuning parameter
#' \eqn{\lambda}, the leave-one-out IRLS weights, the index coefficient, and the bandwidth.
#'  
#' @param theta_init A numeric vector of length p, the initial value of \eqn{\theta}.
#' @param bw_init A positive scalar, the initial value of the bandwidth.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param delta A scalar within (0,1) controling the robustness. See Section 5 of the paper for details.
#' @param abstol A small positive scalar, the stopping criterion based on changes in the index coefficient. Default is 1e-4.
#' @param maxiter An integer, the maximum number of Algorithm 1 iterations. Default is 200.
#'
#' @returns A list containing the following elements:
#' \describe{
#'   \item{theta}{The estimated p-dimensional parameter vector \eqn{\theta}.}
#'   \item{bw}{The estimated bandwidth.}
#'   \item{mu}{An n by d matrix of the estimated local linear regression values.}
#'   \item{mudev}{An n by d matrix of the estimated first derivatives of the local linear regression.}
#'   \item{lambda}{The estimated tuning parameter \eqn{\lambda}.}
#'   \item{psi}{An (n-1) by n matrix of weights used in the local linear regression.}
#'   \item{iter}{Number of Algorithm 1 iterations used.}
#' }
#' 
#' @export
#'
#' @examples
esl_est <- function(theta_init, bw_init, xdat, ydat, delta, abstol=1e-4, maxiter=200){
    dist_ = 1e5
    iter = 0
    n = nrow(ydat); d = ncol(ydat)

    theta_update = sign(theta_init[1]) * SpheNormalize(theta_init) # normalize
    bw_update = bw_init
    w = rep(1, n)
    psi = matrix(1, nrow = n-1, ncol = n)

    while (dist_ >= abstol && iter < maxiter) {
        theta_old = theta_update

        # Step 1: update leave-one-out local fits and lambda
        index_update = xdat %*% theta_update
        llmu_mudev = t(sapply(1:n, function(r_){
            ll_ESL(txdat = index_update[-r_], tydat = ydat[-r_,], exdat = index_update[r_],
                   bw = bw_update, weight = psi[,r_])
        }))
        mu_update = llmu_mudev[,1:d]
        mudev_update = llmu_mudev[,(d+1):(2*d)]
        lambda_update = optimise(cost_lambda_median, lower = 1e-5, upper = 1e5, tol = 1e-10,
                                 delta = delta, ydat = ydat, mu_fit = mu_update)$minimum

        # Step 2: update the leave-one-out IRLS weights
        psi = update_esl_psi(index = index_update, ydat = ydat, mu = mu_update,
                             mudev = mudev_update, lambda = lambda_update, w = w)

        # Step 3: update theta and bandwidth for fixed lambda and IRLS weights
        para_update = optim(c(bw_update, theta_update), cost_index_bw,
                            xdat = xdat, ydat = ydat, lambda = lambda_update, weight = psi)
        theta_update = para_update$par[-1]
        theta_update = sign(theta_update[1]) * SpheNormalize(theta_update)
        bw_update = para_update$par[1]

        dist_ = SpheGeoDist(theta_update, theta_old)
        iter = iter + 1
    }
    if (iter == maxiter) message("ESL fail to converge.")

    # Synchronize Step 1 and Step 2 at the final theta and bandwidth.
    index_update = xdat %*% theta_update
    llmu_mudev = t(sapply(1:n, function(r_){
        ll_ESL(txdat = index_update[-r_], tydat = ydat[-r_,], exdat = index_update[r_],
               bw = bw_update, weight = psi[,r_])
    }))
    mu_update = llmu_mudev[,1:d]
    mudev_update = llmu_mudev[,(d+1):(2*d)]
    lambda_update = optimise(cost_lambda_median, lower = 1e-5, upper = 1e5, tol = 1e-10,
                             delta = delta, ydat = ydat, mu_fit = mu_update)$minimum
    psi = update_esl_psi(index = index_update, ydat = ydat, mu = mu_update,
                         mudev = mudev_update, lambda = lambda_update, w = w)

    list(theta = theta_update, bw = bw_update, mu = mu_update, mudev = mudev_update,
         lambda = lambda_update, psi = psi, iter = iter)
}



# LOO local linear regression to determine the optimal bandwidth
bwConExp <- function(xdat, ydat){
    # given index value (xdat), compute the optimal bandwidth
    n = nrow(ydat); p = ncol(ydat)
    
    objFctn <- function(bw){
        # obtain LOO local linear
        mu = t(sapply(1:n, function(r){
            ll_LS(txdat = xdat[-r], tydat = ydat[-r,], exdat = xdat[r], bw = bw)}))
        mean(rowSums((ydat - mu[,1:p])^2))
    }
    
    bwRange = SetBwRange(xin = xdat, xout = xdat, kernel_type = "gauss")
    result = optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
    
    bw = result$minimum
    fit = t(sapply(1:n, function(r){
        ll_LS(txdat = xdat[-r], tydat = ydat[-r,], exdat = xdat[r], bw = bw)}))
    list(bw = bw, fit = fit)
}
###########################################################
# The implementation of single index Fr\'{e}chet regression 
###########################################################
#' Bandwidth selection using the 5-fold cross-validation for local Fréchet regression on a unit hypersphere. 
#' The implementation is adapted from Bhattacharjee and M\{''}uller (2023). AoS.
#' See also https://github.com/functionaldata/tFrechet
#'
#' @param xin A numeric vector of length n, the input covariates.
#' @param yin An n by d matrix of responses, each row a unit vector on the sphere.
#' @param xout A numeric vector of length k, the evaluation points.
#' @param optns A list of options, including the bandwidth and the kernel function.
#'
#' @returns A numeric value giving the selected bandwidth.
#' @export 
#'
#' @examples
bwCV_sphe <- function(xin, yin, xout, optns) {
    yin <- yin[order(xin),]
    xin <- sort(xin)
    compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5) # omit observations at two boundaries
    
    # k-fold
    objFctn <- function(bw) {
        optns1 <- optns
        optns1$bw <- bw
        folds <- numeric(length(xin))
        n <- sum(compareRange)
        numFolds <- ifelse(n > 30, 5, sum(compareRange))
        
        tmp <- c(sapply(1:ceiling(n/numFolds), function(i)
            sample(x = seq_len(numFolds), size = numFolds, replace = FALSE)))
        tmp <- tmp[1:n]
        repIdx <- which(diff(tmp) == 0)
        for (i in which(diff(tmp) == 0)) {
            s <- tmp[i]
            tmp[i] <- tmp[i-1]
            tmp[i-1] <- s
        }
        
        folds[compareRange] <- tmp
        
        yout <- lapply(seq_len(numFolds), function(foldidx) {
            testidx <- which(folds == foldidx)
            res <- LocSpheGeoReg(xin = xin[-testidx], yin = yin[-testidx,], xout = xin[testidx], optns = optns1)
            res # each row is a spherical vector
        })
        yout <- do.call(rbind, yout)
        yinMatch <- yin[which(compareRange)[order(tmp)],]
        mean(sapply(1:nrow(yout), function(i) SpheGeoDist(yout[i,], yinMatch[i,])^2))
    }
    bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)
    
    res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
    res$minimum
}


#' Set up bandwidth range.
#' The implementation is adapted from https://github.com/functionaldata/tFrechet
#'
#' @param xin A numeric vector of length n, the input covariates.
#' @param xout A numeric vector of length k, the evaluation points.
#' @param kernel_type A character string specifying the kernel function. Currently supports \code{"gauss"} and \code{"gausvar"}.
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{min}{The minimum bandwidth.}
#'   \item{max}{The maximum bandwidth.}
#' }
#' 
#' @export
#'
#' @examples
SetBwRange <- function(xin, xout, kernel_type) {
    xinSt <- unique(sort(xin))
    bw.min <- max(diff(xinSt), xinSt[2] - min(xout), max(xout) - xinSt[length(xinSt)-1])*1.1 / (ifelse(kernel_type == "gauss", 3, 1) *
                                                       ifelse(kernel_type == "gausvar", 2.5, 1))
    bw.max <- diff(range(xin))/3
    if (bw.max < bw.min) {
        if (bw.min > bw.max*3/2) {
            #warning("Data is too sparse.")
            bw.max <- bw.min*1.01
        } else bw.max <- bw.max*3/2
    }
    return(list(min=bw.min, max = bw.max))
}


#' Perform local Fr\'{e}chet regression using trust package and perturbation for initial value.
#' The implementation is adapted from https://github.com/functionaldata/tFrechet
#'
#' @param xin A numeric vector of length n, the input covariates.
#' @param yin An n by d matrix of responses, each row a unit vector on the sphere.
#' @param xout A numeric vector of length k, the evaluation points.
#' @param optns A list of options, including the bandwidth and the kernel function.
#'
#' @returns A k by d matrix, where each row is the estimated regression function at \code{xout[j]}.
#' @export
#'
#' @examples
LocSpheGeoReg <- function(xin, yin, xout, optns = list()) {
    k = length(xout)
    n = length(xin)
    m = ncol(yin)
    
    bw <- optns$bw
    # ker <- kerFctn(optns$kernel)
    ker <- dnorm
    
    yout = sapply(1:k, function(j){
        mu0 = mean(ker((xout[j] - xin) / bw))
        mu1 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j]))
        mu2 = mean(ker((xout[j] - xin) / bw) * (xin - xout[j])^2)
        s = ker((xout[j] - xin) / bw) * (mu2 - mu1 * (xin - xout[j])) /
            (mu0 * mu2 - mu1^2)
        
        # initial guess
        y0 = colMeans(yin*s)
        y0 = y0 / l2norm(y0)
        if (sum(sapply(1:n, function(i) sum(yin[i,]*y0))[ker((xout[j] - xin) / bw)>0] > 1-1e-8)){
            #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) +
            #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
            # return(y0)
            y0 = y0 + rnorm(3) * 1e-3
            y0 = y0 / l2norm(y0)
        }
        
        objFctn = function(y){
            # y <- y / l2norm(y)
            if ( ! isTRUE( all.equal(l2norm(y),1) ) ) {
                return(list(value = Inf))
            }
            f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
            g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * frechet::SpheGeoGrad(yin[i,], y))) * s)
            res = sapply(1:n, function(i){
                grad_i = frechet::SpheGeoGrad(yin[i,], y)
                return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * frechet::SpheGeoHess(yin[i,], y)) * s[i])
            }, simplify = "array")
            h = 2 * apply(res, 1:2, mean)
            return(list(value=f, gradient=g, hessian=h))
        }
        res = trust::trust(objFctn, y0, 0.1, 1e5)
        # res = trust::trust(objFctn, y0, 0.1, 1)
        return(res$argument / l2norm(res$argument))
    })
    return(t(yout))
}


#' Compute the geodesic distance for points on the unit sphere.
#' The implementation is adapted from https://github.com/functionaldata/tFrechet
#'
#' @param y1,y2 Numeric vectors of the same length, representing unit vectors.
#'
#' @returns A numeric value giving the geodesic distance.
#' @export
#'
#' @examples
SpheGeoDist <- function(y1,y2) {
    if (abs(length(y1) - length(y2)) > 0) {
        stop("y1 and y2 should be of the same length.")
    }
    if ( !isTRUE( all.equal(l2norm(y1),1) ) ) {
        stop("y1 is not a unit vector.")
    }
    if ( !isTRUE( all.equal(l2norm(y2),1) ) ) {
        stop("y2 is not a unit vector.")
    }
    y1 = y1 / l2norm(y1)
    y2 = y2 / l2norm(y2)
    if (sum(y1 * y2) > 1){
        return(0)
    } else if (sum(y1*y2) < -1){
        return(pi)
    } else return(acos(sum(y1 * y2)))
}


#' Compute the \eqn{L_2} norm of a vector. 
#' The implementation is adapted from https://github.com/functionaldata/tFrechet
#'
#' @param x A numeric vector of length p.
#'
#' @returns A single numeric value representing the \eqn{L_2} norm of \code{x}.
#' @export
#'
#' @examples
l2norm <- function(x){
    as.numeric(sqrt(crossprod(x)))
}


#' The single-index Fr\'{e}chet regression for spherical data. 
#'
#' This function estimates the index coefficients and bandwidth in the single-index Fr\'{e}chet regression model. 
#' The implementation is adapted from Bhattacharjee and M\{''}uller (2023). AoS.
#' See also https://github.com/functionaldata/tFrechet
#' 
#' @param xdat An n by p matrix of input covariates.
#' @param ydat An n by d matrix of response vectors, each row a unit vector.
#' @param init A vector of length p specifying an initial vector for the index coefficients. 
#' @param kernel A character string specifying the kernel function used in local Fr\'{e}chet regression. Currently only "gauss" (Gaussian) is supported.
#' @param abstol Absolute tolerance for convergence (default = 1e-2).
#' @param maxiter Maximum number of iterations (default = 100).
#'
#' @returns A list with components:
#' \describe{
#'   \item{theta}{Estimated index coefficient vector of length p.}
#'   \item{bw}{Selected bandwidth for the local Fr\'{e}chet regression, using 5-fold cross-validation.}
#'   \item{mu}{An n by d matrix of the leave-one-out single index Fr\'{e}chet estimates.}
#'   \item{iter}{Number of iterations.}
#' }
#' 
#' @export
#'
#' @examples
fsim_est <- function(xdat, ydat, init = NULL, kernel = "gauss"){
    if (!requireNamespace("frechet")) {
        stop("Package 'frechet' is required but not installed.")
    }
    n = nrow(xdat); p = ncol(xdat)
    
    if(is.null(init)){
        theta_init = rnorm(p)
    } else {
        theta_init = init
    }
    
    # select optimal bandwidth using 5-fold CV, given theta_init
    theta_init = sign(theta_init[1]) * SpheNormalize(theta_init)
    U_init = as.vector(xdat %*% theta_init) # index values
    bw_update = bwCV_sphe(xin = U_init, yin = ydat, xout = U_init, optns = list("ker"=kernel))
    
    objFctn <- function(theta){
        theta = sign(theta)[1] * SpheNormalize(theta) # standardize to unit norm and a positive leading component
        index = xdat %*% theta
        
        # compute the empirical Fréchet variance
        Fvar = sapply(seq_len(n), function(j){
            # leave-one-out local linear Frechet estimate
            muhat_j <- LocSpheGeoReg(xin = index[-j], yin = ydat[-j,], xout = index[j], optns=list(kernel=kernel, bw = bw_update))
            SpheGeoDist(as.vector(muhat_j), ydat[j,])^2
        })
        mean(Fvar)
    }
    
    # update theta, given bandwidth
    est = optim(theta_init, fn = objFctn, method = "Nelder-Mead") # objective function could be not differentiable
    theta_update = sign(est$par[1]) * SpheNormalize(est$par)
    
    U_update = as.vector(xdat %*% theta_update) 
    mu_hat = t(sapply(seq_len(n), function(j){LocSpheGeoReg(xin = U_update[-j], yin = ydat[-j,], xout = U_update[j], 
                                                            optns=list(kernel=kernel, bw = bw_update))}))
    return(list(theta = theta_update, bw = bw_update, mu = mu_hat))
}

###########################################################
# The implementation of ESIM with ell_1 loss (SIQR)
###########################################################
#' Compute extrinsic local linear median regression for multivariate responses. Each component of responses are fitted separately.
#'
#' @param txdat Numeric vector of length n, the index values \eqn{\theta^\top X}.
#' @param tydat Numeric n by d matrix, each row a d-dimensional unit vector.
#' @param exdat Numeric scalar, the evaluation point.
#' @param bw Positive numeric, bandwidth for the local linear regression.
#'
#' @returns Numeric vector of length 2*d:
#'  - first d entries: fitted values at exdat for each component
#'  - last d entries: first derivative estimates for each component
#'
#' @export
#'
#' @examples
ll_median <- function(txdat, tydat, exdat, bw){
    d = ncol(tydat)
    z = txdat - exdat
    K = dnorm(z / bw) / bw
    
    mu = numeric(d)
    mudev = numeric(d)
    
    for(j in 1:d){
        fit = rq(tydat[,j] ~ z, tau = 0.5, weights = K)
        mu[j] = fit$coef[1]      
        mudev[j] = fit$coef[2]   
    }
    
    c(mu, mudev)
}


#' Compute a local linear quantile fit at one index value.
#'
#' This helper is used in Step 1 of the SIQR algorithm to estimate the
#' component-wise link function value and derivative at a single evaluation
#' point.
#'
#' @param x Numeric vector of length n, the leave-one-out index values.
#' @param y Numeric vector of length n, the leave-one-out response values for one component.
#' @param h Positive numeric scalar, bandwidth for the local quantile fit.
#' @param tau Numeric scalar in (0, 1), quantile level. Defaults to 0.5.
#' @param x0 Numeric scalar, the evaluation point.
#'
#' @returns A list with entries:
#'  - `x0`: the evaluation point
#'  - `fv`: fitted conditional quantile at `x0`
#'  - `dv`: local linear derivative estimate at `x0`
#'
lprq0<-function (x, y, h, tau = 0.5, x0)  #used in step 1 of the algorithm
{
    fv <- x0
    dv <- x0
    
    z <- x - x0
    wx <- dnorm(z/h)
    r <- rq(y ~ z, weights = wx, tau = tau, ci = FALSE)
    fv <- r$coef[1]
    dv <- r$coef[2]
    list(x0 = x0, fv = fv, dv = dv)
}

#' Compute component-wise local linear quantile fits.
#'
#' Each column of the multivariate response is fitted separately at the same
#' evaluation point. A scalar bandwidth is recycled across all response
#' components; a vector bandwidth uses one bandwidth per component.
#'
#' @param x0 Numeric scalar, the evaluation point.
#' @param x Numeric vector of length n, the index values.
#' @param y Numeric n by d matrix, each column a response component.
#' @param h Positive numeric scalar or length-d vector, bandwidth value(s).
#' @param tau Numeric scalar in (0, 1), quantile level. Defaults to 0.5.
#'
#' @returns Numeric vector of length d, the local linear fitted quantile for each response component.
#'
lprq_mul <- function(x0, x, y, h, tau = 0.5){
    y <- as.matrix(y)
    z <- x - x0
    if(length(h) == 1){
        h <- rep(h, ncol(y))
    }
    fv = NULL
    for (col_ in 1:ncol(y)) {
        wx <- dnorm(z/h[col_]) / h[col_]
        r <- rq(y[,col_] ~ z, weights = wx, tau = tau, ci=FALSE)
        fv = c(fv, as.numeric(r$coef[1]))
    }
    fv
}

#' Estimate the common single-index coefficient using SIQR.
#'
#' This implements the Wu (2010) iterative single-index quantile regression
#' update with a multivariate response. The same index coefficient is shared
#' across all response components, while the local link function estimates are
#' computed component-wise using a common bandwidth. The Step 2 quantile regression stacks
#' the pseudo-observations over all response components so that the objective
#' sums the component-wise quantile losses.
#'
#' @param y Numeric n by d matrix, the multivariate response.
#' @param xx Numeric n by p matrix, the predictors.
#' @param tau Numeric scalar in (0, 1), quantile level.
#' @param gamma0 Numeric vector of length p, initial index coefficient.
#' @param maxiter Positive integer, maximum number of iterations.
#' @param crit Positive numeric scalar, convergence tolerance for squared coefficient change.
#'
#' @returns A list with entries:
#'  - `theta`: normalized estimated common index coefficient
#'  - `bw`: positive numeric scalar, bandwidth from the last successful update
#'  - `flag.conv`: logical value indicating whether the iteration ended without hitting `maxiter`
#'
index.gamma<-function (y, xx, tau, gamma0, maxiter,crit) 
{
    flag.conv<-0 #flag whether maximum iteration is achieved
    
    y <- as.matrix(y)
    xx <- as.matrix(xx)
    
    gamma.new<-gamma0 #starting value
    gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2))
    
    n<-NROW(y); p<-NCOL(xx); d<-NCOL(y)
    
    iter<-1
    gamma.old<-2*gamma.new
    bw.old <- 2
    update.ok <- TRUE
    
    while((iter < maxiter) & (sum((gamma.new-gamma.old)^2)>crit))
    {
        gamma.old<-gamma.new
        iter<-iter+1
        update.ok <- tryCatch({
        x = xx %*% t(t(gamma.old))
        bw.candidate <- sapply(seq_len(d), function(col_){
            tryCatch(
                KernSmooth::dpill(as.vector(x), y[,col_]) *
                    (tau*(1-tau)/(dnorm(qnorm(tau)))^2)^.2,
                error = function(e) NA_real_
            )
        })
        bw.fallback <- stats::sd(as.vector(x)) * n^(-1/5)
        bw.new <- mean(bw.candidate, na.rm = TRUE)
        if(is.na(bw.new)){
            bw.new <- bw.fallback
        }
        
        y.stack <- numeric(n^2*d) # pseudo-observations used for step 2
        x.stack <- matrix(0, nrow = n^2*d, ncol = p) # pseudo-observations used for step 2
        w.stack <- numeric(n^2*d)
        pair.i <- rep(seq_len(n), each = n)
        pair.j <- rep(seq_len(n), times = n)
        
        for (col_ in 1:d) {
            ####################################
            #  step 1: compute a_j,b_j; j=1:n  #
            ####################################
            ab.fit <- sapply(seq_len(n), function(j){
                fit <- lprq0(x[-j], y[-j,col_], bw.new, tau, x[j])
                c(fit$fv, fit$dv)
            })
            a <- ab.fit[1,]
            b <- ab.fit[2,]
            
            ############################################################# 
            # Build Wu's pseudo-observations y*_ij and x*_ij for step 2 #
            #############################################################
            ynew <- y[pair.i, col_] - a[pair.j]
            xnew <- (xx[pair.i,,drop = FALSE] - xx[pair.j,,drop = FALSE]) * b[pair.j]
            
            # Weights are evaluated at the previous gamma estimate and current bandwidth.
            index.diff <- as.vector((xx[pair.i,,drop = FALSE] - xx[pair.j,,drop = FALSE]) %*% gamma.old)
            wts <- dnorm(index.diff/bw.new)/bw.new
            idx <- ((col_-1)*n^2+1):(col_*n^2)
            y.stack[idx] <- ynew
            x.stack[idx,] <- xnew
            w.stack[idx] <- wts
        }
        # filter for valid pseudo-observations
        keep <- is.finite(y.stack) & is.finite(w.stack) & w.stack > 0 &
            apply(is.finite(x.stack), 1, all)
        if(sum(keep) <= p){
            warning("Too few valid pseudo-observations for gamma update.")
        }
        fit <- rq(y.stack[keep]~0+x.stack[keep,], weights = w.stack[keep],
                  tau = tau, method = "fn") ; #fn for very large problems and `0` to exclude intercept
        gamma.new <- as.numeric(fit$coef)
        bw.old <- bw.new
        TRUE
        }, error = function(e){
            print(paste("Error occurred in iteration", iter))
            FALSE
        })
        if(!update.ok){break}
        
        gamma.new<-sign(gamma.new[1])*gamma.new/sqrt(sum(gamma.new^2))   #normalize
        
    }
    
    iter<-iter
    flag.conv <- update.ok && iter < maxiter # = 1 if converge; =0 if not converge
    gamma<-gamma.new
    list(theta=gamma, bw = bw.old, flag.conv=flag.conv)
}


###########################################################
# generate random sample from SvMF 
# The implementations are provided in the supplementary material of Scealy and Wood (2019), JASA. 
###########################################################
simKent=function(kappa,beta,mu,K,n){
    p = length(mu)
    skappa=matrix(kappa,n,1)
    
    sbeta=matrix(0,n,sum(p,-1))
    for (j in 2:sum(p,-1))
    {
        sbeta[,j]=beta[j-1]
    }
    
    #simulate sample
    
    
    cum=1
    rej=0
    rej2=0
    
    zs=matrix(0,n,p)
    zc=matrix(0,1,p)
    sig=matrix(0,n,p)
    sbetas=matrix(0,n,1)
    
    
    for (i in 2:sum(p,-1))
    {
        sbetas=sbetas+sbeta[,i]
        for (j in 1:n)
        {
            if (sbeta[j,i] > 0 ) {sig[j,i]=sqrt(skappa[j]-2*sbeta[j,i])}
            else {sig[j,i]=sqrt(skappa[j])}
            if (sbetas[j] < 0 ) {sig[j,p]=sqrt(skappa[j]+2*sbetas[j])}
            else {sig[j,p]=sqrt(skappa[j])}
            
        }
    }
    
    
    
    
    while (cum <= n)
        
    {
        
        
        for (i in 2:p)
        {
            zc[1,i]=rexp(1,rate=sig[cum,i])
        }
        
        vc=sum(zc^2)/4
        
        if (vc < 1)
        {
            
            
            r=runif(1, min=0, max=1)
            
            bz=0
            for (i in 2:sum(p,-1))
            {
                bz=bz + sbeta[cum,i]*zc[1,i]^2
            }
            ez=0
            for (i in 2:p)
            {
                ez=ez + sig[cum,i]*zc[1,i]
            }
            
            paccept=exp(((p-3)/2)*log(1-vc)-2*vc*skappa[cum]+(1-vc)*(bz-sbetas[cum]*zc[1,p]^2)+ ez-((p-1)/2))
            
            
            if (r < paccept) 
            {
                r2=runif(p, min=0, max=1)
                for (i in 2:p)
                {
                    if (r2[i]< 0.5) {zc[1,i]=-1*zc[1,i]}
                    zs[cum,i]=zc[1,i]
                }
                
                cum=cum+1
            }
            else
            {
                rej2=rej2+1	
            }
            
            
        }
        else
        {
            rej=rej+1	
        }
        
        
    }
    
    reject<<-rej

    H=diag(1,p)
    H[,1]=t(t(mu))
    H[1,]=t(mu)
    mu_L=t(t(mu[2:p]))
    H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
    
    Gamma=H%*%K
    
    vs=matrix(0,n,1)
    ys=matrix(0,n,p)
    y=matrix(0,n,p)
    fold=0
    for (i in 1:n)
    {
        vs[i]=sum(zs[i,]^2)/4
        ys[i,1]=1-2*vs[i]
        ys[i,2:p]=((1-vs[i])^(0.5))*zs[i,2:p]
        y[i,]=ys[i,]%*%t(Gamma)
        fold_count=0
        ##count folding
        for (j in 1:p)
        {
            if (y[i,j] < 0) {fold_count=1}
        }
        fold=fold+fold_count
    }
    fold=fold/n
    
    return(list(y=y,fold=fold,ystd=ys))
}

simProject=function(kappa,V,mu,a1,n){
    p = length(mu)
    kappav=kappa
    av=a1
    
    betav=matrix(0,p-2,1)
    betav[1]=0*kappav
    muv=matrix(0,p,1)
    muv[1]=1
    K=diag(p)
    #simulated sample from von mises fisher distribution:
    sims=simKent(kappav,betav,muv,K,n)
    yv=sims$y
    
    a=eigen(V)$values
    a=rbind(av^2,t(t(a)))
    a=sqrt(a)
    
    yp=matrix(0,n,p)
    for (j in 1:p)
    {
        yp[,j]=yv[,j]*a[j]
        
    }
    
    sumsq=0
    for (j in 1:p)
    {
        sumsq=sumsq+yp[,j]^2
    }
    yp=yp/sqrt(sumsq)
    
    
    H=diag(1,p)
    H[,1]=t(t(mu))
    H[1,]=t(mu)
    mu_L=t(t(mu[2:p]))
    H[2:p,2:p]=(1/(1+H[1,1]))*mu_L%*%t(mu_L)-diag(1,sum(p,-1))
    
    K=diag(1,p)
    K[2:p,2:p]=eigen(V)$vectors
    
    Gamma=H%*%K
    
    ypn=matrix(0,n,p)
    fold=0
    for (i in 1:n)
    {
        ypn[i,]=yp[i,]%*%t(Gamma)
        
        
    }
    
    
    
    #response
    y=ypn
    
    
    
    return(y)
}




# Lin_est <- function(xdat, ydat){
#     # LOOCV
#     dx = ncol(xdat)
#     dy = ncol(ydat)
#     n = nrow(ydat)
#     h0 = rep(0.2, dx)
#     
#     flag = FALSE
#     while (flag!=TRUE) {
#         tryCatch({
#             flag = TRUE
#             optim.bw = optim(par=h0, fn=function(bw){
#                 fit = sapply(1:dy, function(col){
#                     np::npksum(txdat=xdat, tydat = ydat[,col], bws = bw, leave.one.out=T)$ksum/
#                         np::npksum(txdat=xdat, bws = bw, leave.one.out=T)$ksum})
#                 if(any(is.na(fit))){return(9999)}
#                 fit = t(apply(fit, 1, function(x){x/sqrt(sum(x^2))}))
#                 h0 = runif(dx, 0.5, 1.5) * h0
#                 sum(sapply(1:n, function(i){acos(sum(fit[i,]*ydat[i,]))}))
#             },
#             method = "L-BFGS-B", lower = rep(0.01,d), upper = rep(1.5, d)
#             )}, 
#             error=function(e){
#                 flag = FALSE
#                 warning("Optim error!")
#             })
#     }
#     out = sapply(1:dy, function(col){
#         np::npksum(txdat=xdat, tydat = ydat[,col], bws = optim.bw$par, leave.one.out=T)$ksum/
#             np::npksum(txdat=xdat, bws = optim.bw$par, leave.one.out=T)$ksum})
#     out = t(apply(out, 1, function(x){x/sqrt(sum(x^2))}))
#     return(list("bw"=optim.bw$par, "fit"=out))
# }
# 
# Lin_fit <- function(xdat, ydat, xtest, h, K=dnorm){
#     # input: 
#     #       xdat, xtest, ydat: n by dx/dy
#     #       h: dx by 1 bandwidth
#     #       output: n by dy matrix
#     if (ncol(xdat)!=length(h)) {
#         stop("dimensions of xdat should be the same as length of h")
#     }
#     dy = ncol(ydat)
#     fit = matrix(0, nrow = nrow(xtest), ncol = dy)
#     for (i in 1:nrow(xtest)) {
#         Kx = t(apply(xdat, 1, function(x){K((x-xtest[i,])/h)/h}))
#         Ki = apply(Kx, 1, prod)
#         sumK = sum(Ki)
#         Ki = Ki/sumK
#         Ky = t(apply(cbind(ydat,Ki), 1, function(x){x[dy+1]*x[1:dy]}))
#         fit[i,] = colSums(Ky)
#     }
#     ## normalization
#     fit.norm = t(apply(fit, 1,function(x){x/sqrt(sum(x^2))}))
#     return(fit.norm)
# }

# cost_ESIM <- function(param, xdat, ydat){
#  extrinsic single-index model using the local constant.
#     # input:
#     #   param: d+1 * 1, bandwidth and theta
#     #   xdat, ydat: n*p or n*d
#     # output: scalar, ojective value
#     
#     require(np)
#     MaxPanelty = .Machine$double.xmax
#     K = dnorm # kernel function
#     dx = ncol(xdat); dy = ncol(ydat)
#     n = nrow(xdat)
#     
#     if(dx != (length(param)-1)){
#         stop("xdat has different dimensions with param.")
#     }
#     
#     h = param[1] # bandwidth
#     beta = param[2:(dx+1)] # parametric coefficients
#     beta = SpheNormalize(beta)
#     
#     if(param[1]<=1e-4){return(MaxPanelty)}
#     
#     if(sum(is.na(param))==0){
#         # return an infinite penalty for negative h
#         fit = matrix(NA, nrow = nrow(xdat), ncol = dy)
#         index = as.matrix(xdat) %*% beta
#         # matrix of kernel
#         fit = sapply(1:dy, function(col){
#             np::npksum(txdat=index, tydat = ydat[,col], bws = h, leave.one.out=T)$ksum/
#                 np::npksum(txdat=index, bws = h, leave.one.out=T)$ksum
#         })
#         if(any(is.na(fit))){
#             val = MaxPanelty
#         } else {
#             val = sum(apply(ydat-fit, 1, function(x){sum(x^2)}))
#         }
#         
#     } else {val = MaxPanelty}
#     return(val)
# }

# cost_ESIMsep <- function(param, xdat, ydat){
#     # fit separately
#     MaxPanelty = .Machine$double.xmax
#     ## verify dimension of xdat, ydat
#     K=dnorm
#     ydat <- t(t(ydat))
#     xdat <- t(t(xdat))
#     dx = ncol(xdat)
#     dy = ncol(ydat)
#     n = nrow(xdat)
#     
#     h = param[1:dy]
#     beta = param[(dy+1):length(param)]
#     beta = sign(beta[1])*beta/sqrt(sum(beta^2)) # standardize
#     
#     if (dx!=(length(beta))){
#         stop("xdat has different dimensions with param")
#     }
#     if (length(h)!= dy) {
#         stop("ydat has different dimensions with h")
#     }
#     
#     index = xdat %*% t(t(beta))
#     if ((min(h)>1e-2) && (max(h)<=1.5)) {
#         fit = sapply(1:dy, function(col){
#             np::npksum(txdat=index, tydat = ydat[,col], bws = h[col], leave.one.out=T)$ksum/
#                 np::npksum(txdat=index, bws = h[col], leave.one.out=T)$ksum
#         })
#         if(any(is.na(fit))){
#             val = MaxPanelty
#         } else {
#             val = sum(apply(ydat-fit, 1, function(x){sum(x^2)}))
#         }
#     } else {
#         val = MaxPanelty
#     }
#     
#     return(val)
# }
# 
# 
# # input:
# #   xdat: n*1, index value for training
# #   ydat: n*d, response for training
# #   exdat: scalar
# #   bw: scalar, bandwidth
# # output: d*1, fitted value at exdat
# ESIM_fit <- function(xdat, ydat, exdat, bw){
#     d = ncol(ydat)
#     weight = matrix(rep(dnorm((xdat - exdat)/bw)/bw, d), ncol = d)
#     fit = colSums( ydat * weight )/ sum(weight[,1])
#     ## NA value?
#     return(fit)
# }
# 
# ESIMsep_fit <- function(xdat, ydat, xtest, h, K=dnorm, pred=TRUE){
#     # input: 
#     #        xdat, xtest: n by 1 vector, index values
#     #        ydat: n by p matrix
#     #        h: bandwidth for each component of y
#     # output: ntest by p matrix
#     if (ncol(ydat)!= length(h)) {
#         stop("dimension of ydat should be the same as length of h when fitting single index model separately.")
#     }
#     fit = matrix(0, nrow = nrow(xtest), ncol = ncol(ydat))
#     for (i in 1:nrow(xtest)) {
#         if (pred==FALSE) {
#             df = xdat[-i] - xtest[i]
#             Y_ = ydat[-i,]
#         } else {
#             df = xdat - xtest[i]
#             Y_ = ydat
#         }
#         # d = df/max(abs(df))
#         d = df
#         Kx = sapply(h, function(hi){K(d/hi)/hi})
#         Ki = apply(Kx,2,function(x){x/sum(x)})
#         Ky = t(sapply(1:nrow(Y_), function(x){Y_[x,]*Ki[x,]}))
#         fit[i,] = colSums(Ky)
#     }
#     ## normalization
#     fit.norm = t(apply(fit, 1,function(x){x/sqrt(sum(x^2))}))
#     return(fit.norm)
# }

# esl_est <- function(param, xdat, ydat, lambda, abstol=1e-4, maxiter=100){
#     dist_ = 1e5
#     iter = 1
#     d = ncol(ydat)
#     p = ncol(xdat)
#     n = nrow(xdat)
#     weight = matrix(1, nrow = n-1, ncol = n)
#     
#     while (dist_ >= abstol && iter <= maxiter) {
#         # given update theta and bw
#         optim.para = optim(param, cost_index_bw, xdat=xdat, ydat=ydat, lambda=lambda, weight=weight)
#         theta.update = optim.para$par[-1]
#         theta.update = sign(theta.update[1])*theta.update/sqrt(sum(theta.update^2))# normalize
#         bw.update = optim.para$par[1]
#         dist_ = sum((theta.update - param[-1])^2)
#         param = c(bw.update, theta.update)
#         
#         # given theta, bw updates mu and mudev
#         index_update = xdat %*% t(t(theta.update)) # n * 1
#         llmu_mudev = t(sapply(1:nrow(ydat), function(r_){
#             llmu_ESL(txdat = index_update[-r_],
#                      tydat = ydat[-r_,],
#                      exdat = index_update[r_],
#                      bw.update, weight[,r_])
#         })) # n * 2d
#         mu.update = llmu_mudev[,1:d]
#         mudev.update = llmu_mudev[,(d+1):(2*d)]
#         
#         # given theta, bw and mu, update weight
#         weight = sapply(1:nrow(xdat), function(r_){
#             a = matrix(rep((index_update[-r_] - index_update[r_]), d), ncol = d)
#             e = ydat[-r_,] - mu.update[-r_,] - mudev.update[-r_,]*a
#             exp(- rowSums(e^2)/lambda)
#         })
#         iter = iter + 1
#     }
#     
#     list(theta = theta.update, bw = bw.update, mu = mu.update, mudev = mudev.update)
# }

# ISIM_extrinsic <- function(beta0, h0, xdat, ydat, maxiter=30, crit = 1e-4){
#     ydat = t(t(ydat))
#     dx = ncol(xdat)
#     dy = ncol(ydat)
#     n = nrow(xdat)
#     
#     beta0 = sign(beta0[1]) * beta0/ sqrt(sum(beta0^2))
#     I = xdat %*% beta0
#     
#     iter = 0
#     beta.new = beta0
#     beta.old = 2 * beta.new
#     h.optim = h0
#     while ((iter<maxiter)&(sum((beta.old-beta.new)^2)>crit)) {
#         iter = iter + 1
#         beta.old = beta.new
#         ISIM = optim(c(h.optim, beta.old), fn = cost_ISIM_extrinsic, xdat=xdat, ydat=ydat)$par
#         h.optim = ISIM[1]
#         beta.new = ISIM[-1]
#         beta.new = sign(beta.new[1])*beta.new/sqrt(sum(beta.new^2))
#     }
#     
#     I = xdat %*% beta.new
#     fit = t(sapply(1:n, function(i){ISIM_fit_extrinsic(I[i], I[-i], ydat[-i,], h.optim)}))
#     return(list("beta"=beta.new, "bw"=h.optim, "fit"=fit))
# }
# 
# ISIM_fit_extrinsic <- function(xout, xin, yin, bw){
#     mu0 = mean(dnorm((xout - xin) / bw))
#     mu1 = mean(dnorm((xout - xin) / bw) * (xin - xout))
#     mu2 = mean(dnorm((xout - xin) / bw) * (xin - xout)^2)
#     s = dnorm((xout - xin) / bw) * (mu2 - mu1 * (xin - xout)) / (mu0 * mu2 - mu1^2)
#     
#     # y0 = colMeans(yin * s) # initial guess
#     # y0 = y0/sqrt(sum(y0^2))
#     y = colMeans(yin *s)/mean(s)
#     fit = y/sqrt(sum(y^2))
#     return(fit)
# }
# 
# cost_ISIM_extrinsic <- function(param, xdat, ydat){
#     ## verify dimension of xdat, ydat
#     # minimize d(y_i,w) = |y_i - w|^2
#     K = dnorm
#     ydat = t(t(ydat))
#     dx = ncol(xdat)
#     dy = ncol(ydat)
#     n = nrow(xdat)
#     
#     h = param[1]
#     beta = param[-1]
#     beta = sign(beta[1])*beta/sqrt(sum(beta^2)) # standardize
#     
#     if (dx!=(length(beta))){
#         stop("xdat has different dimensions with param")
#     }
#     
#     I = as.vector(xdat %*% t(t(beta)))
#     
#     if((h <= 1.5) & (h>1e-2) ){
#         reg = t(sapply(1:n, function(i){ISIM_fit_extrinsic(I[i], I[-i], ydat[-i,], bw=h)}))
#         if(any(is.nan(reg))){
#             return(1e9)
#         } else {
#             res = sapply(1:n, function(i){SpheGeoDist(ydat[i,], reg[i,])})
#             return(sum(res))
#         }
#         
#     }
#     return(1e9)
# }

# # input:
# #   param0: (p+1)*1, init value for bw and theta
# #   xdat, ydat: n*p n*d, 
# #   bw: bandwidth 
# #   weight: n-1 * n, 
# # output: scalar, minimized optimum
# cost_siqr <- function(param, xdat, ydat, bw){
#     theta = SpheNormalize(param)
#     theta = sign(theta[1])*theta/sqrt(sum(theta^2))
#     
#     index_ = xdat %*% t(t(theta))
#     # leave one out local linear of mu and first derivative
#     llmu = sapply(1:nrow(ydat), function(r_){
#         llmu_median(txdat = index_[-r_], tydat = ydat[-r_,], exdat = index_[r_], bw)[1:3]
#     })
#     llmu = t(llmu)
#     # objective value
#     if (any(is.nan(llmu))) {return(.Machine$double.xmax)}
#     mean(rowSums(abs(ydat - llmu)))
# }

# # input:
# #   param: p*1, init value for bw and theta
# #   xdat, ydat: n*p n*d, 
# # output: list of theta, bw, mu and first derivative.
# siqr_est <- function(param, xdat, ydat){
#     dist_ = 1e5
#     iter = 1
#     d = ncol(ydat)
#     p = ncol(xdat)
#     n = nrow(xdat)
#     param = sign(param[1]) * param / sqrt(sum(param^2))
#     # estimate bw
#     hm = mean(apply(ydat, 2, dpill, x = xdat %*% param))
#     bw.update = hm * (0.5*(1-0.5)/(dnorm(qnorm(0.5)))^2)^.2
#     # given bw update theta
#     theta.update = optim(param, cost_index_siqr, xdat=xdat, ydat=ydat, bw=bw.update)$par
#     theta.update = sign(theta.update[1]) * theta.update / sqrt(sum(theta.update^2))# normalize
#     
#     # 
#     mu.qr = sapply(1:nrow(ydat), function(r_){
#         llmu_median(txdat = (xdat %*% t(t(theta.update)))[-r_], 
#                     tydat = ydat[-r_,], exdat = (xdat %*% t(t(theta.update)))[r_], bw.update)
#     })
#     
#     list(theta = theta.update, bw = bw.update, mu = t(mu.qr[1:3,]), mudev = t(mu.qr[4:6,]))
# }

# # select number of bins
# binsCV_sphe <- function(xin, yin, optns){
# 
#     objFctn <- function(M) {
#         bins <- seq(min(xin), max(xin), length.out = M+1) # split the support into M bins
#         bin_id <- findInterval(xin, bins, rightmost.closed = T) # bin lables for each observation
#         # ensure each bin includes at least one observation
#         min_bin <- min(sapply(seq_len(M), function(b){length(xin[bin_id==b])})) # minimum number of observation in bins
#         if(min_bin >= 1){
#             bin_var = sapply(seq_len(M), function(b){
#                 # find representative data for each bin
#                 xbin <- xin[bin_id == b]
#                 ybin <- yin[bin_id == b,]
#                 ybin <- ybin[order(xbin),]
#                 xbin <- sort(xbin)
#                 mid_idx <- ceiling(length(xbin)/2)
#                 xbin_rep <- xbin[mid_idx]; ybin_rep = ybin[mid_idx,]
#                 # fit local Frechet regression at the representative data
#                 res <- LocSpheGeoReg(xin = xin[bin_id != b], yin = yin[bin_id != b,], xout = xbin_rep, optns = optns)
#                 # obtain empirical Frechet variance
#                 SpheGeoDist(as.vector(res), ybin_rep)^2
#             })
#             mean(bin_var)
#         } else {
#             return(99999)
#         }
#     }
#     
#     result <- optimize(f = objFctn, interval = c(2, n))
#     result$minimum
# }

#' #' The single-index Fr\'{e}chet regression for spherical data. 
#' #' @param xdat 
#' #' @param ydat 
#' #' @param init 
#' #' @param n_theta A positive integer specifying the number of random index coefficients to be generated.
#' #' @param kernel 
#' #'
#' #' @returns 
#' #' 
#' #' @export
#' #'
#' #' @examples
#' fsim_est <- function(xdat, ydat, init = NULL, n_theta = 500, kernel = "gauss"){
#'     if (!requireNamespace("mvtnorm")) {
#'         stop("Package 'mvtnorm' is required but not installed.")
#'     }
#'     if (!requireNamespace("foreach")) {
#'         stop("Package 'foreach' is required but not installed.")
#'     }
#'     if (!requireNamespace("doParallel")) {
#'         stop("Package 'doParallel' is required but not installed.")
#'     }
#'     n <- nrow(xdat); p <- ncol(xdat)
#'     
#'     # generating n_theta random directions as initial values of theta
#'     theta_init = mvtnorm::rmvnorm(n_theta, mean = rep(0, p), sigma = diag(p)) # generate p-dimensional standard Gaussian random vectors
#'     theta_init = t(apply(theta_init, 1, function(r){
#'         ifelse(r[1] > 0, 1, -1) * r / sqrt(sum(r^2)) # change the first element to be positive and normalized 
#'     }))
#'     if(!is.null(init)){theta_init = rbind(init, theta_init)}
#'     
#'     # initialize parallel computing 
#'     num_cores <- parallel::detectCores() - 1
#'     cl <- parallel::makeCluster(num_cores)
#'     parallel::clusterExport(cl, c("bwCV_sphe", "LocSpheGeoReg", "SpheGeoDist", "SetBwRange", "l2norm"))
#'     doParallel::registerDoParallel(cl)
#'     
#'     results <- foreach::foreach(i = 1:nrow(theta_init), .packages = c("frechet"), .combine = rbind) %dopar% {
#'         theta <- theta_init[i,]
#'         
#'         # select optimal bandwidth for each theta_init using 5-fold CV
#'         index <- xdat %*% theta
#'         bw <- bwCV_sphe(xin = index, yin = ydat, xout = index, optns = list("ker"=kernel))
#'         # compute the empirical Fréchet variance
#'         Fdeviation <- sapply(seq_len(n), function(j){
#'             # leave-one-out local linear estimate
#'             muhat_j <- LocSpheGeoReg(xin = index[-j], yin = ydat[-j,], xout = index[j], optns=list(kernel=kernel, bw = bw))
#'             SpheGeoDist(as.vector(muhat_j), ydat[j,])^2
#'         })
#'         c(Fvar = mean(Fdeviation), bw = bw)
#'     }
#'     
#'     parallel::stopCluster(cl)
#'     
#'     min_idx <- which.min(results[,"Fvar"])
#'     thetahat <- theta_init[min_idx,]
#'     bwhat <- results[min_idx, "bw"]
#'     
#'     return(list(theta = thetahat, bw = bwhat, results = results))
#' }
#' # fsim_est <- function(xdat, ydat, init = NULL, kernel = "gauss", abstol = 1e-2, maxiter = 100){
#     if (!requireNamespace("frechet")) {
#         stop("Package 'frechet' is required but not installed.")
#     }
#     n = nrow(xdat); p = ncol(xdat)
#     
#     if(is.null(init)){
#         theta_init = rnorm(p)
#     } else {
#         theta_init = init
#     }
#     
#     # set up before iteration
#     theta_init = sign(theta_init[1]) * SpheNormalize(theta_init)
#     theta_update = theta_init
#     iter = 0 # counter for iteration
#     dist = 9999
#     
#     objFctn <- function(theta){
#         theta = sign(theta)[1] * SpheNormalize(theta) # standardize to unit norm and a positive leading component
#         index = xdat %*% theta
#         
#         # compute the empirical Fréchet variance
#         Fvar = sapply(seq_len(n), function(j){
#             # leave-one-out local linear Frechet estimate
#             muhat_j <- LocSpheGeoReg(xin = index[-j], yin = ydat[-j,], xout = index[j], optns=list(kernel=kernel, bw = bw_update))
#             SpheGeoDist(as.vector(muhat_j), ydat[j,])^2
#         })
#         mean(Fvar)
#     }
#     
#     while ((dist > abstol) && (iter < maxiter)) {
#         # select optimal bandwidth using 5-fold CV, given theta (only update bandwidth when the change in theta is large)
#         if (dist > .5){
#             U_update = as.vector(xdat %*% theta_update) # index values
#             bw_update = bwCV_sphe(xin = U_update, yin = ydat, xout = U_update, optns = list("ker"=kernel))
#         }
#                 
#         # update theta, given bandwidth
#         est = optim(theta_update, fn = objFctn, method = "Nelder-Mead") # objective function could be not differentiable
#         theta_update = sign(est$par[1]) * SpheNormalize(est$par)
#         
#         dist = SpheGeoDist(theta_init, theta_update) # change in the estimate
#         theta_init = theta_update # save theta in the current iteration
#         iter = iter + 1
#     }
#         
#     if (iter == maxiter){
#         message("FSIM fail to converge.")
#     }
#     
#     U_update = as.vector(xdat %*% theta_update) 
#     mu_hat = t(sapply(seq_len(n), function(j){LocSpheGeoReg(xin = U_update[-j], yin = ydat[-j,], xout = U_update[j], 
#                                                           optns=list(kernel=kernel, bw = bw_update))}))
#     return(list(theta = theta_update, bw = bw_update, mu = mu_hat, iter = iter))
# }

#' Estimation of extrinsic single-index model using the exponential squared loss
#' 
#' The estimation iterates the estimation of the tuning parameter \eqn{\lambda} 
#' and the estimation of the parameter coefficients \eqn{\theta} and the bandwidth.
#'  
#' @param theta_init A numeric vector of length p, the initial value of \eqn{\theta}.
#' @param bw_init A positive scalar, the initial value of the bandwidth.
#' @param xdat An n by p matrix, the p-dimensional covariates.
#' @param ydat An n by d matrix of response vectors, each row a d-dimensional unit vector.
#' @param delta A scalar within (0,1) controling the robustness. See Section 5 of the paper for details.
#' @param abstol A small positive scalar, the stopping criterion for the iteration based on parameter changes. Default is 1e-3.
#' @param maxiter An integer, the maximum number of iterations allowed. Default is 100.
#'
#' @returns A list containing the following elements:
#' \describe{
#'   \item{theta}{The estimated p-dimensional parameter vector \eqn{\theta}.}
#'   \item{bw}{The estimated bandwidth.}
#'   \item{mu}{An n by d matrix of the estimated local linear regression values.}
#'   \item{mudev}{An n by d matrix of the estimated first derivatives of the local linear regression.}
#'   \item{lambda}{The estimated tuning parameter \eqn{\lambda}.}
#'   \item{psi}{An (n-1) by n matrix of weights used in the local linear regression.}
#' }
#' 
#' @export
#'
#' @examples
# esl_est <- function(theta_init, bw_init, xdat, ydat, delta, abstol=1e-3, maxiter=100){
#     dist_ = 1e5
#     iter = 0
#     n = nrow(ydat); d = ncol(ydat)
#     
#     theta_update = sign(theta_init[1]) * SpheNormalize(theta_init) # normalize
#     bw_update = bw_init
#     
#     # use the ESIM with least squares to obtain the muhat
#     U = xdat %*% theta_update
#     llmu_mudev = t(sapply(1:n, function(r_){
#         ll_LS(txdat = U[-r_], tydat = ydat[-r_,], exdat = U[r_], bw = bw_update)})) # n * 2d
#     mu_update = llmu_mudev[,1:d]
#     
#     while (dist_ >= abstol && iter < maxiter) {
#         # given mu_update, update the tuning parameter lambda
#         lambda_update = optimise(cost_lambda_median, lower = 1e-5, upper = 3,
#                                  delta = delta, ydat = ydat, mu_fit = mu_update, maximum = FALSE)$minimum
#         
#         # given lambda, update theta and bw
#         result = esl_index_bw(theta_init = theta_update, bw_init = bw_update,
#                               xdat = xdat, ydat = ydat, lambda = lambda_update)
#         
#         dist_ = SpheGeoDist(theta_update, result$theta)
#         iter = iter + 1 # update counter
#         theta_update = result$theta
#         bw_update = result$bw
#         mu_update = result$mu
#         any(is.na(mu_update))
#     }
#     
#     list(theta = result$theta, bw = result$bw, mu = result$mu, mudev = result$mudev,
#          lambda = lambda_update, psi = result$psi, iter = iter)
# }
