## simulation for robust estimaiton of the ESIM with contaminated samples.
## last change: 19.Nov.25

rm(list = ls())
# source("ESLFns.R")
# library(Directional)
# library(mvtnorm)
# library(MASS)
library(foreach)
library(doSNOW)
library(rgl)

B = 500
n = 200 # n = 200 or 500 or 2000
type = "vmf"

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
clusterEvalQ(cl, { source("ESLFns.R") })
registerDoSNOW(cl)

theta_LS_simu <- theta_ESL_simu <- list()
Sigma_LS_simu <- Sigma_ESL_simu <- list()
for (epsilon in c(0, 10, 20)) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    progress <- function(b) setTxtProgressBar(pb, b)
    
    results <- foreach(b = 1:B, .packages = c("mvtnorm","Directional", "MASS"), 
                       .options.snow = list(progress=progress)) %dopar% {
        # generate random sample #####
        p = 3; d = 3; k = 20
        
        set.seed(k+b+epsilon)
        
        theta0 = c(1, -1, 1)
        theta0 = SpheNormalize(theta0)
        
        X = sapply(1:p, function(i) runif(n, -2, 2))
        U = 1/(1 + exp(-X %*% t(t(theta0))))
        mu = cbind(sqrt(1-U^2) * cos(pi*U),
                   sqrt(1-U^2) * sin(pi*U),
                   U)
    
        location = sample(1:n, n*epsilon/100)
        Y = t(apply(mu, 1, function(x){rvmf(1, x, k)})) 
        Y[location,] = -mu[location,] # antipodal contamination
        
        lambda = 3 / k # fixed lambda when fitting ESL
        
        LS = ls_est(c(1, rep(1,p)), xdat = X, ydat = Y)
        ESL = esl_index_bw(theta_init = LS$theta, bw_init = LS$bw,
                           xdat = X, ydat = Y, lambda = lambda)
        
        # obtain asymptotic covariance matrix of beta
        J_LS = rbind(- LS$theta[-1] / LS$theta[1], diag(p-1))  # Jacobian submatrix J
        J_ESL = rbind(- ESL$theta[-1] / ESL$theta[1], diag(p-1))
        
        U_LS = X %*% LS$theta # index value 
        U_ESL = X %*% ESL$theta
        
        Wj_LS = array(NA, dim = c(n, p-1, p-1))
        Wj_ESL = array(NA, dim = c(n, p-1, p-1))
        Mj_LS = array(NA, dim = c(n, p-1, p-1))
        Mj_ESL = array(NA, dim = c(n, p-1, p-1))
        
        for (j in 1:n) {
            mudev_LS = LS$mudev[j,] # local linear estimator of the first derivative
            mudev_ESL = ESL$mudev[j,]
            
            condExp_LS = ll_LS(txdat = U_LS[-j], tydat = X[-j,], exdat = U_LS[j], bw = 0.3)[1:p] # conditional expectation of X, given the index E[X|U_LS]
            condExp_ESL = ll_LS(txdat = U_ESL[-j], tydat = X[-j,], exdat = U_ESL[j], bw = 0.3)[1:p]
            
            Sigma_X_LS = t(J_LS) %*% (X[j,] - condExp_LS) %*% t(X[j,] - condExp_LS) %*% J_LS
            Sigma_X_ESL = t(J_ESL) %*% (X[j,] - condExp_ESL) %*% t(X[j,] - condExp_ESL) %*% J_ESL
            
            G_LS = (Y[j,] - LS$mu[j,]) %*% t(Y[j,] - LS$mu[j,])
            G_ESL = exp(- 2 * sum((Y[j,] - ESL$mu[j,])^2) / lambda) * (Y[j,] - ESL$mu[j,]) %*% t(Y[j,] - ESL$mu[j,])
            F_ESL = exp(- sum((Y[j,] - ESL$mu[j,])^2) / lambda) * (diag(d) - (2 * (Y[j,] - ESL$mu[j,]) %*% t(Y[j,] - ESL$mu[j,])) / lambda)
            
            Wj_LS[j,,] = as.numeric(t(mudev_LS) %*% mudev_LS) *  Sigma_X_LS 
            Wj_ESL[j,,] = as.numeric(t(mudev_ESL) %*% F_ESL %*% mudev_ESL) * Sigma_X_ESL
                
            Mj_LS[j,,] = as.numeric(t(mudev_LS) %*% G_LS %*% mudev_LS) * Sigma_X_LS
            Mj_ESL[j,,] = as.numeric(t(mudev_ESL) %*% G_ESL %*% mudev_ESL) * Sigma_X_ESL 
        }
        W_LS = apply(Wj_LS, 2:3, mean)
        M_LS = apply(Mj_LS, 2:3, mean)
        W_ESL = apply(Wj_ESL, 2:3, mean)
        M_ESL = apply(Mj_ESL, 2:3, mean)
        
        Sigma_LS = solve(W_LS) %*% M_LS %*% solve(W_LS) / n
        Sigma_ESL = solve(W_ESL) %*% M_ESL %*% solve(W_ESL) / n
        
        list(theta_LS = LS$theta, Sigma_LS = Sigma_LS,
             theta_ESL = ESL$theta, Sigma_ESL = Sigma_ESL)}
    
    # save results #####
    theta_ESL_simu[[paste("eps", epsilon, sep="")]] <- t(sapply(results, function(x) x$theta_ESL))
    theta_LS_simu[[paste("eps", epsilon, sep="")]] <- t(sapply(results, function(x) x$theta_LS))

    Sigma_ESL_simu[[paste("eps", epsilon, sep="")]] <- sapply(results, function(x) x$Sigma_ESL)
    Sigma_LS_simu[[paste("eps", epsilon, sep="")]] <- sapply(results, function(x) x$Sigma_LS)
}
save(theta_ESL_simu, theta_LS_simu, Sigma_ESL_simu, Sigma_LS_simu, 
     file = paste("./data/simu_cov_", type,"_n",n,".RData", sep = ""))
close(pb)
stopCluster(cl)

# table ####
rm(list = ls())
n = 2000 # n=200, 500, 2000
load(paste("data/simu_cov_vmf_n",n,".RData", sep=""))
theta0 = c(1, -1, 1) / sqrt(3)
T_LS_eps0 = sapply(1:500, function(i){
    diff = (theta_LS_simu$eps0[i,] - theta0)[-1]
    Sigma = matrix(Sigma_LS_simu$eps0[,i],2)
    t(diff) %*% solve(Sigma) %*% diff
})
mean(T_LS_eps0 > qchisq(0.95, 2))
Sigma_MC_LS = cov(theta_LS_simu$eps0[,-1])
Sigma_asy_LS = matrix(rowMeans(Sigma_LS_simu$eps0),2)
round(norm((Sigma_MC_LS-Sigma_asy_LS) * n, "f"),4)

T_ESL_eps0 = sapply(1:500, function(i){
    diff = (theta_ESL_simu$eps0[i,] - theta0)[-1]
    Sigma = matrix(Sigma_ESL_simu$eps0[,i],2)
    t(diff) %*% solve(Sigma) %*% diff
})
mean(T_ESL_eps0 > qchisq(0.95, 2))
Sigma_MC_ESL = cov(theta_ESL_simu$eps0[,-1])
Sigma_asy_ESL = matrix(rowMeans(Sigma_ESL_simu$eps0),2)
round(norm((Sigma_MC_ESL-Sigma_asy_ESL) * n, "f"),4)

# epsilon = 10
T_LS_eps10 = sapply(1:500, function(i){
    diff = (theta_LS_simu$eps10[i,] - theta0)[-1]
    Sigma = matrix(Sigma_LS_simu$eps10[,i],2)
    t(diff) %*% solve(Sigma) %*% diff
})
mean(T_LS_eps10 > qchisq(0.95, 2))
Sigma_MC_LS = cov(theta_LS_simu$eps10[,-1])
Sigma_asy_LS = matrix(rowMeans(Sigma_LS_simu$eps10),2)
round(norm((Sigma_MC_LS-Sigma_asy_LS) * n, "f"),4)

T_ESL_eps10 = sapply(1:500, function(i){
    diff = (theta_ESL_simu$eps10[i,] - theta0)[-1]
    Sigma = matrix(Sigma_ESL_simu$eps10[,i],2)
    t(diff) %*% solve(Sigma) %*% diff
})
mean(T_ESL_eps10 > qchisq(0.95, 2))
Sigma_MC_ESL = cov(theta_ESL_simu$eps10[,-1])
Sigma_asy_ESL = matrix(rowMeans(Sigma_ESL_simu$eps10),2)
round(norm((Sigma_MC_ESL-Sigma_asy_ESL) * n, "f"),4)

# epsilon = 20
T_LS_eps20 = sapply(1:500, function(i){
    diff = (theta_LS_simu$eps20[i,] - theta0)[-1]
    Sigma = matrix(Sigma_LS_simu$eps20[,i],2)
    t(diff) %*% solve(Sigma) %*% diff
})
mean(T_LS_eps20 > qchisq(0.95, 2))
Sigma_MC_LS = cov(theta_LS_simu$eps20[,-1])
Sigma_asy_LS = matrix(rowMeans(Sigma_LS_simu$eps20),2)
round(norm((Sigma_MC_LS-Sigma_asy_LS) * n, "f"),4)

T_ESL_eps20 = sapply(1:500, function(i){
    diff = (theta_ESL_simu$eps20[i,] - theta0)[-1]
    Sigma = matrix(Sigma_ESL_simu$eps20[,i],2)
    t(diff) %*% solve(Sigma) %*% diff
})
mean(T_ESL_eps20 > qchisq(0.95, 2))
Sigma_MC_ESL = cov(theta_ESL_simu$eps20[,-1])
Sigma_asy_ESL = matrix(rowMeans(Sigma_ESL_simu$eps20),2)
round(norm((Sigma_MC_ESL-Sigma_asy_ESL) * n, "f"),4)
