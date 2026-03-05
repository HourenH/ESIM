rm(list = ls())
source("ESLFns.R")
library(foreach)
library(doSNOW)
library(Directional)
library(MASS)

R = 500 # Monte Carlo runs
B = 500 # Bootstrap runs
epsilon = 20 # epsilon = c(0, 10, 20)

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
clusterEvalQ(cl, { source("ESLFns.R") })
registerDoSNOW(cl)

n = 200
p = 4; d = 3; k = 20
type = "vmf" # "esag"
Wald_stat = list()

# for (C in seq(0, 0.25, 0.01)) {
for (C in c(0, 0.1, 0.2)) {
    
    theta0 = c(1, 1, C, C)
    theta0 = theta0 / sqrt(sum(theta0^2))
    theta_null = c(1, 1, 0, 0)
    theta_null = theta_null / sqrt(sum(theta_null^2))
    A = diag(p-1)[2:3,]
    
    T_stat = rep(NA, R)
    
    for (r in 96:R) {
        set.seed(k+epsilon+r+as.integer(C*100))
        # generate random sample #####
        X = sapply(1:p, function(i) runif(n, -2, 2))
        U = 1/(1 + exp(-X %*% t(t(theta0))))
        mu = cbind(sqrt(1-U^2) * cos(pi*U),
                   sqrt(1-U^2) * sin(pi*U),
                   U)
        
        loc = sample(1:n, n*epsilon/100)
        Y = t(apply(mu, 1, function(x){rvmf(1, x, k)})) 
        Y[loc,] = -mu[loc,] # antipodal contamination
        
        # ESL estimate with original sample #####
        LS = ls_est(c(1, rep(1,p)), xdat = X, ydat = Y) # use LS estimate to initialize ESL
        ESL = esl_est(theta_init = LS$theta, bw_init = LS$bw, xdat = X, ydat = Y, delta = 0.2)
        lambda = ESL$lambda
        
        # multiplier bootstrap
        pb <- txtProgressBar(min = 0, max = B, style = 3)
        progress <- function(b) setTxtProgressBar(pb, b)
        
        results <- foreach(b = 1:B, .options.snow = list(progress=progress)) %dopar% {
            set.seed(k+epsilon+b+r+as.integer(C*100))
            w = rexp(n, rate = 1)
            ESL_boot = esl_index_bw(theta_init = ESL$theta, bw_init = ESL$bw, xdat = X, ydat = Y, lambda = ESL$lambda, w = w)
            theta_boot = ESL_boot$theta[-1]
            list(theta=theta_boot)
        }
        close(pb)
        # Wald-test stats
        Sigma = cov(t(sapply(results, function(x){x$theta})))
        T_stat[r] = t(A %*% (ESL$theta[-1] - theta_null[-1])) %*% solve(A %*% Sigma %*% t(A)) %*% (A %*% (ESL$theta[-1] - theta_null[-1]))
        }
    print(C)
    Wald_stat[[as.character(C)]] = T_stat
}
save(Wald_stat, file = paste("data/simu_power_eps", epsilon, ".RData", sep=""))

stopCluster(cl)


mean(Wald_stat$`0` > qchisq(0.95,2), na.rm = T)
