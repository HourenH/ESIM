## Simulation for contaminated responses in the ESIM paper.
## last change: 30.Sep.2025
##

rm(list = ls())
# source("ESLFns.R")
library(Directional)
library(mvtnorm)
library(foreach)
library(doSNOW)
library(rgl)

B = 500

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
clusterEvalQ(cl, { source("ESLFns.R") })
registerDoSNOW(cl)


for (type in c("perp", "oppo")) {
    for (k in c(5, 20, 50)) {
        para_LS_simu <- para_ESL_simu <- para_FSIM_simu <- para_SIQR_simu <- list()
        bias_LS_simu <- bias_ESL_simu <- bias_FSIM_simu <- bias_SIQR_simu <- list()
        MSE_LS_simu <- MSE_ESL_simu <- MSE_FSIM_simu <- MSE_SIQR_simu <- list()
        MSPE_LS_simu <- MSPE_ESL_simu <- MSPE_FSIM_simu <- MSPE_SIQR_simu <- list()
        for (epsilon in c(0,10,20,30)) {
            pb <- txtProgressBar(min = 0, max = B, style = 3)
            progress <- function(b) setTxtProgressBar(pb, b)
    
            results <- foreach(b = 1:B, .packages = c("mvtnorm","Directional"), .options.snow = list(progress=progress)) %dopar% {
                set.seed(k+b+epsilon)
                
                # generate random sample #####
                n = 100
                p = 3; d = 3
                theta0 = c(1, -1, 1)
                theta0 = SpheNormalize(theta0)
                X = cbind(runif(n, -2, 2),runif(n, -2, 2),runif(n, -2, 2))
                U = 1/(1 + exp(-X %*% t(t(theta0))))
                mu = cbind(sqrt(1-U^2) * cos(pi*U),
                           sqrt(1-U^2) * sin(pi*U),
                           U)
                
                Y = t(apply(mu, 1, function(x){rvmf(1, x, k)}))
                location = sample(1:n, n*epsilon/100)
                if(type == "perp"){
                    mu_perp = t(apply(mu, 1, function(u){
                        u_perp = c(1,0,0) - sum(u* c(1,0,0)) * u
                        SpheNormalize(u_perp)}))
                    Y[location,] = mu_perp[location,] # orthogonal contamination
                } else {
                    Y[location,] = -mu[location,] # antipodal contamination
                }
                
                X_test = cbind(runif(n, -2, 2),runif(n, -2, 2),runif(n, -2, 2))
                U_test = 1/(1 + exp(-X_test %*% t(t(theta0))))
                mu_test = cbind(sqrt(1 - U_test^2) * cos(pi * U_test),
                                sqrt(1 - U_test^2) * sin(pi * U_test),
                                U_test)
                
                para_LS = rep(NA, p+1); bias_LS = NA; MSE_LS = NA; MSPE_LS = NA
                para_SIQR = rep(NA, p+1); bias_SIQR = NA; MSE_SIQR = NA; MSPE_SIQR = NA
                para_ESL = rep(NA, p+2); bias_ESL = NA; MSE_ESL = NA; MSPE_ESL = NA
                para_FSIM = rep(NA, p+1); bias_FSIM = NA; MSE_FSIM = NA; MSPE_FSIM = NA
                
                # 1. LS estimate #####
                tryCatch({
                    LS = ls_est(c(1, rep(1,p)), xdat = X, ydat = Y)
                    para_LS = c(LS$bw, LS$theta)
                    bias_LS = SpheGeoDist(theta0, LS$theta) # bias
                    LS_mu = t(apply(LS$mu, 1, SpheNormalize)) # fitted values
                    MSE_LS = mean(acos(rowSums(mu * LS_mu))^2, na.rm = T)
                    LS_predict = t(sapply(as.vector(X_test %*% LS$theta), function(idx){
                        fit = ll_LS(txdat = X %*% LS$theta, tydat = Y, exdat = idx, bw = LS$bw)
                        SpheNormalize(fit[1:d])})) # predicted values
                    MSPE_LS = mean(acos(rowSums(mu_test * LS_predict))^2, na.rm = T)
                }, error = function(e){
                    message(sprintf("LS estimation failed (iter=%d): %s", b, e$message))
                })
                
                # 2. SIQR for each component of Y ####
                tryCatch({
                    SIQR = index.gamma(y=Y, xx=X, tau=0.5, gamma0 = c(1,1,1), maxiter = 100, crit = 1e-3)
                    para_SIQR = c(SIQR$bw, SIQR$theta)
                    bias_SIQR = SpheGeoDist(theta0, SIQR$theta) # bias
                    # MSE of SIQR
                    SIQR_index = X %*% SIQR$theta
                    SIQR_mu = t(sapply(1:n, function(r){
                        fit = lprq_mul(SIQR_index[r], SIQR_index[-r], Y[-r,], h=SIQR$bw)
                        SpheNormalize(fit)})) # fitted values
                    MSE_SIQR = mean(acos(rowSums(mu * SIQR_mu))^2, na.rm = T)
                    ## MSPE of SIQR
                    SIQR_predict = sapply(as.vector(X_test %*% SIQR$theta),
                                          lprq_mul, x = SIQR_index, y = Y, h = SIQR$bw)
                    SIQR_predict = t(apply(SIQR_predict, 2, SpheNormalize))
                    MSPE_SIQR = mean(acos(rowSums(mu_test * SIQR_predict))^2, na.rm = T)
                }, error = function(e){
                    message(sprintf("SIQR estimation failed (iter=%d): %s", b, e$message))
                })
                
                # 3. ESL estimate ####
                tryCatch({
                    ESL = esl_est(theta_init = LS$theta, bw_init = LS$bw, 
                                  xdat = X, ydat = Y, delta = 0.4)
                    para_ESL = c(ESL$theta, ESL$bw, ESL$lambda)
                    bias_ESL = SpheGeoDist(theta0, ESL$theta)
                    # MSE of ESL
                    ESL_mu = ESL$mu
                    ESL_mu = t(apply(ESL_mu, 1, SpheNormalize))
                    MSE_ESL = mean(acos(rowSums(mu * ESL_mu))^2, na.rm = T)
                    # MSPE of ESL
                    ESL_index_test = as.vector(X_test %*% ESL$theta)
                    ESL_predict = t(sapply(seq_len(nrow(X_test)), function(idx){
                        init_para = c(LS_predict[idx,], LS_predict[idx,])
                        fit = ll_ESL_predict(txdat = as.vector(X %*% ESL$theta), tydat = Y, exdat = ESL_index_test[idx], 
                                             bw = ESL$bw, lambda = ESL$lambda, init_param = init_para)
                        SpheNormalize(fit$mu)})) # predicted values
                    MSPE_ESL = mean(acos(rowSums(mu_test * ESL_predict))^2, na.rm = T)
                }, error = function(e){
                    message(sprintf("ESL estimation failed (iter=%d): %s", b, e$message))
                })
                
                
                # 4. FSIM ####
                tryCatch({
                    FSIM = fsim_est(xdat = X, ydat = Y, init = theta0)
                    para_FSIM = c(FSIM$bw, FSIM$theta)
                    bias_FSIM = SpheGeoDist(theta0, FSIM$theta)
                    # MSE of FSIM
                    FSIM_index = as.vector(X %*% FSIM$theta)
                    FSIM_mu = t(sapply(1:n, function(r){LocSpheGeoReg(FSIM_index[-r], Y[-r,], FSIM_index[r], optns = list(bw=FSIM$bw, kernel="gauss"))}))
                    MSE_FSIM = mean(acos(rowSums(mu * FSIM_mu))^2, na.rm = T)
                    # MSPE of FSIM
                    FSIM_predict =t(sapply(as.vector(X_test %*% FSIM$theta), LocSpheGeoReg,
                                           xin =FSIM_index, yin = Y, optns = list(bw=FSIM$bw,kernel = "gauss")))
                    MSPE_FSIM = mean(acos(rowSums(mu_test * FSIM_predict))^2, na.rm = T)
                }, error = function(e){
                    message(sprintf("FSIM estimation failed (iter=%d): %s", b, e$message))
                })
                

                list(
                    para_LS = para_LS, bias_LS = bias_LS, MSE_LS = MSE_LS, MSPE_LS = MSPE_LS,
                    para_SIQR = para_SIQR, bias_SIQR = bias_SIQR, MSE_SIQR = MSE_SIQR, MSPE_SIQR = MSPE_SIQR,
                    para_ESL = para_ESL, bias_ESL = bias_ESL, MSE_ESL = MSE_ESL, MSPE_ESL = MSPE_ESL,
                    para_FSIM = para_FSIM, bias_FSIM = bias_FSIM, MSE_FSIM = MSE_FSIM, MSPE_FSIM = MSPE_FSIM)
            }
            
            # save results #####
            para_ESL_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- t(sapply(results, function(x) x$para_ESL))
            para_LS_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- t(sapply(results, function(x) x$para_LS))
            para_FSIM_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- t(sapply(results, function(x) x$para_FSIM))
            para_SIQR_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- t(sapply(results, function(x) x$para_SIQR))
            
            bias_ESL_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$bias_ESL)
            bias_LS_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$bias_LS)
            bias_FSIM_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$bias_FSIM)
            bias_SIQR_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$bias_SIQR)
            
            MSE_ESL_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSE_ESL)
            MSE_LS_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSE_LS)
            MSE_FSIM_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSE_FSIM)
            MSE_SIQR_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSE_SIQR)
            
            MSPE_ESL_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSPE_ESL)
            MSPE_LS_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSPE_LS)
            MSPE_FSIM_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSPE_FSIM)
            MSPE_SIQR_simu[[paste("k", k, "_eps", epsilon, sep="")]] <- sapply(results, function(x) x$MSPE_SIQR)
        }
        save(para_LS_simu, para_ESL_simu, para_SIQR_simu, para_FSIM_simu,
             bias_LS_simu, bias_ESL_simu, bias_SIQR_simu, bias_FSIM_simu,
             MSE_LS_simu, MSE_ESL_simu, MSE_SIQR_simu, MSE_FSIM_simu,
             MSPE_LS_simu, MSPE_ESL_simu, MSPE_SIQR_simu, MSPE_FSIM_simu,
             file = paste("./data/simu_contaminate_", type, "_k_", k, ".RData", sep = ""))
}
}


# figures ####
rm(list = ls())
# k; concentration parameter, 5, 20, 50.
# type; type of outliers, either "perp" or "oppo"
for (k in c(5,20,50)){
    for (type in c("oppo","perp")) {
        load(paste("./data/simu_contaminate_", type, "_k_", k, ".RData", sep = ""))
        epsilon = c(0,10,20)
        
        ## boxplot for vmf bias ----
        bias = data.frame(value = numeric(0),
                          epsilon = character(0),
                          method = character(0))
        ## import data
        for (i in 1:length(epsilon)) {
            bias_ESL = bias_ESL_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
            bias_LS = bias_LS_simu[[paste("k", k, "_eps", epsilon[i], sep = "")]]
            bias_SIQR = bias_SIQR_simu[[paste("k", k, "_eps", epsilon[i], sep = "")]]
            bias_FSIM = bias_FSIM_simu[[paste("k", k, "_eps", epsilon[i], sep = "")]]
            df1 = data.frame(
                value = log(bias_ESL * sqrt(k)),
                epsilon = rep(paste(epsilon[i]), length(bias_ESL)),
                method = rep("ESL", length(bias_ESL))
            )
            df2 = data.frame(
                value = log(bias_LS * sqrt(k)),
                epsilon = rep(paste(epsilon[i]), length(bias_LS)),
                method = rep("LS", length(bias_LS))
            )
            df3 = data.frame(
                value = log(bias_SIQR * sqrt(k)),
                epsilon = rep(paste(epsilon[i]), length(bias_SIQR)),
                method = rep("SIQR", length(bias_SIQR))
            )
            df4 = data.frame(
                value = log(bias_FSIM * sqrt(k)),
                epsilon = rep(paste(epsilon[i]), length(bias_FSIM)),
                method = rep("FSIM", length(bias_FSIM))
            ) 
            bias = rbind(bias, df1, df2, df3, df4)
        }
        bias$method = factor(bias$method, levels = c("LS", "SIQR", "ESL", "FSIM"))

        
        png(paste("figures/contaminate/bias_", type, "_k_", k,".png", sep=""), width = 375, height = 275)
        par(pin=c(6,4), mar=c(4,5,1,1))
        boxplot(value~method+epsilon, data=bias, col = rainbow(4, s=0.35),
                names = rep("", 12), at = c(1:4, 6:9, 11:14),
                xlab = "", ylab = "Log(bias)",
                las = 1, cex.lab = 1.2, cex.axis = 1.2, ylim=c(-6.5, 4), xaxt = "n" )
        at_pos <- c(mean(1:4), mean(6:9), mean(11:14))
        axis(1, at = at_pos, labels = c(expression(epsilon == 0 * "%"),
                                        expression(epsilon == 10 * "%"),
                                        expression(epsilon == 20 * "%")), 
             tick = F, line = 0, cex.axis = 1.35)
        mtext(bquote(kappa == .(k)), side = 1, line = 2.2, cex = 1.2)
        legend("topleft", legend = c("LS", "SIQR", "ESL", "FSIM"), fill = rainbow(4, s=0.35), ncol = 4, cex=0.9)
        dev.off()
        
        ## boxplot for vmf MSE ####
        MSE = data.frame(value = numeric(0),
                         epsilon = character(0),
                         method = character(0))
        ## import data
        for (i in 1:length(epsilon)) {
            MSE_ESL = MSE_ESL_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
            MSE_LS = MSE_LS_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
            MSE_SIQR = MSE_SIQR_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
            MSE_FSIM = MSE_FSIM_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
            
            df1 = data.frame(
                value = log(MSE_ESL * k),
                epsilon = rep(paste(epsilon[i]), length(MSE_ESL)),
                method = rep("ESL", length(MSE_ESL))
            )
            df2 = data.frame(
                value = log(MSE_LS * k),
                epsilon = rep(paste(epsilon[i]), length(MSE_LS)),
                method = rep("LS", length(MSE_LS))
            )
            df3 = data.frame(
                value = log(MSE_SIQR * k),
                epsilon = rep(paste(epsilon[i]), length(MSE_SIQR)),
                method = rep("SIQR", length(MSE_SIQR))
            )
            df4 = data.frame(
                value = log(MSE_FSIM * k),
                epsilon = rep(paste(epsilon[i]), length(MSE_FSIM)),
                method = rep("FSIM", length(MSE_FSIM))
            )
            MSE = rbind(MSE, df1, df2, df3, df4)
        }
        MSE$method = factor(MSE$method, levels=c("LS", "SIQR", "ESL", "FSIM"))
        
        png(paste("figures/contaminate/mse_", type, "_k_", k,".png", sep=""), width = 375, height = 275)
        par(pin=c(6,4), mar=c(4,5,1,1))
        boxplot(value~method+epsilon, data=MSE, col = rainbow(4, s=0.35),
                names = rep("", 12), at = c(1:4, 6:9, 11:14),
                xlab = "", ylab = "Log(MSE)",
                las = 1, cex.lab = 1.2, ylim=c(-3.5,4), cex.axis = 1.2)
        axis(1, at = at_pos, labels = c(expression(epsilon == 0 * "%"),
                                        expression(epsilon == 10 * "%"),
                                        expression(epsilon == 20 * "%")), 
             tick = F, line = 0, cex.axis = 1.35)
        mtext(bquote(kappa == .(k)), side = 1, line = 2.2, cex = 1.2)
        legend("topleft", legend = c("LS", "SIQR", "ESL", "FSIM"), fill = rainbow(4, s=0.35), ncol = 4, cex=.85)
        dev.off()
    }
}

## summary stat for vmf bias ####
epsilon = c(0, 10, 20)
# k; concentration parameter, 5, 20, 50.
# type; type of outliers, either "perp" or "oppo"
type = "oppo"
k = 50 # 5, 20, 50
load(paste("data/simu_contaminate_", type, "_k_", k, ".RData", sep=""))

bias = data.frame(value = numeric(0),
                  epsilon = character(0),
                  method = character(0))
## import data
for (i in 1:length(epsilon)) {
    bias_ESL = bias_ESL_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
    bias_LS = bias_LS_simu[[paste("k", k, "_eps", epsilon[i], sep = "")]]
    bias_SIQR = bias_SIQR_simu[[paste("k", k, "_eps", epsilon[i], sep = "")]]
    bias_FSIM = bias_FSIM_simu[[paste("k", k, "_eps", epsilon[i], sep = "")]]
    df1 = data.frame(
        value = bias_ESL,
        epsilon = rep(paste(epsilon[i]), length(bias_ESL)),
        method = rep("ESL", length(bias_ESL))
    )
    df2 = data.frame(
        value = bias_LS,
        epsilon = rep(paste(epsilon[i]), length(bias_LS)),
        method = rep("LS", length(bias_LS))
    )
    df3 = data.frame(
        value = bias_SIQR,
        epsilon = rep(paste(epsilon[i]), length(bias_SIQR)),
        method = rep("SIQR", length(bias_SIQR))
    )
    df4 = data.frame(
        value = bias_FSIM,
        epsilon = rep(paste(epsilon[i]), length(bias_FSIM)),
        method = rep("FSIM", length(bias_FSIM))
    ) 
    bias = rbind(bias, df1, df2, df3, df4)
}

bias$method = factor(bias$method, levels = c("LS", "SIQR", "ESL", "FSIM"))
bias$value = bias$value * sqrt(k)

bias_ave = aggregate(value~epsilon+method, data = bias, FUN = mean)
bias_sd = aggregate(value~epsilon+method, data = bias, FUN = sd)
bias_ave; bias_sd

MSE = data.frame(value = numeric(0),
                 epsilon = character(0),
                 method = character(0))
## import data
for (i in 1:length(epsilon)) {
    MSE_ESL = MSE_ESL_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
    MSE_LS = MSE_LS_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
    MSE_SIQR = MSE_SIQR_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
    MSE_FSIM = MSE_FSIM_simu[[paste("k", k,"_eps", epsilon[i], sep = "")]]
    
    df1 = data.frame(
        value = MSE_ESL,
        epsilon = rep(paste(epsilon[i]), length(MSE_ESL)),
        method = rep("ESL", length(MSE_ESL))
    )
    df2 = data.frame(
        value = MSE_LS,
        epsilon = rep(paste(epsilon[i]), length(MSE_LS)),
        method = rep("LS", length(MSE_LS))
    )
    df3 = data.frame(
        value = MSE_SIQR,
        epsilon = rep(paste(epsilon[i]), length(MSE_SIQR)),
        method = rep("SIQR", length(MSE_SIQR))
    )
    df4 = data.frame(
        value = MSE_FSIM,
        epsilon = rep(paste(epsilon[i]), length(MSE_FSIM)),
        method = rep("FSIM", length(MSE_FSIM))
    )
    MSE = rbind(MSE, df1, df2, df3, df4)
}
MSE$method = factor(MSE$method, levels=c("LS", "SIQR", "ESL", "FSIM"))
MSE$value = MSE$value * k

MSE_ave = aggregate(value~epsilon+method, data = MSE, FUN = mean)
MSE_sd = aggregate(value~epsilon+method, data = MSE, FUN = sd)
MSE_ave; MSE_sd
