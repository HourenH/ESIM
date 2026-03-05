# heavy-tailed and elliptical distributions ####
rm(list = ls())
source("ESLFns.R")
library(Directional)
library(foreach)
library(doSNOW)

## scatter plot
library(rgl)
set.seed(5)

n = 500
mu = c(0,0.5,1)
mu = mu / sqrt(sum(mu^2))
Y_SvMF = simProject(kappa = 1, V = diag(c(5,1)), mu = mu, a1 = 10, n=n)
Y_Cauchy = rsespc(n, 5 * c(1,1,1), theta = c(1,0.5))
Y_Cauchy = t(apply(Y_Cauchy, 1, function(r){rotation(c(1,1,1)/sqrt(3), mu) %*% r}))

options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
points3d(Y_SvMF, alpha=0.5, col="black", size=3)
points3d(mu, alpha=0.8, col="red", size=8)
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000 ,   1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60) # 500*500

options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
points3d(Y_Cauchy, alpha=0.5, col="black", size=3)
points3d(mu, alpha=0.8, col="red", size=8)
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000 ,   1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60) # 500*500
##################
rm(list = ls())
source("ESLFns.R")

B = 500
n = 150
p = 3
theta0 = SpheNormalize(c(1,-1,1))
type = "svmf" 
# type = "cauchy"

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
clusterEvalQ(cl, { source("ESLFns.R") })
registerDoSNOW(cl)

pb <- txtProgressBar(min = 0, max = B, style = 3)
progress <- function(b) setTxtProgressBar(pb, b)

results <- foreach(b = 1:B, .packages = c("mvtnorm","Directional"),
                   .options.snow = list(progress=progress),
                   .errorhandling = "pass") %dopar% {
    
    if (type == "svmf") {
        set.seed(b+1000)
        X = cbind(runif(n, -2, 2), runif(n, -2, 2), runif(n, -2, 2))
        U = 1/(1 + exp(-X %*% theta0))
        mu = cbind(sqrt(1-U^2) * cos(pi*U),
                   sqrt(1-U^2) * sin(pi*U),
                   U)
        ## generate sample from SvMF
        Y_SvMF = simProject(kappa = 1, V = diag(c(5,1)), mu = c(1,0,0), a1 = 10, n=n)
        ## rotate
        Y = t(sapply(1:n, function(r_){
            H = rotation(c(1,0,0), mu[r_,])
            H %*% t(t(Y_SvMF[r_,]))
        }))
    } else {
        set.seed(b+2000)
        X = cbind(runif(n, -2, 2), runif(n, -2, 2), runif(n, -2, 2))
        U = 1/(1 + exp(-X %*% theta0))
        mu = cbind(sqrt(1-U^2) * cos(pi*U),
                   sqrt(1-U^2) * sin(pi*U),
                   U)
        ## generate sample from Cauchy
        Y_Cauchy = rsespc(n, 5 * c(1,1,1), theta = c(1,0.5))
        ## rotate
        Y = t(sapply(1:n, function(r_){
            H = rotation(c(1,1,1)/sqrt(3), mu[r_,])
            H %*% t(t(Y_Cauchy[r_,]))
        }))
    }
        
    para_LS = rep(NA, p+1); bias_LS = NA; MSE_LS = NA
    para_SIQR = rep(NA, p+1); bias_SIQR = NA; MSE_SIQR = NA
    para_ESL = rep(NA, p+2); bias_ESL = NA; MSE_ESL = NA
    para_FSIM = rep(NA, p+1); bias_FSIM = NA; MSE_FSIM = NA
    
    # 1. LS estimate #####
    tryCatch({
        LS = ls_est(c(1, rep(1,p)), xdat = X, ydat = Y)
        para_LS = c(LS$bw, LS$theta)
        bias_LS = SpheGeoDist(theta0, LS$theta) # bias
        LS_mu = t(apply(LS$mu, 1, SpheNormalize)) # fitted values
        MSE_LS = mean(acos(rowSums(mu * LS_mu))^2, na.rm = T)
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
    }, error = function(e){
        message(sprintf("SIQR estimation failed (iter=%d): %s", b, e$message))
    })
    
    # 3. ESL estimate ####
    tryCatch({
        ESL = esl_est(theta_init = LS$theta, bw_init = LS$bw, 
                      xdat = X, ydat = Y, delta = 0.2)
        para_ESL = c(ESL$theta, ESL$bw, ESL$lambda)
        bias_ESL = SpheGeoDist(theta0, ESL$theta)
        # MSE of ESL
        ESL_mu = ESL$mu
        ESL_mu = t(apply(ESL_mu, 1, SpheNormalize))
        MSE_ESL = mean(acos(rowSums(mu * ESL_mu))^2, na.rm = T)
    }, error = function(e){
        message(sprintf("ESL estimation failed (iter=%d): %s", b, e$message))
    })
    
    # 4. FSIM ####
    tryCatch({
        FSIM = fsim_est(xdat = X, ydat = Y, init = rep(1,p))
        para_FSIM = c(FSIM$bw, FSIM$theta)
        bias_FSIM = SpheGeoDist(theta0, FSIM$theta)
        # MSE of FSIM
        FSIM_index = as.vector(X %*% FSIM$theta)
        FSIM_mu = t(sapply(1:n, function(r){LocSpheGeoReg(FSIM_index[-r], Y[-r,], FSIM_index[r], optns = list(bw=FSIM$bw, kernel="gauss"))}))
        MSE_FSIM = mean(acos(rowSums(mu * FSIM_mu))^2, na.rm = T)
    }, error = function(e){
        message(sprintf("FSIM estimation failed (iter=%d): %s", b, e$message))
    })
    
    res = list(para_LS = para_LS, bias_LS = bias_LS, MSE_LS = MSE_LS,
               para_SIQR = para_SIQR, bias_SIQR = bias_SIQR, MSE_SIQR = MSE_SIQR, 
               para_ESL = para_ESL, bias_ESL = bias_ESL, MSE_ESL = MSE_ESL, 
               para_FSIM = para_FSIM, bias_FSIM = bias_FSIM, MSE_FSIM = MSE_FSIM)
    save(res, file = paste("data/", type, "_", "replicate", b,".RData", sep = ""))
}
close(pb)

rm(list = ls())
load("simu_heavy.RData")
## bias ####
bias = data.frame(value = numeric(0),
                  type = character(0),
                  method = character(0))

for (type in c("svmf", "cauchy")) {
    dfLS = data.frame(
        value = log(bias_LS[[type]]),
        type = rep(type, length(bias_LS)),
        method = rep("LS", length(bias_LS)))
    
    dfESL = data.frame(
        value = log(bias_ESL[[type]]),
        type = rep(type, length(bias_ESL)),
        method = rep("ESL", length(bias_ESL)))
    
    dfSIQR = data.frame(
        value = log(bias_SIQR[[type]]),
        type = rep(type, length(bias_SIQR)),
        method = rep("SIQR", length(bias_SIQR)))
    
    dfFSIM = data.frame(
        value = log(bias_FSIM[[type]]),
        type = rep(type, length(bias_FSIM)),
        method = rep("FSIM", length(bias_FSIM)))
        
    bias = rbind(bias, dfLS, dfESL, dfSIQR, dfFSIM)
}
    
bias$method = factor(bias$method, levels = c("LS", "ESL", "SIQR","FSIM"))
bias$type = factor(bias$type, levels = c("svmf", "cauchy"))
    
png("../figures/heavy/bias.png", width = 375, height = 275)
par(pin=c(6,4), mar=c(4,5,1,1))
boxplot(value~method+type, data=bias, col = rainbow(4, s=0.35),
        at = c(1:4, 6:9),
        xlab = "", ylab = expression(paste("Log(bias) of ",hat(beta))), xaxt = "n",
        las = 1, cex.lab = 1.2, cex.axis = 1.5)
axis(1, at = c(2, 7.5), labels = c("SvMF", "Cauchy"), tick = FALSE, cex.axis = 1.5)
legend("bottomleft", legend = levels(bias$method), fill = rainbow(4, s=0.35), ncol = 4, cex=.85)
dev.off()

## MSE ####
MSE = data.frame(value = numeric(0),
                  type = character(0),
                  method = character(0))

for (type in c("svmf", "cauchy")) {
    dfLS = data.frame(
        value = log(MSE_LS[[type]]),
        type = rep(type, length(MSE_LS)),
        method = rep("LS", length(MSE_LS)))
    
    dfESL = data.frame(
        value = log(MSE_ESL[[type]]),
        type = rep(type, length(MSE_ESL)),
        method = rep("ESL", length(MSE_ESL)))
    
    dfSIQR = data.frame(
        value = log(MSE_SIQR[[type]]),
        type = rep(type, length(MSE_SIQR)),
        method = rep("SIQR", length(MSE_SIQR)))
    
    dfFSIM = data.frame(
        value = log(MSE_FSIM[[type]]),
        type = rep(type, length(MSE_FSIM)),
        method = rep("FSIM", length(MSE_FSIM)))
    
    MSE = rbind(MSE, dfLS, dfESL, dfSIQR, dfFSIM)
}

MSE$method = factor(MSE$method, levels = c("LS", "ESL", "SIQR", "FSIM"))
MSE$type = factor(MSE$type, levels = c("svmf", "cauchy"))

png("../figures/heavy/MSE.png", width = 375, height = 275)
par(pin=c(6,4), mar=c(4,5,1,1))
boxplot(value~method+type, data=MSE, col = rainbow(4, s=0.35),
        at = c(1:4, 6:9),
        xlab = "", ylab = expression(paste("Log(MSE) of ",hat(beta))), xaxt = "n",
        las = 1, cex.lab = 1.2, cex.axis = 1.5)
axis(1, at = c(2, 7.5), labels = c("SvMF", "Cauchy"), tick = FALSE, cex.axis = 1.5)
legend("bottomleft", legend = levels(MSE$method), fill = rainbow(4, s=0.35), ncol = 4, cex=.85)
dev.off()
