## Simulations for shape of mean functions in the ESIM paper. 
## last change: Nov.23.2024
##
rm(list = ls())
source("ESLFns.R")

library(Directional)
library(mvtnorm)
library(foreach)
library(doSNOW)

# 1.shape of link function ----
B = 500
p = 4
n= 150
theta0 = SpheNormalize(c(1,-1,1,-1))

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
clusterEvalQ(cl, { source("ESLFns.R") })
registerDoSNOW(cl)

para_LS_simu <- para_ESL_simu <- para_FSIM_simu <- list()
bias_LS_simu <- bias_ESL_simu <- bias_FSIM_simu <- list()
MSE_LS_simu <- MSE_ESL_simu <- MSE_FSIM_simu <- list()
time_LS_simu <- time_ESL_simu <- time_FSIM_simu <- list()

for (mu_type in c("mu1", "mu2", "mu3")) {
for (type in c("vmf", "esag")) {
    pb <- txtProgressBar(min = 0, max = B, style = 3)
    progress <- function(b) setTxtProgressBar(pb, b)
    
    results <- foreach(b = 1:B, .packages = c("mvtnorm","Directional"), .options.snow = list(progress=progress),
                       .errorhandling = "remove") %dopar% {
        set.seed(20 + b)
        ## generate random sample  ####
        X = sapply(1:p, function(i) runif(n, -2, 2))
        I = X %*% t(t(theta0)) # index value 
        
        if (mu_type == "mu1") {
            ## mean curve 1
            mu = cbind(2*I/(I^2+2), -2*I/(I^2+2), (I^2-2)/(I^2+2))
            mu = t(apply(mu, 1, SpheNormalize))
        } else if (mu_type == "mu2") {
            ## mean curve 2
            tt = pnorm(I)
            mu = cbind(sqrt(1-tt^2)*cos(pi*tt),
                       sqrt(1-tt^2)*sin(pi*tt),
                       tt)
            mu = t(apply(mu, 1, SpheNormalize))
        } else {
            ## mean curve 3
            mu = cbind(sin(2*pi*I),
                       cos(pi*I),
                       2*I/sqrt(I^2+2))
            mu = t(apply(mu, 1, SpheNormalize))
        }
        
        if (type == "vmf") {
            Y = t(apply(mu, 1, function(x){rvmf(1, x, k=20)}))
        } else {
            Y = t(apply(mu, 1, function(x){resag(1, x * 5, c(1,0.5))}))
        }
        
        # LS estimate
        LS_time  = bias_LS = MSE_LS = NA
        para_LS = rep(NA, p+1)
        tryCatch({
            LS_time = system.time({LS = ls_est(c(0.5, rep(1,p)), xdat = X, ydat = Y)})[3]
            para_LS = c(LS$bw, LS$theta)
            bias_LS = SpheGeoDist(theta0, LS$theta)
            
            LS_mu = t(apply(LS$mu, 1, SpheNormalize)) # fitted values
            MSE_LS = mean(acos(rowSums(mu * LS_mu))^2, na.rm = T)
        })
        
        # ESL estimate
        ESL_time = bias_ESL = MSE_ESL = NA
        para_ESL = rep(NA, p+2)
        tryCatch({
            ESL_time = system.time({ESL = esl_est(theta_init = LS$theta, bw_init = LS$bw, 
                                                  xdat = X, ydat = Y, delta = 0.2)})[3]
            para_ESL = c(ESL$theta, ESL$bw, ESL$lambda)
            bias_ESL = SpheGeoDist(theta0, ESL$theta)
            
            ESL_mu = t(apply(ESL$mu, 1, SpheNormalize))
            MSE_ESL = mean(acos(rowSums(mu * ESL_mu))^2, na.rm = T)
        })
        
        # FSIM estimate
        FSIM_time = bias_FSIM = MSE_FSIM = NA
        para_FSIM = rep(NA, p+1)
        tryCatch({
            FSIM_time = system.time({FSIM = fsim_est(xdat = X, ydat = Y, init = rep(1,p))})[3]
            para_FSIM = c(FSIM$bw, FSIM$theta)
            bias_FSIM = SpheGeoDist(theta0, FSIM$theta)
            
            FSIM_index = as.vector(X %*% FSIM$theta)
            FSIM_mu = t(sapply(1:n, function(r){LocSpheGeoReg(FSIM_index[-r], Y[-r,], FSIM_index[r], optns = list(bw=FSIM$bw, kernel="gauss"))}))
            MSE_FSIM = mean(acos(rowSums(mu * FSIM_mu))^2, na.rm = T)
            
        })
                
        list(
            para_LS = para_LS, bias_LS = bias_LS, MSE_LS = MSE_LS, LS_time = LS_time,
            para_ESL = para_ESL, bias_ESL = bias_ESL, MSE_ESL = MSE_ESL, ESL_time = ESL_time,
            para_FSIM = para_FSIM, bias_FSIM = bias_FSIM, MSE_FSIM = MSE_FSIM, FSIM_time = FSIM_time
        )
    }
    ## save results 
    para_ESL_simu[[paste(mu_type, type, sep = "_")]] <- t(sapply(results, function(x) x$para_ESL))
    para_LS_simu[[paste(mu_type, type, sep = "_")]] <- t(sapply(results, function(x) x$para_LS))
    para_FSIM_simu[[paste(mu_type, type, sep = "_")]] <- t(sapply(results, function(x) x$para_FSIM))
    
    bias_ESL_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$bias_ESL)
    bias_LS_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$bias_LS)
    bias_FSIM_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$bias_FSIM)
    
    MSE_ESL_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$MSE_ESL)
    MSE_LS_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$MSE_LS)
    MSE_FSIM_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$MSE_FSIM)
    
    time_ESL_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$ESL_time)
    time_LS_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$LS_time)
    time_FSIM_simu[[paste(mu_type, type, sep = "_")]] <- sapply(results, function(x) x$FSIM_time)
}
}
save(para_LS_simu, para_ESL_simu, para_FSIM_simu,
     bias_LS_simu, bias_ESL_simu, bias_FSIM_simu,
     MSE_LS_simu, MSE_ESL_simu, MSE_FSIM_simu,
     time_LS_simu, time_ESL_simu, time_FSIM_simu,
     file = "./data/simu_shape.RData")


# figures ----
rm(list = ls())
load("./data/simu_shape.RData")
## estimated theta ####
for(type in c("vmf", "esag")){
    for (mu in c("mu1","mu2","mu3")) {
        d1 = data.frame(
            theta=c(para_LS_simu[[paste(mu, "_", type, sep="")]][,2], 
                    para_LS_simu[[paste(mu, "_", type, sep="")]][,3], 
                    para_LS_simu[[paste(mu, "_", type, sep="")]][,4], 
                    para_LS_simu[[paste(mu, "_", type, sep="")]][,5]),
            theta_inx = rep(c("1","2","3","4"), each=500),
            index = rep("LS",500*4))
        d2 = data.frame(
            theta=c(para_ESL_simu[[paste(mu, "_", type, sep="")]][,1], 
                    para_ESL_simu[[paste(mu, "_", type, sep="")]][,2], 
                    para_ESL_simu[[paste(mu, "_", type, sep="")]][,3], 
                    para_ESL_simu[[paste(mu, "_", type, sep="")]][,4]),
            theta_inx = rep(c("1","2","3","4"), each=500),
            index = rep("ESL",500*4))
        d3 = data.frame(
            theta=c(para_FSIM_simu[[paste(mu, "_", type, sep="")]][,2], 
                    para_FSIM_simu[[paste(mu, "_", type, sep="")]][,3], 
                    para_FSIM_simu[[paste(mu, "_", type, sep="")]][,4], 
                    para_FSIM_simu[[paste(mu, "_", type, sep="")]][,5]),
            theta_inx = rep(c("1","2","3","4"), each=500),
            index = rep("FSIM",500*4))
        d_ = rbind(d1,d2,d3)
        d_$index = factor(d_$index, levels = c("LS", "ESL", "FSIM"))
        d_$theta_inx = factor(d_$theta_inx, levels = c("1","2","3","4"))
        
        png(paste("./figures/shape/", mu,"_", type, ".png", sep = ""), width = 375, height = 275)
        par(pin=c(6,4), mar=c(4,5,1,1))
        boxplot(theta~index+theta_inx, data=d_, col=rainbow(3,s=0.35),
                at=c(1:3, 5:7, 9:11, 13:15),
                names=c("",expression(theta[1]),"",
                        "",expression(theta[2]),"",
                        "",expression(theta[3]),"",
                        "",expression(theta[4]),""),
                xlab = "", ylab = expression(hat(theta)),
                las = 1, cex.lab=1.2, cex.axis=1.2, ylim=c(-1,1)
        )
        abline(h=0.5, col="blue", lty=2, lwd=1)
        abline(h=-0.5, col="blue", lty=2, lwd=1)
        legend("bottomright", legend = levels(d_$index), fill = rainbow(3,s=0.35), cex = 1, ncol=3)
        dev.off()
    }
}

## bias ####
theta0 = c(.5, -.5, .5, -.5)
for (type in c("vmf", "esag")) {
    bias = data.frame(value = numeric(0),
                      link = character(0),
                      method = character(0))
    
    for (mu in c("mu1", "mu2", "mu3")) {
        bias_LS = apply(para_LS_simu[[paste(mu, "_", type, sep="")]], 1, function(x){acos(sum(x[2:5] * theta0))^2}) 
        bias_ESL = apply(para_ESL_simu[[paste(mu, "_", type, sep="")]], 1, function(x){acos(sum(x[1:4] * theta0))^2})
        bias_FSIM = apply(para_FSIM_simu[[paste(mu, "_", type, sep="")]], 1, function(x){acos(sum(x[2:5] * theta0))^2}) 
        
        df1 = data.frame(
            value = log(bias_LS),
            link = rep(mu, length(bias_LS)),
            method = rep("LS", length(bias_LS))
        )
        df2 = data.frame(
            value = log(bias_ESL),
            link = rep(mu, length(bias_ESL)),
            method = rep("ESL", length(bias_ESL))
        )
        df3 = data.frame(
            value = log(bias_FSIM),
            link = rep(mu, length(bias_FSIM)),
            method = rep("FSIM", length(bias_FSIM))
        )
        
        bias = rbind(bias, df1, df2, df3)
    }
    
    bias$method = factor(bias$method, levels = c("LS", "ESL", "FSIM"))
    bias$link = factor(bias$link, levels = c("mu1", "mu2", "mu3"))
    
    png(paste("./figures/shape/bias_",type,".png", sep=""), width = 375, height = 275)
    par(pin=c(6,4), mar=c(4,5,1,1))
    boxplot(value~method+link, data=bias, col = rainbow(3, s=0.35),
            at = c(1:3, 5:7, 9:11),
            names = c("", expression(tilde(mu)[1]), "", 
                      "", expression(tilde(mu)[2]), "",
                      "", expression(tilde(mu)[3]), ""),
            xlab = "", ylab = expression(paste("Log(bias) of ",hat(beta))),
            las = 1, cex.lab = 1.2, cex.axis = 1.5, ylim = c(-15, 2)
    )
    legend("bottomleft", legend = levels(bias$method), fill = rainbow(3, s=0.35), ncol = 3, cex=.85)
    dev.off()
}

## MSE ####
for (type in c("vmf", "esag")) {
    MSE = data.frame(value = numeric(0),
                      link = character(0),
                      method = character(0))
    
    for (mu in c("mu1", "mu2", "mu3")) {
        MSE_LS = MSE_LS_simu[[paste(mu, "_", type, sep="")]]
        MSE_ESL = MSE_ESL_simu[[paste(mu, "_", type, sep="")]]
        MSE_FSIM = MSE_FSIM_simu[[paste(mu, "_", type, sep="")]]
        
        df1 = data.frame(
            value = log(MSE_LS),
            link = rep(mu, length(MSE_LS)),
            method = rep("LS", length(MSE_LS))
        )
        df2 = data.frame(
            value = log(MSE_ESL),
            link = rep(mu, length(MSE_ESL)),
            method = rep("ESL", length(MSE_ESL))
        )
        df3 = data.frame(
            value = log(MSE_FSIM),
            link = rep(mu, length(MSE_FSIM)),
            method = rep("FSIM", length(MSE_FSIM))
        )
        
        MSE = rbind(MSE, df1, df2, df3)
    }
    
    MSE$method = factor(MSE$method, levels = c("LS", "ESL", "FSIM"))
    MSE$link = factor(MSE$link, levels = c("mu1", "mu2", "mu3"))
    
    png(paste("./figures/shape/MSE_",type,".png", sep=""), width = 375, height = 275)
    par(pin=c(6,4), mar=c(4,5,1,1))
    boxplot(value~method+link, data=MSE, col = rainbow(3, s=0.35),
            at = c(1:3, 5:7, 9:11),
            names = c("", expression(tilde(mu)[1]), "", 
                      "", expression(tilde(mu)[2]), "",
                      "", expression(tilde(mu)[3]), ""),
            xlab = "", ylab = expression(paste("Log(MSE)")),
            las = 1, cex.lab = 1.2, cex.axis = 1.5, ylim = c(-6.5, 1.2)
    )
    legend("bottomright", legend = levels(MSE$method), fill = rainbow(3, s=0.35), ncol = 3, cex=.85)
    dev.off()
}

## time ####
type = "vmf" # type = "esag"
mean_mat = matrix(0, nrow = 3, ncol = 3)
rownames(mean_mat) = c("LS", "ESL", "FSIM")
colnames(mean_mat) = c("mu1", "mu2", "mu3")

for (j in 1:3) {
    mu = c("mu1", "mu2", "mu3")[j]
    
    time_LS   = time_LS_simu[[paste0(mu, "_", type)]]
    time_ESL  = time_ESL_simu[[paste0(mu, "_", type)]]
    time_FSIM = time_FSIM_simu[[paste0(mu, "_", type)]]
    
    mean_mat["LS",   j] = mean(time_LS)
    mean_mat["ESL",  j] = mean(time_ESL)
    mean_mat["FSIM", j] = mean(time_FSIM, na.rm = T)
}
round(mean_mat, 2)

for (type in c("vmf", "esag")) {
    mean_mat = matrix(0, nrow = 3, ncol = 3)
    rownames(mean_mat) = c("LS", "ESL", "FSIM")
    colnames(mean_mat) = c("mu1", "mu2", "mu3")
    
    for (j in 1:3) {
        mu = c("mu1", "mu2", "mu3")[j]
        
        time_LS   = time_LS_simu[[paste0(mu, "_", type)]]
        time_ESL  = time_ESL_simu[[paste0(mu, "_", type)]]
        time_FSIM = time_FSIM_simu[[paste0(mu, "_", type)]]
        
        mean_mat["LS",   j] = mean(log(time_LS))
        mean_mat["ESL",  j] = mean(log(time_ESL))
        mean_mat["FSIM", j] = mean(log(time_FSIM), na.rm = T)
    }
    
    cols = rainbow(3, s = 0.35)
    
    png(paste0("./figures/shape/time_", type, ".png"),
        , width = 375, height = 275)
    par(pin = c(6,4), mar = c(4,5,1,1))
    
    barplot(mean_mat, beside = TRUE,
            col = cols,
            ylim = range(mean_mat) * c(0, 1.1),
            names.arg = c(expression(tilde(mu)[1]),
                          expression(tilde(mu)[2]),
                          expression(tilde(mu)[3])),
            ylab = "Time in seconds",
            cex.axis = 1.4, cex.lab = 1.2)
    
    legend("topleft", legend = rownames(mean_mat),
           fill = cols, ncol = 3, cex = 0.9)
    dev.off()
}



