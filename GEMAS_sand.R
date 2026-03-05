## application of ESIM to GEMAS sand compositions
## last revision: 14.Feb.2025

rm(list = ls())
gemas = read.csv("data/gemas_merge.csv")

colnames(gemas)
table(gemas$COUNTRY) 
## consider SWE
gemas.SWE = gemas[gemas$COUNTRY == "SWE", -3]
# gemas.SPA = gemas[gemas$COUNTRY == "SPA", -3]
rm(gemas)
# convert compositional data
Y = apply(gemas.SWE[,c("sand", "silt", "clay")], 1, function(x){sqrt(x/100)})
Y = t(apply(Y,2,function(y){y/sqrt(sum(y^2))}))
rownames(Y) = NULL

# 1.EDA ----
## 1.1 Soil textural composition ----
png("./figures/gemas/sand.png", width = 500, height = 450)
par(mfrow=c(2,2), mar=c(4,5,2,1))
plot(Y[,1], ylim = c(0,1), xlab = "index of data", ylab = "Y", pch=1, cex.lab=1.5, cex.axis = 1.5)
points(Y[,2], col="blue", pch=2)
points(Y[,3], col="red", pch=3)
plot(Y[,1], ylim = c(0,1), xlab = "index of data", ylab = expression(Y[1]), pch=1, cex.lab=1.5, cex.axis = 1.5)
symbols(104, Y[104,1], circles = 0.08, inches = 0.08, add = TRUE, fg="red")
text(104, Y[104,1]-0.07, labels = c(104))
plot(Y[,2], ylim = c(0,1), xlab = "index of data", ylab = expression(Y[2]), pch=1, cex.lab=1.5, cex.axis = 1.5)
symbols(100, Y[100,2], circles =0.08, inches = 0.08, add = TRUE, fg="red")
text(100, Y[100,2]-0.07, labels = c(100))
plot(Y[,3], ylim = c(0,1), xlab = "index of data", ylab = expression(Y[3]), pch=1, cex.lab=1.5, cex.axis = 1.5)
symbols(99, Y[99,3], circles = 0.08, inches = 0.08, add = TRUE, fg="red")
text(99, Y[99,3]+0.07, labels = c(99))
dev.off()

## 1.2 categorical X ----
gemas.SWE$soilclass = factor(gemas.SWE$soilclass, levels = c("ll", "l", "m","s"))
table(gemas.SWE$soilclass)
which(gemas.SWE$soilclass == "s")
gemas.SWE = gemas.SWE[-99,]; Y = Y[-99,]
dim(gemas.SWE); dim(Y)
gemas.SWE$soilclass = factor(gemas.SWE$soilclass)
table(gemas.SWE$soilclass)
# require(robCompositions)
# par(mfrow=c(1,1))
# ternaryDiag(gemas.SWE[,c("sand", "silt", "clay", "soilclass")], 
#             col = as.numeric(gemas.SWE$soilclass), 
#             pch = as.numeric(gemas.SWE$soilclass)) 
# legend("topleft", c("l","ll","m"), col = 1:3, pch = 1:3)

for (i in c("sand", "silt", "clay")) {
    png(paste("./figures/gemas/box_", i, ".png",sep = ""), width = 400, height = 350)
    par(mfrow=c(1,1), pin=c(6,4), mar=c(4,5,1,1))
    boxplot(gemas.SWE[,i]~gemas.SWE$soilclass, 
            ylab = i, xlab = "soilclass", col=rainbow(3, s=0.35),
            cex.lab=2, cex.axis=1.5)
    dev.off()
}

## 1.3 continuous X ---- 
# "longitude" "latitude" "MeanTemp" "AnnPrec" "cec" "ph_cacl2" "toc"
for (i in c("MeanTemp", "AnnPrec", "cec", "ph_cacl2", "toc")) {
    png(paste("./figures/gemas/hist_", i, ".png",sep = ""), width = 300, height = 200)
    par(pin=c(6,4), mar=c(4,5,1,1))
    hist(gemas.SWE[,i], main="", xlab=i, cex.lab=1.5, cex.axis=1.2)
    dev.off()
}

gemas.SWE$lAnnPrec = log(gemas.SWE$AnnPrec)
gemas.SWE$ltoc = log(gemas.SWE$toc)

for (i in c("lAnnPrec", "ltoc")) {
    png(paste("./figures/gemas/hist_", i, ".png",sep = ""), width = 300, height = 200)
    par(pin=c(6,4), mar=c(4,5,1,1))
    hist(gemas.SWE[,i], main="", xlab=i, cex.lab=1.5, cex.axis=1.2)
    dev.off()
}

# pairs(cbind(Y,gemas.SWE[,c("MeanTemp", "cec", "ph_cacl2", "lAnnPrec","ltoc")]))

for (i in c("ph_cacl2", "ltoc")) {
    png(paste("./figures/gemas/scatter_", i, ".png",sep = ""), width = 300, height = 200)
    par(pin=c(6,4), mar=c(4,5,1,1))
    plot(gemas.SWE[,i], gemas.SWE[,"cec"], xlab = i, ylab = "cec",
         cex.lab=1.5, cex.axis=1.2, pch=19)
    dev.off()
}
cor(gemas.SWE[,"ph_cacl2"], gemas.SWE[,"cec"])
cor(gemas.SWE[,"ltoc"], gemas.SWE[,"cec"])

## standardize continuous X
X = apply(gemas.SWE[, c("MeanTemp", "lAnnPrec", "ltoc", "ph_cacl2")], 2,
          function(x){ (x - mean(x))/sd(x)})
# pairs(X)

# interactive effect
# for (j in c("MeanTemp", "cec", "ph_cacl2", "lAnnPrec","ltoc")) {
#     for (i in c("sand", "silt", "clay")) {
#         par(mfrow=c(1,1), pin=c(6,4), mar=c(4,5,1,1))
#         plot(gemas.SWE[,j], gemas.SWE[,i], xlab=j, ylab=i, type = "n")
#         points(gemas.SWE[gemas.SWE$soilclass=="l",j],
#                gemas.SWE[gemas.SWE$soilclass=="l",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="l","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="l","soilclass"]))
#         points(gemas.SWE[gemas.SWE$soilclass=="ll",j],
#                gemas.SWE[gemas.SWE$soilclass=="ll",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="ll","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="ll","soilclass"]))
#         points(gemas.SWE[gemas.SWE$soilclass=="m",j],
#                gemas.SWE[gemas.SWE$soilclass=="m",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="m","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="m","soilclass"]))
#         points(gemas.SWE[gemas.SWE$soilclass=="s",j],
#                gemas.SWE[gemas.SWE$soilclass=="s",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="s","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="s","soilclass"]))
#         
#     }
# }
# 
# for (j in c("MeanTemp", "cec", "ph_cacl2", "lAnnPrec","ltoc")) {
#     for (i in 1:3) {
#         par(mfrow=c(1,1), pin=c(6,4), mar=c(4,5,1,1))
#         plot(gemas.SWE[,j], Y[,i], xlab=j, ylab=i, type = "n")
#         points(gemas.SWE[gemas.SWE$soilclass=="l",j],
#                Y[gemas.SWE$soilclass=="l",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="l","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="l","soilclass"]))
#         points(gemas.SWE[gemas.SWE$soilclass=="ll",j],
#                Y[gemas.SWE$soilclass=="ll",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="ll","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="ll","soilclass"]))
#         points(gemas.SWE[gemas.SWE$soilclass=="m",j],
#                Y[gemas.SWE$soilclass=="m",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="m","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="m","soilclass"]))
#         points(gemas.SWE[gemas.SWE$soilclass=="s",j],
#                Y[gemas.SWE$soilclass=="s",i],
#                col=as.numeric(gemas.SWE[gemas.SWE$soilclass=="s","soilclass"]),
#                pch=as.numeric(gemas.SWE[gemas.SWE$soilclass=="s","soilclass"]))
#         
#     }
# }

# encode categorical X, "l" is the reference level
gemas.SWE$soilclass = factor(gemas.SWE$soilclass, levels = c("l","ll","m"))
X = cbind(model.matrix(~gemas.SWE$soilclass)[,-1], X)
colnames(X) = c("ll","m","MeanTemp", "lAnnPrec","ltoc","ph_cacl2")
head(X)
# create interactive term

rm(gemas.SWE,i)

# 2.Fit model ----
source("ESLFns.R")

## 2.1 10-fold CV for MSE and MSPE ----
MSE_LS_cv <- MSPE_LS_cv <- MSE_ESL_cv <- MSPE_ESL_cv <- MSE_FSIM_cv <- MSPE_FSIM_cv <- MSE_SIQR_cv <- MSPE_SIQR_cv <- rep(0, 10)
theta_LS.cv <- theta_ESL.cv <- theta_FSIM.cv <- theta_SIQR.cv <- matrix(0, nrow = 10, ncol = ncol(X))
h_LS.cv <- h_ESL.cv <- h_FSIM.cv <- h_SIQR.cv <- lambda.cv <- rep(0,10)
p = dim(X)[2]
n = nrow(X)
d = ncol(Y)

set.seed(1)

id = sample(1:n)
for (fold in 1:10) {
    test_id <- id[((fold-1) * (n%/%10) + 1):(fold * (n%/%10))]
    X_test = X[test_id,]
    Y_test = Y[test_id,]
    X_train = X[-test_id,]
    Y_train = Y[-test_id,]
    
    # fit ESIM(LS)
    LS_cv = ls_est(c(0.5, rep(1,p)), xdat = X_train, ydat = Y_train)
    h_LS.cv[fold] = LS_cv$bw
    theta_LS.cv[fold,] = LS_cv$theta
    ## MSE(LS)
    LS_mu = t(apply(LS_cv$mu, 1, SpheNormalize))
    MSE_LS_cv[fold] = mean(acos(rowSums(LS_mu * Y_train))^2, na.rm=T)
    ## MSPE(LS)
    LS_index = X_train %*% t(t(LS_cv$theta))
    LS_pred_index = X_test %*% t(t(LS_cv$theta))
    LS_pred = t(sapply(LS_pred_index, ll_LS, txdat = LS_index, tydat = Y_train, bw = LS_cv$bw))
    LS_mu_pred = t(apply(LS_pred[,1:3], 1, SpheNormalize))
    MSPE_LS_cv[fold] = mean(acos(rowSums(LS_mu_pred * Y_test))^2, na.rm = T)
    
    # fit ESIM(ESL)
    ESL_cv = esl_est(theta_init = LS_cv$theta, bw_init = LS_cv$bw, 
                  xdat = X_train, ydat = Y_train, delta = 0.2)
    h_ESL.cv[fold] = ESL_cv$bw
    theta_ESL.cv[fold,] = ESL_cv$theta
    ## MSE(ESL)
    ESL_mu = t(apply(ESL_cv$mu, 1, SpheNormalize))
    MSE_ESL_cv[fold] = mean(acos(rowSums(ESL_mu * Y_train))^2, na.rm = T)
    ## MSPE(ESL)
    ESL_index = as.vector(X_train %*% ESL_cv$theta)
    ESL_pred_index = as.vector(X_test %*% ESL_cv$theta)
    ESL_pred = sapply(1:length(test_id), function(i){
        ll_ESL_predict(txdat = ESL_index, tydat = Y_train, exdat = ESL_pred_index[i], bw = ESL_cv$bw, lambda = ESL_cv$lambda, 
                       init_param = LS_pred[i,])$mu
    }) 
    ESL_mu_pred = t(apply(ESL_pred, 2, SpheNormalize))                 
    MSPE_ESL_cv[fold] =  mean(acos(rowSums(ESL_mu_pred * Y_test))^2, na.rm = T)
    
    # fit FSIM
    FSIM_cv = fsim_est(xdat = X_train, ydat = Y_train, init = rep(1,p))
    h_FSIM.cv[fold] = FSIM_cv$bw
    theta_FSIM.cv[fold,] = FSIM_cv$theta
    ## MSE(FSIM)
    MSE_FSIM_cv[fold] = mean(acos(rowSums(FSIM_cv$mu * Y_train))^2, na.rm = T)
    ## MSPE(FSIM)
    FSIM_index = as.vector(X_train %*% FSIM_cv$theta)
    FSIM_pred_index = as.vector(X_test %*% FSIM_cv$theta)
    FSIM_mu_pred =t(sapply(FSIM_pred_index, LocSpheGeoReg,
                           xin =FSIM_index, yin = Y_train, optns = list(bw=FSIM_cv$bw,kernel = "gauss")))
    MSPE_FSIM_cv[fold] = mean(acos(rowSums(FSIM_mu_pred * Y_test))^2, na.rm = T)
    
    # fit SIQR
    SIQR.fit = index.gamma(y=Y_train, xx=X_train, tau=0.5, gamma0 = theta_LS.cv[fold,], maxiter = 100, crit = 1e-3)
    # SIQR.fit = siqr_est(param = c(1,1,1), xdat = X, ydat = Y)
    theta_SIQR.cv[fold,] = SIQR.fit$theta
    h_SIQR.cv[fold] = SIQR.fit$bw
    ## MSE of SIQR
    SIQR_index = as.vector(X_train %*% theta_SIQR.cv[fold,])
    SIQR_mu = sapply(1:length(SIQR_index), function(r){
        lprq_mul(SIQR_index[r], SIQR_index[-r], Y_train[-r,], h=h_SIQR.cv[fold])
    })
    SIQR_mu = t(apply(SIQR_mu, 2, SpheNormalize))
    MSE_SIQR_cv[fold] = mean(acos(rowSums(SIQR_mu * Y_train))^2, na.rm = T)
    ## MSPE of SIQR
    SIQR.predict = sapply(as.vector(X_test %*% theta_SIQR.cv[fold,]),
                          lprq_mul, x = SIQR_index, y = Y_train, h = h_SIQR.cv[fold])
    SIQR.predict = t(apply(SIQR.predict, 2, SpheNormalize))
    MSPE_SIQR_cv[fold] = mean(acos(rowSums(SIQR.predict * Y_test))^2, na.rm = T)
}


## save cv results
colnames(theta_LS.cv)= colnames(theta_ESL.cv) = colnames(theta_FSIM.cv) = colnames(theta_SIQR.cv) = colnames(X)
# save.image("./data/GEMAS_sand.Rdata")

## 2.2 model with full samples ----
# fit ESIM(LS)
## init
init.LS = colMeans(theta_LS.cv)
init.LS = SpheNormalize(init.LS)
init.LS = c(mean(h_LS.cv),init.LS)
LS = ls_est(init.LS, xdat = X, ydat = Y)
h_LS = LS$bw
theta_LS = SpheNormalize(LS$theta)
LS_mu_full = t(apply(LS$mu, 1, SpheNormalize))
## MSE(LS)
MSE_LS = mean(acos(rowSums(LS_mu_full * Y))^2, na.rm = T)

# fit ESIM(ESL)
## init
init.ESL = colMeans(theta_ESL.cv)
init.ESL = SpheNormalize(init.ESL)
ESL = esl_est(init.ESL, bw_init = mean(h_ESL.cv), xdat = X, ydat = Y, delta = 0.2)
h_ESL= ESL$bw
theta_ESL = SpheNormalize(ESL$theta)
ESL_mu_full = t(apply(ESL$mu, 1, SpheNormalize))
## MSE(ESL)
MSE_ESL = mean(acos(rowSums(ESL_mu_full * Y))^2, na.rm = T)

# fit FSIM
init.FSIM = colMeans(theta_FSIM.cv)
init.FSIM = SpheNormalize(init.FSIM)
FSIM = fsim_est(xdat = X, ydat = Y, init = init.FSIM)
h_FSIM = FSIM$bw
theta_FSIM = SpheNormalize(FSIM$theta)
## MSE(FSIM)
MSE_FSIM = mean(acos(rowSums(FSIM$mu * Y))^2, na.rm = T)

# fit SIQR
SIQR = index.gamma(y=Y, xx=X, tau=0.5, gamma0 = colMeans(theta_SIQR.cv), maxiter = 100, crit = 1e-3)
theta_SIQR = SpheNormalize(SIQR$theta)
h_SIQR = SIQR$bw
## MSE(SIQR)
SIQR_index = X %*% t(t(theta_SIQR))
SIQR_mu = sapply(1:length(SIQR_index), function(r){
    lprq_mul(SIQR_index[r], SIQR_index[-r], Y[-r,], h=h_SIQR)
})
SIQR_mu = t(apply(SIQR_mu, 2, SpheNormalize))
MSE_SIQR = mean(acos(rowSums(SIQR_mu * Y))^2, na.rm = T)


result = matrix(c(theta_LS, h_LS, MSE_LS, 
                  theta_ESL, h_ESL, MSE_ESL,
                  theta_FSIM, h_FSIM, MSE_FSIM,
                  theta_SIQR, h_SIQR, MSE_SIQR), byrow = T, nrow = 4)
colnames(result) = c(colnames(X), "h", "MSE")
rownames(result) = c("LS", "ESL", "FSIM", "SIQR")
result
# save.image("./data/GEMAS_sand.Rdata")

## 2.3 Bootstrap for SE of LS and ESL ----
nIter = 500
b_ESL.boot = b_LS.boot = matrix(0, nrow = nIter, ncol = ncol(X))
h_ESL.boot = h_LS.boot = rep(0, nIter)

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

set.seed(100)
for (b in 1:nIter) {
    # weighted bootstrap for LS
    w = rexp(n, rate = 1)
    w = n * w / sum(w)
    LS_boot = ls_est_weight(c(LS$bw, LS$theta), xdat = X, ydat = Y, weight = w)
    b_LS.boot[b,] = LS_boot$theta; h_LS.boot[b] = LS_boot$bw
    # weighted bootstrap for ESL
    w = rexp(n, rate = 1)
    w = n * w / sum(w)
    ESL_boot = esl_index_bw(theta_init = ESL$theta, bw_init = ESL$bw, xdat = X, ydat = Y, lambda = ESL$lambda, w = w)
    b_ESL.boot[b,] = ESL_boot$theta; h_ESL.boot[b] = ESL_boot$bw
}
colnames(b_LS.boot) = colnames(b_ESL.boot) = colnames(X)
# save.image("./data/GEMAS_sand.Rdata")


## 2.4 summary stats ----
rm(list = ls())
load("GEMAS_sand.Rdata")
round(result, 4)

round(mean(MSPE_LS_cv), 4) 
round(mean(MSPE_ESL_cv), 4) 
round(mean(MSPE_FSIM_cv, na.rm = T), 4)
round(mean(MSPE_SIQR_cv), 4)

round(apply(b_LS.boot, 2, sd),3)
round(apply(b_ESL.boot, 2, sd),3)

# round(sd(h_LS.boot, na.rm = T), 3)
# round(sd(h_ESL.boot[which(h_ESL.boot<1)], na.rm = T), 3) # remove unusual runs
# round(sd(h_FSIM.boot, na.rm = T), 3)
# round(sd(h_SIQR.boot, na.rm = T), 3)

A = diag(p-1)[3:5,]
theta_null = c(1,1,1,0,0,0); theta_null = SpheNormalize(theta_null)
Sigma_LS = cov(b_LS.boot[,-1])
Sigma_ESL = cov(b_ESL.boot[,-1])
T_LS = t(A %*% (LS$theta[-1] - theta_null[-1])) %*% solve(A %*% Sigma_LS %*% t(A)) %*% (A %*% (LS$theta[-1] - theta_null[-1]))
T_ESL = t(A %*% (ESL$theta[-1] - theta_null[-1])) %*% solve(A %*% Sigma_ESL %*% t(A)) %*% (A %*% (ESL$theta[-1] - theta_null[-1]))

# 3. Diagnostic figure ----
rm(list = ls())
load("GEMAS_sand.Rdata")
# rotated residuals for sampling, tangent space at the north pole c(1,0,0)
e_rot_LS = t(sapply(1:n, function(i){
    ## crude residuals
    u = LS_mu_full[i,]
    v = Y[i,]
    u = u / sqrt(sum(u^2))
    v = v / sqrt(sum(v^2))
    eps = v - sum(v * u) * u
    ## Rotation 
    R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
    
    return(R %*% t(t(eps))) 
}))

e_rot_ESL = t(sapply(1:n, function(i){
    ## project y at the tangent space of mu
    u = ESL_mu_full[i,]
    v = Y[i,]
    u = u / sqrt(sum(u^2))
    v = v / sqrt(sum(v^2))
    eps = v - sum(v * u) * u
    ## Rotation
    R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
    
    return(R %*% t(t(eps))) 
}))

e_rot_FSIM = t(sapply(1:n, function(i){
    ## crude residuals
    u = FSIM$mu[i,]
    v = Y[i,]
    u = u / sqrt(sum(u^2))
    v = v / sqrt(sum(v^2))
    
    eps = v - sum(v * u) * u
    ## Rotation
    R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
    
    return(R %*% t(t(eps))) 
}))

e_rot_SIQR = t(sapply(1:n, function(i){
    u = SIQR_mu[i,]
    v = Y[i,]
    u = u / sqrt(sum(u^2))
    v = v / sqrt(sum(v^2))
    eps = v - sum(v * u) * u
    ## Rotation
    R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
    
    return(R %*% t(t(eps))) 
}))

group_colors <- c("#1b9e77", "#d95f02", "#7570b3")
## 3.1 residuals plot for LS ----
png("../figures/gemas/sand_LS_res.png", width = 250, height = 250)
par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
plot(e_rot_LS[,2], e_rot_LS[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
     xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
     col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
dev.off()

LS_index_full = X %*% LS$theta
png("../figures/gemas/sand_LS_res2.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(LS_index_full, e_rot_LS[,2], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
points(LS_index_full[X[,"ll"]==1,], e_rot_LS[X[,"ll"]==1,2], 
       col = group_colors[1], pch=16, cex=1.5)
points(LS_index_full[X[,"m"]==1,], e_rot_LS[X[,"m"]==1,2], 
       col = group_colors[2], pch=17, cex=1.5)
points(LS_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
       e_rot_LS[(X[,"m"]==0)&(X[,"ll"]==0),2], 
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()

png("../figures/gemas/sand_LS_res3.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(LS_index_full, e_rot_LS[,3], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
points(LS_index_full[X[,"ll"]==1,], e_rot_LS[X[,"ll"]==1,3],
       col = group_colors[1], pch=16, cex=1.5)
points(LS_index_full[X[,"m"]==1,], e_rot_LS[X[,"m"]==1,3], 
       col = group_colors[2], pch=17, cex=1.5)
points(LS_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
       e_rot_LS[(X[,"m"]==0)&(X[,"ll"]==0),3],
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()


## 3.2 residuals plot for ESL ----
ESL_index_full = X %*% theta_ESL

png("../figures/gemas/sand_ESL_res.png", width = 250, height = 250)
par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
plot(e_rot_ESL[,2], e_rot_ESL[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
     xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
     col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
dev.off()

png("../figures/gemas/sand_ESL_res2.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(ESL_index_full, e_rot_ESL[,2], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
points(ESL_index_full[X[,"ll"]==1,], e_rot_ESL[X[,"ll"]==1,2], 
       col = group_colors[1], pch=16, cex=1.5)
points(ESL_index_full[X[,"m"]==1,], e_rot_ESL[X[,"m"]==1,2], 
       col = group_colors[2], pch=17, cex=1.5)
points(ESL_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
       e_rot_ESL[(X[,"m"]==0)&(X[,"ll"]==0),2], 
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()

png("../figures/gemas/sand_ESL_res3.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(ESL_index_full, e_rot_ESL[,3], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
points(ESL_index_full[X[,"ll"]==1,], e_rot_ESL[X[,"ll"]==1,3],
       col = group_colors[1], pch=16, cex=1.5)
points(ESL_index_full[X[,"m"]==1,], e_rot_ESL[X[,"m"]==1,3],
       col = group_colors[2], pch=17, cex=1.5)
points(ESL_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
       e_rot_ESL[(X[,"m"]==0)&(X[,"ll"]==0),3], 
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()

## 3.3 residuals plot for FSIM ----
FSIM_index_full = as.vector(X %*% FSIM$theta)
png("../figures/gemas/sand_FSIM_res.png", width = 250, height = 250)
par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
plot(e_rot_FSIM[,2], e_rot_FSIM[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
     xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
     col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
dev.off()

png("../figures/gemas/sand_FSIM_res2.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(FSIM_index_full, e_rot_FSIM[,2], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
points(FSIM_index_full[X[,"ll"]==1], e_rot_FSIM[X[,"ll"]==1,2], 
       col = group_colors[1], pch=16, cex=1.5)
points(FSIM_index_full[X[,"m"]==1], e_rot_FSIM[X[,"m"]==1,2], 
       col = group_colors[2], pch=17, cex=1.5)
points(FSIM_index_full[(X[,"m"]==0)&(X[,"ll"]==0)], 
       e_rot_FSIM[(X[,"m"]==0)&(X[,"ll"]==0),2], 
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()

png("../figures/gemas/sand_FSIM_res3.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(FSIM_index_full, e_rot_FSIM[,3], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
points(FSIM_index_full[X[,"ll"]==1], e_rot_FSIM[X[,"ll"]==1,3], 
       col = group_colors[1], pch=16, cex=1.5)
points(FSIM_index_full[X[,"m"]==1], e_rot_FSIM[X[,"m"]==1,3], 
       col = group_colors[2], pch=17, cex=1.5)
points(FSIM_index_full[(X[,"m"]==0)&(X[,"ll"]==0)], 
       e_rot_FSIM[(X[,"m"]==0)&(X[,"ll"]==0),3], 
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()

## 3.4 residuals plot for SIQR ----
SIQR_index_full = X %*% SIQR$theta

png("../figures/gemas/sand_SIQR_res.png", width = 250, height = 250)
par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
plot(e_rot_SIQR[,2], e_rot_SIQR[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
     xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
     col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
dev.off()

png("../figures/gemas/sand_SIQR_res2.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(SIQR_index_full, e_rot_SIQR[,2], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
points(SIQR_index_full[X[,"ll"]==1,], e_rot_SIQR[X[,"ll"]==1,2],
       col = group_colors[1], pch=16, cex=1.5)
points(SIQR_index_full[X[,"m"]==1,], e_rot_SIQR[X[,"m"]==1,2],
       col = group_colors[2], pch=17, cex=1.5)
points(SIQR_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
       e_rot_SIQR[(X[,"m"]==0)&(X[,"ll"]==0),2], 
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()

png("../figures/gemas/sand_SIQR_res3.png", width = 375, height = 250)
par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
plot(SIQR_index_full, e_rot_SIQR[,3], ylim = c(-.5,.5), 
     xlab = "", ylab = "", cex.lab = 1.5, type = "n")
mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
points(SIQR_index_full[X[,"ll"]==1,], e_rot_SIQR[X[,"ll"]==1,3], 
       col = group_colors[1], pch=16, cex=1.5)
points(SIQR_index_full[X[,"m"]==1,], e_rot_SIQR[X[,"m"]==1,3], 
       col = group_colors[2], pch=17, cex=1.5)
points(SIQR_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
       e_rot_SIQR[(X[,"m"]==0)&(X[,"ll"]==0),3],
       col = group_colors[3], pch=18, cex=1.5)
legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
dev.off()
