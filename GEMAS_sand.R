## application of ESIM to GEMAS sand compositions
## last revision: 2.Jun.2026

rm(list = ls())
gemas = read.csv("data/gemas_merge.csv")

colnames(gemas)
table(gemas$COUNTRY) 
## consider FRA
gemas.FRA = gemas[gemas$COUNTRY == "FRA", -3]
rm(gemas)
# convert compositional data
Y = apply(gemas.FRA[,c("sand", "silt", "clay")], 1, function(x){sqrt(x/100)})
Y = t(apply(Y,2,function(y){y/sqrt(sum(y^2))}))
rownames(Y) = NULL

# 1.EDA ----
## 1.1 Soil textural composition ----
png("./figures/gemas/lon_comp.png", width = 500, height = 450)
par(mfrow=c(2,2), mar=c(4,5,2,1))
plot(gemas.FRA$longitude, Y[,1], ylim = c(0,1), xlab = "longitude", ylab = "Y", pch=1, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$longitude, Y[,2], col="blue", pch=2)
points(gemas.FRA$longitude, Y[,3], col="red", pch=3)
plot(gemas.FRA$longitude, Y[,1], ylim = c(0,1), xlab = "longitude", ylab = expression(Y[1]), pch=1, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$longitude[207], Y[207,1], col = "red", pch = 1, cex = 2)
text(gemas.FRA$longitude[207], Y[207,1], labels = "207", col = "red", font = 1.5, pos = 1, cex = 1.2)

plot(gemas.FRA$longitude, Y[,2], ylim = c(0,1), xlab = "longitude", ylab = expression(Y[2]), pch=2, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$longitude[c(22,205)], Y[c(22,205),2], col = "red", pch = 1, cex = 2)
text(gemas.FRA$longitude[c(22,205)], Y[c(22,205),2], labels = c("22","205"), col = "red", font = 1.5, pos = 1, cex = 1.2)

plot(gemas.FRA$longitude, Y[,3], ylim = c(0,1), xlab = "longitude", ylab = expression(Y[3]), pch=3, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$longitude[c(22,23,41)], Y[c(22,23,41),3], col = "red", pch = 1, cex = 2)
text(gemas.FRA$longitude[c(23)], Y[c(23),3], labels = c("23"), col = "red", font = 1.5, pos = 4, cex = 1.2)
text(gemas.FRA$longitude[c(22,41)], Y[c(22,41),3], labels = c("22","41"), col = "red", font = 1.5, pos = 1, cex = 1.2)
dev.off()

png("./figures/gemas/lat_comp.png", width = 500, height = 450)
par(mfrow=c(2,2), mar=c(4,5,2,1))
plot(gemas.FRA$latitude, Y[,1], ylim = c(0,1), xlab = "latitude", ylab = "Y", pch=1, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$latitude, Y[,2], col="blue", pch=2)
points(gemas.FRA$latitude, Y[,3], col="red", pch=3)
plot(gemas.FRA$latitude, Y[,1], ylim = c(0,1), xlab = "latitude", ylab = expression(Y[1]), pch=1, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$latitude[207], Y[207,1], col = "red", pch = 1, cex = 2)
text(gemas.FRA$latitude[207], Y[207,1], labels = "207", col = "red", font = 1.5, pos = 1, cex = 1.2)

plot(gemas.FRA$latitude, Y[,2], ylim = c(0,1), xlab = "latitude", ylab = expression(Y[2]), pch=2, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$latitude[22], Y[22,2], col = "red", pch = 1, cex = 2)
text(gemas.FRA$latitude[22], Y[22,2], labels = "22", col = "red", font = 1.5, pos = 1, cex = 1.2)

plot(gemas.FRA$latitude, Y[,3], ylim = c(0,1), xlab = "latitude", ylab = expression(Y[3]), pch=3, cex.lab=1.5, cex.axis = 1.5)
points(gemas.FRA$latitude[41], Y[41,3], col = "red", pch = 1, cex = 2)
text(gemas.FRA$latitude[41], Y[41,3], labels = "41", col = "red", font = 1.5, pos = 1, cex = 1.2)
dev.off()

png("./figures/gemas/MeanTemp_comp.png", width = 500, height = 450)
par(mfrow=c(2,2), mar=c(4,5,2,1))
plot(gemas.FRA$MeanTemp, Y[,1], ylim = c(0,1), xlab = "MeanTemp", ylab = expression(Y[1]), pch=1, cex.lab=1.5, cex.axis = 1.5)
plot(gemas.FRA$MeanTemp, Y[,2], ylim = c(0,1), xlab = "MeanTemp", ylab = expression(Y[2]), pch=2, cex.lab=1.5, cex.axis = 1.5)
plot(gemas.FRA$MeanTemp, Y[,3], ylim = c(0,1), xlab = "MeanTemp", ylab = expression(Y[3]), pch=3, cex.lab=1.5, cex.axis = 1.5)
dev.off()

png("./figures/gemas/ph_comp.png", width = 500, height = 450)
par(mfrow=c(2,2), mar=c(4,5,2,1))
plot(gemas.FRA$ph_cacl2, Y[,1], ylim = c(0,1), xlab = "ph", ylab = expression(Y[1]), pch=1, cex.lab=1.5, cex.axis = 1.5)
plot(gemas.FRA$ph_cacl2, Y[,2], ylim = c(0,1), xlab = "ph", ylab = expression(Y[2]), pch=2, cex.lab=1.5, cex.axis = 1.5)
plot(gemas.FRA$ph_cacl2, Y[,3], ylim = c(0,1), xlab = "ph", ylab = expression(Y[3]), pch=3, cex.lab=1.5, cex.axis = 1.5)
dev.off()

## map of Swedish sampling locations
# lon_rng = range(gemas.SWE$longitude, na.rm = TRUE)
# lat_rng = range(gemas.SWE$latitude, na.rm = TRUE)
# lon_pad = 0.5
# lat_pad = 0.8
# 
# png("./figures/gemas/map_SWE_locations.png", width = 600, height = 700)
# par(mfrow = c(1,1), mar = c(2,2,2,1), mgp = c(2.2,0.8,0))
# maps::map("world", regions = "FRA",
#           xlim = lon_rng + c(-.5, .8),
#           ylim = lat_rng + c(-.5, .8),
#           fill = TRUE, col = "grey95", border = "grey60")
# box()
# points(gemas.SWE$longitude[which(gemas.SWE$latitude>60)], gemas.SWE$latitude[which(gemas.SWE$latitude>60)],
#        pch = 16, col = "darkblue", cex = 0.8)
# points(gemas.SWE$longitude[which(gemas.SWE$latitude<=60)], gemas.SWE$latitude[which(gemas.SWE$latitude<=60)],
#        pch = 17, col = "darkgreen", cex = 0.8)
# points(gemas.SWE$longitude[c(99,100)], gemas.SWE$latitude[c(99,100)],
#        pch = 1, col = "red", cex = 1.5, lwd = 2)
# text(gemas.SWE$longitude[99] + 1, gemas.SWE$latitude[99],
#      labels = "99", col = "red", cex = 0.9)
# text(gemas.SWE$longitude[100] + 1.1, gemas.SWE$latitude[100] ,
#      labels = "100", col = "red", cex = 0.9)
# title(main = "GEMAS sand observations in Sweden", line = 1)
# mtext("Longitude", side = 1, line = 1)
# mtext("Latitude", side = 2, line = 1)
# dev.off()


## 1.2 categorical X ----
gemas.FRA$soilclass = factor(gemas.FRA$soilclass, levels = c("ll", "l", "m","s"))
table(gemas.FRA$soilclass)

require(robCompositions)
pch.vec <-c(3,15:17)
png(paste("./figures/gemas/ternary_diag.png",sep = ""), width = 400, height = 350)
par(mfrow=c(1,1))
ternaryDiag(gemas.FRA[,c("sand", "silt", "clay")],
            col = as.numeric(gemas.FRA$soilclass),
            pch = pch.vec[as.numeric(gemas.FRA$soilclass)])
legend("topleft", c("ll","l","m","s"), col = 1:4, pch = pch.vec, pt.cex = 1.2)
dev.off()

# for (i in c("sand", "silt", "clay")) {
#     png(paste("./figures/gemas/box_", i, ".png",sep = ""), width = 400, height = 350)
#     par(mfrow=c(1,1), pin=c(6,4), mar=c(4,5,1,1))
#     boxplot(gemas.FRA[,i]~gemas.FRA$soilclass,
#             ylab = i, xlab = "soilclass", col=rainbow(3, s=0.35),
#             cex.lab=2, cex.axis=1.5)
#     dev.off()
# }

## 1.3 continuous X ---- 
# "longitude" "latitude" "MeanTemp" "AnnPrec" "cec" "ph_cacl2" "toc"
# for (i in c("MeanTemp", "AnnPrec", "cec", "ph_cacl2", "toc")) {
#     png(paste("./figures/gemas/hist_", i, ".png",sep = ""), width = 300, height = 200)
#     par(pin=c(6,4), mar=c(4,5,1,1))
#     hist(gemas.FRA[,i], main="", xlab=i, cex.lab=1.5, cex.axis=1.2)
#     dev.off()
# }

gemas.FRA$lAnnPrec = log(gemas.FRA$AnnPrec)
gemas.FRA$ltoc = log(gemas.FRA$toc)

# pairs(gemas.FRA[,-c(4,5,6,7,8,11)])

# for (i in c("lAnnPrec", "ltoc")) {
#     png(paste("./figures/gemas/hist_", i, ".png",sep = ""), width = 300, height = 200)
#     par(pin=c(6,4), mar=c(4,5,1,1))
#     hist(gemas.SWE[,i], main="", xlab=i, cex.lab=1.5, cex.axis=1.2)
#     dev.off()
# }


# for (i in c("ph_cacl2", "ltoc")) {
#     png(paste("./figures/gemas/scatter_", i, ".png",sep = ""), width = 300, height = 200)
#     par(pin=c(6,4), mar=c(4,5,1,1))
#     plot(gemas.FRA[,i], gemas.FRA[,"cec"], xlab = i, ylab = "cec",
#          cex.lab=1.5, cex.axis=1.2, pch=19)
#     dev.off()
# }
cor(gemas.FRA[,"ph_cacl2"], gemas.FRA[,"cec"])
cor(gemas.FRA[,"ltoc"], gemas.FRA[,"cec"])
cor(gemas.FRA[,"ltoc"], gemas.FRA[,"ph_cacl2"])
cor(gemas.FRA[,"ltoc"], gemas.FRA[,"MeanTemp"])

## standardize continuous X
X = apply(gemas.FRA[, c("longitude", "latitude", "MeanTemp", "lAnnPrec", "ltoc", "ph_cacl2")], 2,
          function(x){ (x - mean(x))/sd(x)})
## encode categorical X, "ll" is the reference level
X = cbind(model.matrix(~gemas.FRA$soilclass)[,-1], X)
colnames(X) = c("l","m","s", "lon","lat","MeanTemp", "lAnnPrec","ltoc","ph_cacl2")

# pairs(X)
rm(gemas.FRA)

# 2.Fit model ----
source("ESLFns.R")

## 2.1 10-fold CV for MSE and MSPE ----
MSE_LS_cv <- MSPE_LS_cv <- MSE_ESL_cv <- MSPE_ESL_cv <- MSE_FSIM_cv <- MSPE_FSIM_cv <- MSE_SIQR_cv <- MSPE_SIQR_cv <- rep(0, 10)
theta_LS.cv <- theta_ESL.cv <- theta_FSIM.cv <- theta_SIQR.cv <- matrix(0, nrow = 10, ncol = ncol(X))
colnames(theta_LS.cv)= colnames(theta_ESL.cv) = colnames(theta_FSIM.cv) = colnames(theta_SIQR.cv) = colnames(X)
h_LS.cv <- h_ESL.cv <- h_FSIM.cv <- h_SIQR.cv <- lambda.cv <- rep(0,10)
p = dim(X)[2]
n = nrow(X)
d = ncol(Y)

set.seed(1)

id = sample(1:n)
fold_id = split(id, rep(1:10, length.out = n))

for (fold in 1:10) {
    test_id <- fold_id[[fold]]
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
    ESL_cv = esl_est(theta_init = rep(1,p), bw_init = 0.5,
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
    try({FSIM_cv = fsim_est(xdat = X_train, ydat = Y_train, init = rep(1,p))
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
    })
    
    # fit SIQR
    SIQR.fit = index.gamma(y=Y_train, xx=X_train, tau=0.5, gamma0 = rep(1,p), maxiter = 100, crit = 1e-3)
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
    
    ## save cv results
    save.image("./data/GEMAS_sand.Rdata")
}



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
save.image("./data/GEMAS_sand.Rdata")


## 2.3 model with full samples after removing suspectable outliers ----
# remove suspectable outliers 22 23 41 205 207
outlier_id = c(22,23,41,205,207)
X_sub = X[-outlier_id,]
Y_sub = Y[-outlier_id,]
n_sub = nrow(X_sub)

# fit ESIM(LS)
LS_sub = ls_est(c(LS$bw, LS$theta), xdat = X_sub, ydat = Y_sub)
h_LS_sub = LS_sub$bw
theta_LS_sub = SpheNormalize(LS_sub$theta)
LS_mu_sub = t(apply(LS_sub$mu, 1, SpheNormalize))
## MSE(LS)
MSE_LS_sub = mean(acos(rowSums(LS_mu_sub * Y_sub))^2, na.rm = T)

# fit ESIM(ESL)
ESL_sub = esl_est(theta_init = ESL$theta, bw_init = ESL$bw,
                  xdat = X_sub, ydat = Y_sub, delta = 0.2)
h_ESL_sub= ESL_sub$bw
theta_ESL_sub = SpheNormalize(ESL_sub$theta)
ESL_mu_sub = t(apply(ESL_sub$mu, 1, SpheNormalize))
## MSE(ESL)
MSE_ESL_sub = mean(acos(rowSums(ESL_mu_sub * Y_sub))^2, na.rm = T)

# fit FSIM
FSIM_sub = fsim_est(xdat = X_sub, ydat = Y_sub, init = FSIM$theta)
h_FSIM_sub = FSIM_sub$bw
theta_FSIM_sub = SpheNormalize(FSIM_sub$theta)
## MSE(FSIM)
MSE_FSIM_sub = mean(acos(rowSums(FSIM_sub$mu * Y_sub))^2, na.rm = T)

# fit SIQR
SIQR_sub = index.gamma(y=Y_sub, xx=X_sub, tau=0.5, gamma0 = SIQR$theta, maxiter = 100, crit = 1e-3)
theta_SIQR_sub = SpheNormalize(SIQR_sub$theta)
h_SIQR_sub = SIQR_sub$bw
## MSE(SIQR)
SIQR_index_sub = X_sub %*% t(t(theta_SIQR_sub))
SIQR_mu_sub = sapply(1:length(SIQR_index_sub), function(r){
    lprq_mul(SIQR_index_sub[r], SIQR_index_sub[-r], Y_sub[-r,], h=h_SIQR_sub)
})
SIQR_mu_sub = t(apply(SIQR_mu_sub, 2, SpheNormalize))
MSE_SIQR_sub = mean(acos(rowSums(SIQR_mu_sub * Y_sub))^2, na.rm = T)

result_sub = matrix(c(theta_LS_sub, h_LS_sub, MSE_LS_sub,
                      theta_ESL_sub, h_ESL_sub, MSE_ESL_sub,
                      theta_FSIM_sub, h_FSIM_sub, MSE_FSIM_sub,
                      theta_SIQR_sub, h_SIQR_sub, MSE_SIQR_sub), byrow = T, nrow = 4)
colnames(result_sub) = c(colnames(X_sub), "h", "MSE")
rownames(result_sub) = c("LS", "ESL", "FSIM", "SIQR")
result_sub
save.image("./data/GEMAS_sand.Rdata")

## 2.4 Bootstrap for SE of LS and ESL ----
nIter = 500
b_ESL.boot = b_LS.boot = matrix(0, nrow = nIter, ncol = ncol(X))
colnames(b_LS.boot) = colnames(b_ESL.boot) = colnames(X)
h_ESL.boot = h_LS.boot = rep(0, nIter)
b_ESL_sub.boot = b_LS_sub.boot = matrix(0, nrow = nIter, ncol = ncol(X_sub))
colnames(b_LS_sub.boot) = colnames(b_ESL_sub.boot) = colnames(X_sub)
h_ESL_sub.boot = h_LS_sub.boot = rep(0, nIter)

set.seed(100)
for (b in 1:nIter) {
    # weighted bootstrap for LS with full sample
    w = rexp(n, rate = 1)
    w = n * w / sum(w)
    LS_boot = ls_est_weight(c(LS$bw, LS$theta), xdat = X, ydat = Y, weight = w)
    b_LS.boot[b,] = LS_boot$theta; h_LS.boot[b] = LS_boot$bw
    # weighted bootstrap for ESL with full sample
    w = rexp(n, rate = 1)
    w = n * w / sum(w)
    ESL_boot = esl_index_bw(theta_init = ESL$theta, bw_init = ESL$bw,
                            xdat = X, ydat = Y, lambda = ESL$lambda, w = w)
    b_ESL.boot[b,] = ESL_boot$theta; h_ESL.boot[b] = ESL_boot$bw
    # weighted bootstrap for LS with subset
    w = rexp(n_sub, rate = 1)
    w = n_sub * w / sum(w)
    LS_sub_boot = ls_est_weight(c(LS_sub$bw, LS_sub$theta), xdat = X_sub, ydat = Y_sub, weight = w)
    b_LS_sub.boot[b,] = LS_sub_boot$theta; h_LS_sub.boot[b] = LS_sub_boot$bw
    # weighted bootstrap for ESL with subset
    w = rexp(n_sub, rate = 1)
    w = n_sub * w / sum(w)
    ESL_sub_boot = esl_index_bw(theta_init = ESL_sub$theta, bw_init = ESL_sub$bw,
                                xdat = X_sub, ydat = Y_sub, lambda = ESL_sub$lambda, w = w)
    b_ESL_sub.boot[b,] = ESL_sub_boot$theta; h_ESL_sub.boot[b] = ESL_sub_boot$bw
    save.image("./data/GEMAS_sand.Rdata")
}

## 2.5 summary stats ----
# rm(list = ls())
load("./data/GEMAS_sand.Rdata")
round(result, 4)
round(result_sub, 4)
# 
round(mean(MSPE_LS_cv), 4)
round(mean(MSPE_ESL_cv), 4)
round(mean(MSPE_FSIM_cv, na.rm = T), 4)
round(mean(MSPE_SIQR_cv), 4)
# 
round(apply(b_LS.boot, 2, sd),3)
round(apply(b_ESL.boot, 2, sd),3)
round(apply(b_LS_sub.boot, 2, sd),3)
round(apply(b_ESL_sub.boot, 2, sd),3)

# round(sd(h_LS.boot, na.rm = T), 3)
# round(sd(h_ESL.boot[which(h_ESL.boot<1)], na.rm = T), 3) # remove unusual runs
# round(sd(h_FSIM.boot, na.rm = T), 3)
# round(sd(h_SIQR.boot, na.rm = T), 3)

# A = diag(p-1)[3:5,]
# theta_null = c(1,1,1,0,0,0); theta_null = SpheNormalize(theta_null)
# Sigma_LS = cov(b_LS.boot[,-1])
# Sigma_ESL = cov(b_ESL.boot[,-1])
# T_LS = t(A %*% (LS$theta[-1] - theta_null[-1])) %*% solve(A %*% Sigma_LS %*% t(A)) %*% (A %*% (LS$theta[-1] - theta_null[-1]))
# T_ESL = t(A %*% (ESL$theta[-1] - theta_null[-1])) %*% solve(A %*% Sigma_ESL %*% t(A)) %*% (A %*% (ESL$theta[-1] - theta_null[-1]))

# 3. Diagnostic figure ----
# rm(list = ls())
# load("GEMAS_sand.Rdata")
# # rotated residuals for sampling, tangent space at the north pole c(1,0,0)
# e_rot_LS = t(sapply(1:n, function(i){
#     ## crude residuals
#     u = LS_mu_full[i,]
#     v = Y[i,]
#     u = u / sqrt(sum(u^2))
#     v = v / sqrt(sum(v^2))
#     eps = v - sum(v * u) * u
#     ## Rotation 
#     R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
#     
#     return(R %*% t(t(eps))) 
# }))
# 
# e_rot_ESL = t(sapply(1:n, function(i){
#     ## project y at the tangent space of mu
#     u = ESL_mu_full[i,]
#     v = Y[i,]
#     u = u / sqrt(sum(u^2))
#     v = v / sqrt(sum(v^2))
#     eps = v - sum(v * u) * u
#     ## Rotation
#     R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
#     
#     return(R %*% t(t(eps))) 
# }))
# 
# e_rot_FSIM = t(sapply(1:n, function(i){
#     ## crude residuals
#     u = FSIM$mu[i,]
#     v = Y[i,]
#     u = u / sqrt(sum(u^2))
#     v = v / sqrt(sum(v^2))
#     
#     eps = v - sum(v * u) * u
#     ## Rotation
#     R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
#     
#     return(R %*% t(t(eps))) 
# }))
# 
# e_rot_SIQR = t(sapply(1:n, function(i){
#     u = SIQR_mu[i,]
#     v = Y[i,]
#     u = u / sqrt(sum(u^2))
#     v = v / sqrt(sum(v^2))
#     eps = v - sum(v * u) * u
#     ## Rotation
#     R = t(t(u + c(1, 0, 0))) %*% t(u + c(1, 0, 0)) / (1 + u[1]) - diag(d)
#     
#     return(R %*% t(t(eps))) 
# }))
# 
# group_colors <- c("#1b9e77", "#d95f02", "#7570b3")
# ## 3.1 residuals plot for LS ----
# png("../figures/gemas/sand_LS_res.png", width = 250, height = 250)
# par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
# plot(e_rot_LS[,2], e_rot_LS[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
#      xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
#      col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
# dev.off()
# 
# LS_index_full = X %*% LS$theta
# png("../figures/gemas/sand_LS_res2.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(LS_index_full, e_rot_LS[,2], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
# points(LS_index_full[X[,"ll"]==1,], e_rot_LS[X[,"ll"]==1,2], 
#        col = group_colors[1], pch=16, cex=1.5)
# points(LS_index_full[X[,"m"]==1,], e_rot_LS[X[,"m"]==1,2], 
#        col = group_colors[2], pch=17, cex=1.5)
# points(LS_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
#        e_rot_LS[(X[,"m"]==0)&(X[,"ll"]==0),2], 
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# png("../figures/gemas/sand_LS_res3.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(LS_index_full, e_rot_LS[,3], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
# points(LS_index_full[X[,"ll"]==1,], e_rot_LS[X[,"ll"]==1,3],
#        col = group_colors[1], pch=16, cex=1.5)
# points(LS_index_full[X[,"m"]==1,], e_rot_LS[X[,"m"]==1,3], 
#        col = group_colors[2], pch=17, cex=1.5)
# points(LS_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
#        e_rot_LS[(X[,"m"]==0)&(X[,"ll"]==0),3],
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# 
# ## 3.2 residuals plot for ESL ----
# ESL_index_full = X %*% theta_ESL
# 
# png("../figures/gemas/sand_ESL_res.png", width = 250, height = 250)
# par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
# plot(e_rot_ESL[,2], e_rot_ESL[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
#      xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
#      col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
# dev.off()
# 
# png("../figures/gemas/sand_ESL_res2.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(ESL_index_full, e_rot_ESL[,2], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
# points(ESL_index_full[X[,"ll"]==1,], e_rot_ESL[X[,"ll"]==1,2], 
#        col = group_colors[1], pch=16, cex=1.5)
# points(ESL_index_full[X[,"m"]==1,], e_rot_ESL[X[,"m"]==1,2], 
#        col = group_colors[2], pch=17, cex=1.5)
# points(ESL_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
#        e_rot_ESL[(X[,"m"]==0)&(X[,"ll"]==0),2], 
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# png("../figures/gemas/sand_ESL_res3.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(ESL_index_full, e_rot_ESL[,3], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
# points(ESL_index_full[X[,"ll"]==1,], e_rot_ESL[X[,"ll"]==1,3],
#        col = group_colors[1], pch=16, cex=1.5)
# points(ESL_index_full[X[,"m"]==1,], e_rot_ESL[X[,"m"]==1,3],
#        col = group_colors[2], pch=17, cex=1.5)
# points(ESL_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
#        e_rot_ESL[(X[,"m"]==0)&(X[,"ll"]==0),3], 
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# ## 3.3 residuals plot for FSIM ----
# FSIM_index_full = as.vector(X %*% FSIM$theta)
# png("../figures/gemas/sand_FSIM_res.png", width = 250, height = 250)
# par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
# plot(e_rot_FSIM[,2], e_rot_FSIM[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
#      xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
#      col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
# dev.off()
# 
# png("../figures/gemas/sand_FSIM_res2.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(FSIM_index_full, e_rot_FSIM[,2], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
# points(FSIM_index_full[X[,"ll"]==1], e_rot_FSIM[X[,"ll"]==1,2], 
#        col = group_colors[1], pch=16, cex=1.5)
# points(FSIM_index_full[X[,"m"]==1], e_rot_FSIM[X[,"m"]==1,2], 
#        col = group_colors[2], pch=17, cex=1.5)
# points(FSIM_index_full[(X[,"m"]==0)&(X[,"ll"]==0)], 
#        e_rot_FSIM[(X[,"m"]==0)&(X[,"ll"]==0),2], 
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# png("../figures/gemas/sand_FSIM_res3.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(FSIM_index_full, e_rot_FSIM[,3], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
# points(FSIM_index_full[X[,"ll"]==1], e_rot_FSIM[X[,"ll"]==1,3], 
#        col = group_colors[1], pch=16, cex=1.5)
# points(FSIM_index_full[X[,"m"]==1], e_rot_FSIM[X[,"m"]==1,3], 
#        col = group_colors[2], pch=17, cex=1.5)
# points(FSIM_index_full[(X[,"m"]==0)&(X[,"ll"]==0)], 
#        e_rot_FSIM[(X[,"m"]==0)&(X[,"ll"]==0),3], 
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# ## 3.4 residuals plot for SIQR ----
# SIQR_index_full = X %*% SIQR$theta
# 
# png("../figures/gemas/sand_SIQR_res.png", width = 250, height = 250)
# par(mar=c(3.,3.,.5,0.5), mgp=c(1.6,0.25,0))
# plot(e_rot_SIQR[,2], e_rot_SIQR[,3], xlim = c(-.5, .5), ylim = c(-.5, .5),
#      xlab = expression(epsilon[2]), ylab = expression(epsilon[3]), cex.lab=2, asp=1,
#      col = rgb(0.2, 0.4, 0.8, 0.6), pch = 16, cex=1.25)
# dev.off()
# 
# png("../figures/gemas/sand_SIQR_res2.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(SIQR_index_full, e_rot_SIQR[,2], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[2]), side = 2, line = 1, cex = 2)
# points(SIQR_index_full[X[,"ll"]==1,], e_rot_SIQR[X[,"ll"]==1,2],
#        col = group_colors[1], pch=16, cex=1.5)
# points(SIQR_index_full[X[,"m"]==1,], e_rot_SIQR[X[,"m"]==1,2],
#        col = group_colors[2], pch=17, cex=1.5)
# points(SIQR_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
#        e_rot_SIQR[(X[,"m"]==0)&(X[,"ll"]==0),2], 
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
# 
# png("../figures/gemas/sand_SIQR_res3.png", width = 375, height = 250)
# par(mar=c(3,3,.5,0.5), mgp=c(1.6,0.25,0))
# plot(SIQR_index_full, e_rot_SIQR[,3], ylim = c(-.5,.5), 
#      xlab = "", ylab = "", cex.lab = 1.5, type = "n")
# mtext(expression(hat(beta)^T * X), side = 1, line = 2, cex = 1.3)
# mtext(expression(epsilon[3]), side = 2, line = 1, cex = 2)
# points(SIQR_index_full[X[,"ll"]==1,], e_rot_SIQR[X[,"ll"]==1,3], 
#        col = group_colors[1], pch=16, cex=1.5)
# points(SIQR_index_full[X[,"m"]==1,], e_rot_SIQR[X[,"m"]==1,3], 
#        col = group_colors[2], pch=17, cex=1.5)
# points(SIQR_index_full[(X[,"m"]==0)&(X[,"ll"]==0),], 
#        e_rot_SIQR[(X[,"m"]==0)&(X[,"ll"]==0),3],
#        col = group_colors[3], pch=18, cex=1.5)
# legend("bottomright", legend = c("ll","m","l"), pch=16:18, col=group_colors, cex=1.2, ncol = 3, pt.cex = 1.5)
# dev.off()
