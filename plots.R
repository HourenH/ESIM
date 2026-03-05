# Figures for ESIM

# Plots for asymptotic relative efficiency and SGES ----
rm(list=ls())

ARE = NULL
SGES = NULL

delta = seq(0.001,0.999,0.001)
for (d in c(3, 4, 5, 8, 10)) {
    cdelta = 2/(1 - (1 - delta)^(2/(d-1))) - 2 
    ARE = cbind(ARE, cdelta^((d+1)/2) * (cdelta+4)^((d+1)/2) * (cdelta + 2)^(-(d+1)))
    SGES = cbind(SGES, cdelta^(-(d-1)/4) * (cdelta + 4)^((d+1)/4))
}

png(paste("./figures/ARE.png", sep=""), width = 375, height = 275)
par(mfrow=c(1,1), pin=c(6,4), mar=c(4,5,1,1))
plot(delta, (ARE[,1]),
     type = "n", xlab = expression(delta),
     ylab = expression(plain("ARE")(hat(theta)[plain("ESL")])),
     ylim=c(0, 1.3), cex.lab=1.5, cex.axis=1.2)
for (i in 1:5) {
    lines(delta, ARE[,i], col=i, lty=i, lwd = 2.5)
}
legend(0.25, 1.3, legend = c("d=3", "d=4","d=5", "d=8", "d=10"), col = 1:5, lty = 1:5, cex = 0.95, ncol=3, lwd=2.5)
dev.off()

png(paste("./figures/SGES.png", sep=""), width = 375, height = 275)
par(mfrow=c(1,1), pin=c(6,4), mar=c(4,5,1,1))
plot(delta, SGES[,1],
     type = "n", xlab = expression(delta),
     ylab = expression(paste(textstyle("Q"), "(", delta, ")")),
     ylim=c(0, 30), cex.lab=1.5, cex.axis=1.2)
for (i in 1:5) {
    lines(delta, SGES[,i], col=i, lty=i, lwd = 2.5)
}
legend(0.1, 30, legend = c("d=3", "d=4","d=5", "d=8", "d=10"), col = 1:5, lty = 1:5, cex = 0.95, ncol=3, lwd=2.5)
dev.off()

# Illustrative plots for concentration ----
rm(list = ls())
library(rgl)
library(Directional)

k = 50 # k = 5, 20, 50
set.seed(k)
mu = c(0,0.5,1)
mu = mu / sqrt(sum(mu^2))
Y = rvmf(500, mu, k)

options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
points3d(Y, alpha=0.5, col="black", size=3)
points3d(mu, alpha=0.8, col="red", size=8)
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000 ,   1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60)

## ESAG
set.seed(5)

Y = resag(500, mu*5, gam = c(1,0.5))

options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
points3d(Y, alpha=0.5, col="black", size=3)
points3d(mu, alpha=0.8, col="red", size=8)
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000, 1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60)


# Plots for shape of mean functions ----
rm(list = ls())

k = 20
I = seq(-3,3,by=0.02) # index value 
tt = pnorm(I)

mu1 = cbind(2*I/(I^2+2), -2*I/(I^2+2), (I^2-2)/(I^2+2))
mu1 = t(apply(mu1, 1, function(x){x/sqrt(sum(x^2))}))

mu2 = cbind(sqrt(1-tt^2)*cos(pi*tt),
            sqrt(1-tt^2)*sin(pi*tt),
            tt)
mu2 = t(apply(mu2, 1, function(x){x/sqrt(sum(x^2))}))

mu3 = cbind(sin(2*pi*I),
            cos(pi*I),
            2*I/sqrt(I^2+2))
mu3 = t(apply(mu3, 1, function(x){x/sqrt(sum(x^2))}))

Y1 <- Y2 <- Y3 <- matrix(0, nrow = length(tt), ncol = 3)
for (j in 1:length(tt)) {
    Y1[j,] = rvmf(1, mu1[j,], k)
    Y2[j,] = rvmf(1, mu2[j,], k)
    Y3[j,] = rvmf(1, mu3[j,], k)
}

options(rgl.useNULL = TRUE, rgl.printRglwidget = TRUE)
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
lines3d(mu1[order(I),1], mu1[order(I),2], mu1[order(I),3], lwd=3, col="red")
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000 ,   1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60)

clear3d()
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
lines3d(mu2[order(I),1], mu2[order(I),2], mu2[order(I),3], lwd=3, col="red")
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000 ,   1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60)

clear3d()
open3d()
bg3d("white")
par3d(windowRect = c(0,0,300,300))
spheres3d(0,0,0,radius=1.00,color="grey90",alpha=.7, lit=F)
lines3d(mu3[order(I),1], mu3[order(I),2], mu3[order(I),3], lwd=3, col="red")
axes3d()
viewMatrix = matrix(c(-0.999993324, 0.003601942, -0.0006337985, 0,
                      -0.003027634, -0.717994153, 0.6960424185,   0,
                      0.002051935,  0.696039855 , 0.7180004120 ,  0,
                      0.0000000,  0.00000000,  0.00000000 ,   1), 4, byrow = T)
view3d(zoom = 0.9, userMatrix = viewMatrix, fov = 60)
