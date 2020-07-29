################################################################################
## linear combinations: Multivariate P-splines models                         ##
################################################################################
## Distribucion a posteriori del spline temporal
lc.temp<- inla.make.lincombs(idy = kronecker(diag(k),Bt))
names(lc.temp)<- paste("temporal.",as.character(seq(1,k*t,1)),sep="")


all.lc <- c(lc.temp)
################################################################################
################################################################################
