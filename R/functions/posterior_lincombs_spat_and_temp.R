################################################################################
## linear combinations: Multivariate P-splines models                         ##
################################################################################
#############
## with temporal correlations
#############
## posteriori distribution -  spatial spline
lc.spat<- inla.make.lincombs(idx = kronecker(diag(k),Bs))
names(lc.spat)<- paste("spatial.",as.character(seq(1,k*n,1)),sep="")

## posteriori distribution - temporal spline
lc.temp<- inla.make.lincombs(idy = kronecker(diag(k),Bt))
names(lc.temp)<- paste("temporal.",as.character(seq(1,k*t,1)),sep="")

all.lc <- c(lc.spat,lc.temp)

#############
## without temporal correlations
#############
## posteriori distribution - temporal spline (without temporal correlations)
lc.temp.1<- inla.make.lincombs(idy.1 = Bt)
names(lc.temp.1)<- paste("temporal.1.",as.character(seq(1,t,1)),sep="")
lc.temp.2<- inla.make.lincombs(idy.2 = Bt)
names(lc.temp.2)<- paste("temporal.2.",as.character(seq(1,t,1)),sep="")
lc.temp.3<- inla.make.lincombs(idy.3 = Bt)
names(lc.temp.3)<- paste("temporal.3.",as.character(seq(1,t,1)),sep="")
lc.temp.4<- inla.make.lincombs(idy.4 = Bt)
names(lc.temp.4)<- paste("temporal.4.",as.character(seq(1,t,1)),sep="")


all.lc.ind <- c(lc.spat,lc.temp.1,lc.temp.2,lc.temp.3,lc.temp.4)

################################################################################
## linear combinations: Univariate P-splines models                           ##
################################################################################
## posteriori distribution - espatial spline
lc.spat.univ<- inla.make.lincombs(idx = Bs)
names(lc.spat.univ)<- paste("spatial.",as.character(seq(1,n,1)),sep="")

## posteriori distribution - temporal spline
lc.temp.univ<- inla.make.lincombs(idy = Bt)
names(lc.temp.univ)<- paste("temporal.",as.character(seq(1,t,1)),sep="")

all.lc.univ <- c(lc.spat.univ,lc.temp.univ)

################################################################################
################################################################################