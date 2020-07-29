################################################################################
######            Univariate spatio-temporal P-splines models             ######
######                      with centered coefficients                    ######
################################################################################
rm(list=ls())

## libraries
library(INLA); library(spdep); library(spData); library(splines)

### save results
if(!file.exists("resul")) {dir.create("resul")}

################################################################################
## Selections
################################################################################
############
## Selected prior distributions for the coefficients of the
## spatial and temporal P-splines (RW1 or RW2)
############
ord.spat<- 1  # 1: RW1; 2: RW2 
ord.temp<- 1  # 1: RW1; 2: RW2

############
## strategy
############
strategy <- "simplified.laplace" # see 'control.inla' (?control.inla)

################################################################################
## data loading and organization                                              ##
################################################################################
## Load data and Maharashtra SpatialPolygonsDataFrame
load("data_Psplines.RData")


## Number of areas and number of time periods
n<- length(unique(data$ID_area))
t<- length(unique(data$ID_year))

## selected crimes
crimes<- c("rape", "assault", "cruelty", "kidnapping")
e_crimes<- paste0("e_",crimes)

## Parameter of the model (number of crimes)
k <- length(crimes)


## data.frame ID_disease
data<- data[order(data$ID_year, data$ID_area),]

################################################################################
## Bs and Bt                                                                  ##
################################################################################
## Order of the splines
bdeg <- 3

## function for Bs= (B2 kronecker 1) (.) (1 kronecker B1)
Rten <- function(X1,X2){
  one1 <- matrix(1,1,ncol(X1))
  one2 <- matrix(1,1,ncol(X2))
  kronecker(X1,one2)*kronecker(one1,X2)
}

## Rescale the longitude and latitude covariate into the [0,1] interval
x1 <- coordinates(carto)[,1]
x1 <- (x1-min(x1))/(max(x1)-min(x1))
ndx1 <- ceiling(min(length(unique(x1))/4, 40))
dis1 <- (max(x1)-min(x1))/ndx1
x1l <- min(x1)-dis1*0.05
x1r <- max(x1)+dis1*0.05
dx1 <- (x1r-x1l)/ndx1
knots1 <- seq(x1l-bdeg*dx1, x1r+bdeg*dx1, by=dx1)
B1 <- spline.des(knots1,x1,bdeg+1, 0*x1)$design
k1 <- dim(B1)[2]


x2 <- coordinates(carto)[,2]
x2 <- (x2-min(x2))/(max(x2)-min(x2))
ndx2 <- ceiling(min(length(unique(x2))/4, 40))
dis2 <- (max(x2)-min(x2))/ndx2
x2l <- min(x2)-dis2*0.05
x2r <- max(x2)+dis2*0.05
dx2 <- (x2r-x2l)/ndx2
knots2 <- seq(x2l-bdeg*dx2, x2r+bdeg*dx2, by=dx2)
B2 <- spline.des(knots2,x2,bdeg+1, 0*x2)$design
k2 <- dim(B2)[2]
Bs <- Rten(B2,B1)
ks<- dim(Bs)[2]


## Rescale the longitude and latitude covariate into the [0,1] interval ##
xt <- unique(data$year)
xt <- (xt-min(xt))/(max(xt)-min(xt))
ndxt <- ceiling(min(length(unique(xt))/4, 40))
dis.t <- (max(xt)-min(xt))/ndxt
xtl <- min(xt)- dis.t*0.05
xtr <- max(xt)+ dis.t*0.05
dxt <- (xtr-xtl)/ndxt
knots.t <- seq(xtl- bdeg*dxt, xtr+bdeg*dxt, by=dxt)
Bt <- spline.des(knots.t,xt,bdeg+1, 0*xt)$design
kt <- dim(Bt)[2]

################################################################################
## Design matrices for random effects                                         ##
################################################################################
## D.inter 
D.inter<- matrix(1, nrow = n*t, ncol=1)

## D.spat
D.spat<- kronecker(matrix(1, nrow = t, ncol=1), Bs)

## D.temp
D.temp<- kronecker(Bt, matrix(1, nrow=n, ncol=1))

## D.spat.temp
D.spat.temp<- diag(1,t*n)

################################################################################
## Adjacency and precision matrices                                           ##
################################################################################
############
## Spatial structure matrices
############
P1<- crossprod(diff(diag(k1),differences=ord.spat))
P2<- crossprod(diff(diag(k2),differences=ord.spat))

## Spatial structure matrices
R1 <- kronecker(Diagonal(nrow(P2),1), P1)
R2 <- kronecker(P2, Diagonal(nrow(P1),1))
Cmat.s <- list(inla.as.sparse(R1),inla.as.sparse(R2))

############
## Temporal structure matrix
############
Pt<- crossprod(diff(diag(kt),differences=ord.temp))
Cmat.t <- list(inla.as.sparse(Pt))

############
## Precision matrices for spatio-temporal interaction
############
spdep::nb2INLA("carto_nb.graph", spdep::poly2nb(carto))
g<- INLA::inla.read.graph("carto_nb.graph")
Qs<- matrix(0, g$n, g$n)
for(i in 1:g$n){
  Qs[i,i]=g$nnbs[[i]]
  Qs[i,g$nbs[[i]]]=-1
}

Qt<- crossprod(diff(diag(t), differences=ord.temp))

## Precision matrices for spatio-temporal interaction (generic0) = prec todos los crimes
R.1<- diag(n*t)
r.def.1<- 0
R.2<- kronecker(Qt,diag(n))
r.def.2<- ord.temp*n
R.3<- kronecker(diag(t),Qs)
r.def.3<- t
R.4<- kronecker(Qt, Qs)
r.def.4<- ord.temp*n+t-1

################################################################################
## Constraints                                                                ##
################################################################################
## spatio-temporal
A.constr.st.1<- matrix(1, nrow=1, ncol=n*t)
e.st.1= rep(0,dim(A.constr.st.1)[1])

aux1<- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
aux2<- rbind(aux1, kronecker(matrix(1:t,nrow=1, ncol=t), matrix(c(1, rep(0,n-1)),nrow=1, ncol=n)))
if(ord.spat==1 & ord.temp==1){
  A.constr.st.2<- aux1
} else{
  A.constr.st.2<- aux2
}
e.st.2= rep(0,dim(A.constr.st.2)[1])

A.constr.st.3<- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
e.st.3=rep(0,dim(A.constr.st.3)[1])

A.constr.st.4<- rbind(A.constr.st.2, A.constr.st.3)
e.st.4=rep(0,dim(A.constr.st.4)[1])

################################################################################
## Define appropriate hyperprior distributions                                ##
################################################################################
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

################################################################################
## Fit model for each crime                                                   ##
################################################################################
## To compute patterns (as linear combinations of the log-risks)
source("functions/posterior_lincombs_spat_and_temp.R")

for(mm in 1:k){
  ######################################
  ## Data by crime
  ######################################
  d.data_frame<- data[,c(crimes[mm], e_crimes[mm], "ID_area", "ID_year" )]
  colnames(d.data_frame)<- c("OBS", "EXP", "ID_area", "ID_year")
  
  ## For models with spatio-temporal interactions random effects
  id.inter<-c(1, rep(NA, ks + kt + n*t))
  id.idx<- c(rep(NA,1), 1:ks, rep(NA, kt + n*t))
  id.idy<- c(rep(NA, 1 + ks), 1:kt, rep(NA, n*t))
  id.idxy<- c(rep(NA, 1 + ks + kt), 1:(n*t))
  
  Data<- list(OBS = d.data_frame$OBS, EXP = d.data_frame$EXP,
              intercept=id.inter, idx=id.idx, idy=id.idy, idxy=id.idxy)
  rm(list = c("id.inter","id.idx","id.idy","id.idxy"))
  
  ######################################
  ## Fit models                       ##
  ######################################
  resulta<- list()
  for(type in 1:4){
    eval(parse(text = paste0("R.st<-R.",type)))
    eval(parse(text = paste0("r.def.st<-r.def.",type) ))
    eval(parse(text = paste0("A.constr.st<- A.constr.st.",type) ))
    eval(parse(text = paste0("e.st<- e.st.",type) ))
    
    ## formula
    formula.inla<- OBS ~ -1 + intercept +
      f(idx, model = "generic3", Cmatrix=Cmat.s, constr=TRUE, hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)) ) +
      f(idy, model = "generic3", Cmatrix=Cmat.t, constr=TRUE, hyper=list(prec1=list(prior=sdunif))) +
      f(idxy, model="generic0", Cmatrix=R.st, rankdef=r.def.st,
        constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st), hyper=list(prec=list(prior=sdunif)) )
    
    ## model
    resulta[[type]]<- inla(formula.inla, data = Data, E = EXP, family = "poisson",
                      control.predictor= list(compute = TRUE, A= cBind(D.inter, D.spat, D.temp,D.spat.temp), link=1, cdf=c(log(1))),
                      control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                      lincomb=all.lc.univ,
                      control.inla = list(strategy=strategy),
                      verbose = FALSE)
    
    ## rm
    rm(list = c("R.st", "r.def.st","A.constr.st","e.st","formula.inla"))
  }
  names(resulta)<- c("univ.t1", "univ.t2", "univ.t3", "univ.t4")
  
  
  ## save
  save(resulta, file = paste0("resul/resul_univ_Psplines_ps",ord.spat,"_pt",ord.temp,"_",crimes[[mm]],"_",gsub("\\.", "_", strategy),".RData"))
  
  ## rm
  rm(list = c("d.data_frame","Data"))
}
################################################################################
################################################################################