################################################################################
######            Multivariate spatio-temporal P-splines models           ######
######                   (without temporal correlations)                  ######
################################################################################
rm(list=ls())

## libraries
library(INLA); library(spdep); library(spData); library(splines); library(blockmatrix)

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
## selected type of interaction
############
type<- 2 # 1: Type I; 2: Type II; 3: Type III; 4: Type IV 

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
d.data_frame<- NULL
for(i in 1:t){
  for(j in 1:length(crimes)){
    d1<- data[data$ID_year==i, c(crimes[j], e_crimes[j], "ID_area", "ID_year" )]
    d1$ID_disease<- j
    colnames(d1)<- c("OBS", "EXP", "ID_area", "ID_year", "ID_disease")
    d.data_frame<- rbind(d.data_frame, d1)
    rm(d1)
  }
}

## orden: disease - year - area
d.data_frame<- d.data_frame[order(d.data_frame$ID_disease,d.data_frame$ID_year,d.data_frame$ID_area),]

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
## D.inter.j
D.inter<- kronecker(diag(1, nrow = k), matrix(1, nrow = n*t, ncol=1))
eval(parse(text = paste0("D.inter.",1:k,"<- matrix(D.inter[,",1:k,"], ncol=1)")))

## D.spat
D.spat<- kronecker(diag(1, nrow = k), kronecker(matrix(1, nrow = t, ncol=1), Bs))

## D.temp
D.temp<- kronecker(diag(1, nrow = k), kronecker(Bt, matrix(1, nrow=n, ncol=1)))
for(i in 1:k){
  eval(parse(text = paste0("D.temp.",i,"<- D.temp[,(",i,"-1)*kt + 1:kt]")))
}

## D.spat.temp.j
for(i in 1:k){
  aux<- matrix(NA, t*n*k, t*n)
  aux.block<- suppressWarnings(as.blockmatrix(aux, nrowe=t*n, ncole=t*n, nrow=k, ncol=1))
  aux.block[[i]]<- diag(1, t*n)
  aux<- as.matrix(aux.block)
  eval(parse(text = paste0("D.spat.temp.",i,"<-aux") ))
  rm(list = c("aux", "aux.block"))
}

################################################################################
## Data for INLA                                                              ##
################################################################################
## id.idx
id.idx<- c(rep(NA,k), 1:(ks*k), rep(NA, kt*k + n*t*k))

## id.inter, id.idxy, and id.idy 
for(i in 1:k){
  ## id.inter
  aux<- rep(NA, k*(1+ks+kt+n*t))
  aux[i]<-1
  eval(parse(text = paste0("id.inter.",i,"<-aux") ))
  rm(aux)
  
  ## id.idxy
  aux<- rep(NA, k*(1+ks+kt+n*t))
  aux[(k*(1+ks+kt)+(i-1)*(n*t))+1:(n*t)]<- 1:(n*t)
  eval(parse(text = paste0("id.idxy.",i,"<-aux") ))
  rm(aux)
  
  ## id.idy
  aux<- rep(NA,k + ks*k + kt*k + n*t*k)
  aux[(k*(1+ks)+(i-1)*kt) + 1:kt]<- 1:kt
  eval(parse(text = paste0("id.idy.",i,"<-aux") ))
}

## Data
Data<- list(OBS = d.data_frame$OBS, EXP = d.data_frame$EXP, idx=id.idx)
for (i in 1:k){
  eval(parse(text = paste0("Data$intercept.",i,"<-id.inter.",i) ))
  eval(parse(text = paste0("Data$idy.",i,"<- id.idy.",i) ))
  eval(parse(text = paste0("Data$idxy.",i,"<- id.idxy.",i) ))
}


## rm
rm(list = c("id.idx", paste0(rep(c("id.inter.","id.idy.","id.idxy."),each=k),1:k)))

################################################################################
## Adjacency and precision matrices                                           ##
################################################################################
############
## W: adjacency matrix (spatial)
############
P1<- crossprod(diff(diag(k1),differences=ord.spat))
P2<- crossprod(diff(diag(k2),differences=ord.spat))

## Ps
Ps<- kronecker(P2, Diagonal(nrow(P1),1)) + kronecker(Diagonal(nrow(P2),1), P1)
W<- inla.as.sparse(Diagonal(n=nrow(Ps),diag(Ps))- Ps)

## I_{k2} (x) P1 and P2 (x) I_{k1}
I.P1<- kronecker(Diagonal(nrow(P2),1), P1)
W1<- inla.as.sparse(Diagonal(n=nrow(I.P1),diag(I.P1))- I.P1)
P2.I<- kronecker(P2, Diagonal(nrow(P1),1))
W2<- inla.as.sparse(Diagonal(n=nrow(P2.I),diag(P2.I))- P2.I)

############
## W: adjacency matrix (temporal)
############
Pt<- crossprod(diff(diag(kt),differences=ord.temp))
Wt<- inla.as.sparse(Diagonal(n=nrow(Pt),diag(Pt))- Pt)

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

############
## Precision matrices for spatio-temporal interaction (generic0) = prec todos los crimes
############
R.1<- diag(n*t)
r.def.1<- 0
R.2<- kronecker(Qt,diag(n))
r.def.2<- ord.temp*n
R.3<- kronecker(diag(t),Qs)
r.def.3<- t
R.4<- kronecker(Qt, Qs)
r.def.4<- ord.temp*n+t-1

eval(parse(text = paste0("R.st<-R.",type)))
eval(parse(text = paste0("r.def.st<-r.def.",type) ))

################################################################################
## Constraints                                                                ##
################################################################################
#############
## Constraints
#############
## spatial
A.constr.s<- kronecker(diag(1,k), matrix(apply(Bs,2,sum), nrow=1, ncol=ks))
e.s<- rep(0,dim(A.constr.s)[1])

## temporal
A.constr.t<- kronecker(diag(1,k), matrix(apply(Bt,2,sum), nrow=1, ncol=kt)) 
e.t<- rep(0,dim(A.constr.t)[1])

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

eval(parse(text = paste0("A.constr.st<- A.constr.st.",type) ))
eval(parse(text = paste0("e.st<- e.st.",type) ))


################################################################################
## Fit model                                                                  ##
################################################################################
## To compute patterns (as linear combinations of the log-risks)
source("functions/posterior_lincombs_spat_and_temp.R")

## spalial P-spline
source("functions/mult_spat_Psplines.R")
model.s <- inla.rgeneric.define(mult.spat.Psplines, debug = TRUE, k = k, W = W, W1 = W1, W2 = W2)


## Define appropriate hyperprior distributions (temporal and spatio-temporal)
sdunif="expression:
logdens=-log_precision/2;
return(logdens)"

## formula
formula.inla<- OBS ~ -1 +  intercept.1 + intercept.2 + intercept.3 + intercept.4 +
  f(idx, model = model.s, constr=FALSE, extraconstr=list(A=A.constr.s, e=e.s )) +
  f(idy.1, model = "generic3", Cmatrix=Cmat.t, constr=TRUE, hyper=list(prec1=list(prior=sdunif))) +
  f(idy.2, model = "generic3", Cmatrix=Cmat.t, constr=TRUE, hyper=list(prec1=list(prior=sdunif))) +
  f(idy.3, model = "generic3", Cmatrix=Cmat.t, constr=TRUE, hyper=list(prec1=list(prior=sdunif))) +
  f(idy.4, model = "generic3", Cmatrix=Cmat.t, constr=TRUE, hyper=list(prec1=list(prior=sdunif))) +
  f(idxy.1, model="generic0", Cmatrix=R.st, rankdef=r.def.st,
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st), hyper=list(prec=list(prior=sdunif)) ) +
  f(idxy.2, model="generic0", Cmatrix=R.st, rankdef=r.def.st,
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st), hyper=list(prec=list(prior=sdunif)) ) +
  f(idxy.3, model="generic0", Cmatrix=R.st, rankdef=r.def.st,
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st), hyper=list(prec=list(prior=sdunif)) ) +
  f(idxy.4, model="generic0", Cmatrix=R.st, rankdef=r.def.st,
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st), hyper=list(prec=list(prior=sdunif)) )


## model
model<- inla(formula.inla, data = Data, E = EXP, family = "poisson",
                  control.predictor= list(compute = TRUE, 
                                          A= cBind(D.inter.1,D.inter.2,D.inter.3,D.inter.4, 
                                                   D.spat, D.temp.1,D.temp.2,D.temp.3,D.temp.4,
                                                   D.spat.temp.1, D.spat.temp.2, D.spat.temp.3, D.spat.temp.4),
                                          link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  lincomb=all.lc.ind,
                  control.inla = list(strategy=strategy),
                  verbose = FALSE)

## save
save(model, file= paste0("resul/resul_Centered_Psplines_ps",ord.spat,"_pt",ord.temp,"_type",type,"_Indep_",gsub("\\.", "_", strategy),".RData"))

################################################################################
################################################################################