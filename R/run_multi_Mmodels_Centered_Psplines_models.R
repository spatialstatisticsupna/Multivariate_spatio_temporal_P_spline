################################################################################
######                    Multivariate spatio-temporal                    ######
######                        M-models / P-splines                        ######
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
m.model<- "fe"  # "fe": FE; "re": RE 

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

## D.spat = I_{J} (x) 1_{t} (x) I_{n}
D.spat<- kronecker(diag(1, nrow = k), kronecker(matrix(1, nrow = t, ncol=1), diag(1, nrow = n) ))

## D.temp = I_{J} (x) Bt (x) 1_{n}
D.temp<- kronecker(diag(1, nrow = k), kronecker(Bt, matrix(1, nrow=n, ncol=1)))

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
id.idx<- c(rep(NA,k), 1:(n*k), rep(NA, kt*k + n*t*k))

## id.idy
id.idy<- c(rep(NA,k + n*k), 1:(kt*k), rep(NA, n*t*k))

## id.inter and ## id.idxy 
for(i in 1:k){
  ## id.inter
  aux<- rep(NA, k*(1+n+kt+n*t))
  aux[i]<-1
  eval(parse(text = paste0("id.inter.",i,"<-aux") ))
  rm(aux)
  
  ## id.idxy
  aux<- rep(NA, k*(1+n+kt+n*t))
  aux[(k*(1+n+kt)+(i-1)*(n*t))+1:(n*t)]<- 1:(n*t)
  eval(parse(text = paste0("id.idxy.",i,"<-aux") ))
  rm(aux)
}

## Data
Data<- list(OBS = d.data_frame$OBS, EXP = d.data_frame$EXP, idx=id.idx, idy=id.idy)
for (i in 1:k){
  eval(parse(text = paste0("Data$intercept.",i,"<-id.inter.",i) ))
  eval(parse(text = paste0("Data$idxy.",i,"<- id.idxy.",i) ))
}

## rm
rm(list = c("id.idx","id.idy",paste0(rep(c("id.inter.","id.idxy."),each=k),1:k)))

################################################################################
## Adjacency and precision matrices                                           ##
################################################################################

################################
## W: adjacency matrix (spatial)
################################
adj<-spdep::poly2nb(carto)
W <- as(nb2mat(adj, style = "B"), "Matrix")

############
## W: adjacency matrix (temporal)
############
Pt<- crossprod(diff(diag(kt),differences=ord.temp))
Wt<- inla.as.sparse(Diagonal(n=nrow(Pt),diag(Pt))- Pt)


############
## Precision matrices for spatio and temporal interaction
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
## spatial
A.constr.s<- kronecker(diag(1,k), matrix(1, nrow=1, ncol=n))
e.s<- rep(0,dim(A.constr.s)[1])

## temporal
A.constr.t<- kronecker(diag(1,k), matrix(apply(Bt,2,sum), nrow=1, ncol=kt)) 
e.t<- rep(0,dim(A.constr.t)[1])

## spatio-temporal
A.constr.st.1<- matrix(1, nrow=1, ncol=n*t)
e.st.1= rep(0,dim(A.constr.st.1)[1])

aux1<- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
aux2<- rbind(aux1, kronecker(matrix(1:t,nrow=1, ncol=t), matrix(c(1, rep(0,n-1)),nrow=1, ncol=n)))
if(ord.temp==1){
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
source("functions/posterior_lincombs_temp.R")

## spalial M-models
source(paste0("functions/Mmodel_icar_",m.model,".R"))
eval(parse(text = paste0("funct.s<-inla.rgeneric.Mmodel.icar.",m.model) ))
model.s <- inla.rgeneric.define(funct.s, debug = TRUE, k = k, W = W)


## temporal P-spline
source("functions/mult_temp_Psplines.R")
model.t <- inla.rgeneric.define(mult.temp.Psplines, debug = TRUE, k = k, W = Wt)

## Define appropriate hyperprior distributions (spatio-temporal)
sdunif="expression:
logdens=-log_precision/2;
return(logdens)"

## formula
formula.inla<- OBS ~ -1 +  intercept.1 + intercept.2 + intercept.3 + intercept.4 +
  f(idx, model = model.s, constr=FALSE, extraconstr=list(A=A.constr.s, e=e.s )) +
  f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A.constr.t, e=e.t)) +
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
                                                   D.spat, D.temp,
                                                   D.spat.temp.1, D.spat.temp.2, D.spat.temp.3, D.spat.temp.4), 
                                          link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  lincomb=all.lc,
                  control.inla = list(strategy=strategy),
                  verbose = FALSE)


## save
save(model, file= paste0("resul/resul_Centered_Psplines_M_",m.model,"_pt",ord.temp,"_type",type,"_",gsub("\\.", "_", strategy),".RData") )

################################################################################
################################################################################