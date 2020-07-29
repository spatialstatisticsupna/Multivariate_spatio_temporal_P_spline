################################################################################
######                Multivariate spatio-temporal M-models               ######
################################################################################
rm(list=ls())

## libraries
library(INLA); library(spdep); library(spData)

### save results
if(!file.exists("resul")) {dir.create("resul")}

################################################################################
## Selections
################################################################################
effects<- "fe"  # "fe": fixed effects M-models; "re": random effects M-models

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
# alpha.min <- 0
# alpha.max <- 1

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

## Intercept for each crime
intercepts<- paste0("I",1:k)
d.data_frame[intercepts]<- NA
for(i in 1:k){
  d.data_frame[d.data_frame$ID_disease==i, intercepts[i]]<- 1
}

## Define idx (different ID_area for different crimes)
d.data_frame$idx<- (d.data_frame$ID_disease-1)*n + d.data_frame$ID_area

## Define idy  (different ID_year for different crimes)
d.data_frame$idy<- (d.data_frame$ID_disease-1) *t +  d.data_frame$ID_year

## Define idxy.j (different ID_area_year for different crimes)
# d.data_frame$idxy<- d.data_frame$ID_area + (d.data_frame$ID_year-1)*n
# ## idxy.
idxy.<- paste0("idxy.", 1:k)
d.data_frame[idxy.]<- NA
for(i in 1:k){
  d.data_frame[d.data_frame$ID_disease==i, idxy.[i]]<- (d.data_frame[d.data_frame$ID_disease==i, c("ID_year")]-1)*n+d.data_frame[d.data_frame$ID_disease==i, c("ID_area")]
}

################################################################################
## Data for INLA                                                              ##
################################################################################
d <- list(OBS = d.data_frame$OBS, EXP = d.data_frame$EXP,
          idx= d.data_frame$idx, idy= d.data_frame$idy)
for(i in 1:k){
  eval(parse(text = paste0("d$idxy.",i,"<-d.data_frame$idxy.",i) ))
  eval(parse(text = paste0("d$I",i,"<-d.data_frame$I",i) ))
}

################################################################################
## Adjacency and precision matrices                                           ##
################################################################################
############
## W: adjacency matrix (spatial)
############
adj<-spdep::poly2nb(carto)
W <- as(nb2mat(adj, style = "B"), "Matrix")

############
## Precision matrices for spatial random effects
############
## Spatial neighborhood matrix (Q_{xi})
spdep::nb2INLA("carto_nb.graph", spdep::poly2nb(carto))
g <- INLA::inla.read.graph("carto_nb.graph")
Qs <- matrix(0, g$n, g$n)
for (i in 1:g$n){
  Qs[i,i]=g$nnbs[[i]]
  Qs[i,g$nbs[[i]]]=-1
}

############
## W: adjacency matrix (temporal)
############
## W.t: adjacency matrix (temporal)
adj.t<-list(); adj.t[[1]]<-c(2); adj.t[[t]]<-c(t-1)
for(l in 2:(t-1)){adj.t[[l]]<-c(l-1,l+1)}
for(x in 1:length(adj.t)){adj.t[[x]]<-as.integer(adj.t[[x]])}
names(adj.t)<- 1:length(adj.t)
b=list(class="nb", call="nuestro_nb", cell=TRUE, rook=TRUE, sym=TRUE)
attributes(adj.t)<-b
W.t <- as(nb2mat(adj.t, style = "B"), "Matrix")


############
## Precision matrices for temporal random effects (RW1)
############
Qt<- crossprod(diff(diag(t), differences=1))

############
## Precision matrices for spatio-temporal interaction (generic0) = prec todos los crimes
############
R.1<- diag(n*t)
r.def.1<- 0
R.2<- kronecker(Qt,diag(n))
r.def.2<- n
R.3<- kronecker(diag(t),Qs)
r.def.3<- t
R.4<- kronecker(Qt, Qs)
r.def.4<- n+t-1

eval(parse(text = paste0("R.st<-R.",type)))
eval(parse(text = paste0("r.def.st<-r.def.",type) ))

################################################################################
## Constraints                                                                ##
################################################################################
## spatial
A.constr.s<- kronecker(diag(k), matrix(1,1,n))
e.s<- rep(0,dim(A.constr.s)[1])

## temporal
A.constr.t<- kronecker(diag(k), matrix(1,1,t)) 
e.t<- rep(0,dim(A.constr.t)[1])

## spatio-temporal
A.constr.st.1<- matrix(1, nrow=1, ncol=n*t)
e.st.1= rep(0,dim(A.constr.st.1)[1])

A.constr.st.2<- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
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
## Define appropriate hyperprior distributions (spatio-temporal)
sdunif="expression:
logdens=-log_precision/2;
return(logdens)"

## spalio-temporal M-models
source(paste0("functions/Mmodel_icar_",effects,".R"))
eval(parse(text = paste0("funct.model<- inla.rgeneric.Mmodel.icar.", effects) ))
model.s <- inla.rgeneric.define(funct.model, debug = TRUE, k = k, W = W)
model.t <- inla.rgeneric.define(funct.model, debug = TRUE, k = k, W = W.t)

## formula
formula.inla<- OBS ~ -1 + I1 + I2 + I3 + I4 +
  f(idx, model = model.s, constr=FALSE, extraconstr=list(A=A.constr.s, e=e.s)) +
  f(idy, model = model.t, constr=FALSE, extraconstr=list(A=A.constr.t, e=e.t) ) +
  f(idxy.1, model="generic0", Cmatrix=R.st, rankdef=r.def.st, hyper=list(prec=list(prior=sdunif)),
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st ))+
  f(idxy.2, model="generic0", Cmatrix=R.st, rankdef=r.def.st, hyper=list(prec=list(prior=sdunif)),
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st ))+
  f(idxy.3, model="generic0", Cmatrix=R.st, rankdef=r.def.st, hyper=list(prec=list(prior=sdunif)),
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st ))+
  f(idxy.4, model="generic0", Cmatrix=R.st, rankdef=r.def.st, hyper=list(prec=list(prior=sdunif)),
    constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st ))

## model
model<- inla(formula.inla,
             data = d, E = EXP,
             family = "poisson", control.predictor = list(compute = TRUE, cdf=c(log(1))),
             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
             control.fixed(list(mean=list(I1=0, I2=0, I3=0, I4=0), prec=list(I1=0.01, I2=0.01, I3=0.01, I4=0.01))),
             control.inla = list(strategy=strategy),
             verbose = FALSE)

## save
save(model, file=paste0("resul/resul_Mmodel_",effects,"_type",type,"_",gsub("\\.", "_", strategy),".RData"))

################################################################################
################################################################################