################################################################################
######          Univariate spatio-temporal models (iCAR - RW1)            ######
################################################################################
rm(list=ls())

## libraries
library(INLA); library(spdep); library(spData)

### save results
if(!file.exists("resul")) {dir.create("resul")}

################################################################################
## Selection strategy
################################################################################
############
## Selected prior distributions for the coefficients of the temporal random effects (RW1 or RW2)
############
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
## Precision matrices for spatial, temporal and spatio-temporal interaction   ##
## random effects                                                             ##
################################################################################
## Spatial
spdep::nb2INLA("carto_nb.graph", spdep::poly2nb(carto))
g<- INLA::inla.read.graph("carto_nb.graph")
Qs<- matrix(0, g$n, g$n)
for(i in 1:g$n){
  Qs[i,i]=g$nnbs[[i]]
  Qs[i,g$nbs[[i]]]=-1
}

## Temporal
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

A.constr.st.2<- kronecker(matrix(1,1,t),diag(n))
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
for(mm in 1:k){
  ######################################
  ## Data by crime
  ######################################
  d.data_frame<- data[,c(crimes[mm],e_crimes[mm],"ID_area","ID_year")]
  d.data_frame <- d.data_frame[order(d.data_frame$ID_year, d.data_frame$ID_area),]
  rownames(d.data_frame)<-NULL
  d.data_frame$ID_area_year<- seq(1,length(unique(d.data_frame$ID_year))*length(unique(d.data_frame$ID_area)))
  colnames(d.data_frame)<- c("OBS", "EXP", "ID_area", "ID_year", "ID_area_year")
  
  ## data.inla
  Data<- data.frame(O=d.data_frame$OBS, E=d.data_frame$EXP,
                    ID.area= d.data_frame$ID_area,
                    ID.year=d.data_frame$ID_year,
                    ID.area.year=d.data_frame$ID_area_year)
  
  ######################################
  ## Fit models                       ##
  ######################################
  # type<- 2
  resulta<- list()
  for(type in 1:4){
    eval(parse(text = paste0("R.st<-R.",type)))
    eval(parse(text = paste0("r.def.st<-r.def.",type) ))
    eval(parse(text = paste0("A.constr.st<- A.constr.st.",type) ))
    eval(parse(text = paste0("e.st<- e.st.",type) ))
    
    ## formula
    formula.inla<- O ~  f(ID.area, model="besag", graph=g, constr=TRUE, hyper=list(prec=list(prior=sdunif)) ) +
                        f(ID.year, model=paste0("rw",ord.temp), constr=TRUE, hyper=list(prec=list(prior=sdunif))) +
                        f(ID.area.year, model="generic0", Cmatrix=R.st, rankdef=r.def.st, hyper=list(prec=list(prior=sdunif)),
                          constr=FALSE, extraconstr=list(A=A.constr.st, e=e.st))
    
    ## model
    resulta[[type]]<- INLA::inla(formula.inla, family="poisson", data=Data, E=E,
                                 control.predictor=list(compute=TRUE, cdf=c(log(1))),
                                 control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                                 control.inla=list(strategy=strategy, npoints=21),
                                 verbose = FALSE)
    
    ## rm
    rm(list = c("R.st", "r.def.st","A.constr.st","e.st","formula.inla"))
  }
  names(resulta)<- c("univ.t1", "univ.t2", "univ.t3", "univ.t4")
  
  
  ## save
  save(resulta, file = paste0("resul/resul_univ_iCAR_rw",ord.temp,"_",crimes[[mm]],"_",gsub("\\.", "_", strategy),".RData"))
  
  ## rm
  rm(list = c("d.data_frame","Data"))
}
################################################################################
################################################################################