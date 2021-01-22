########################################################################
## Fits several Bayesian spatio-temporal P-spline models              ##
########################################################################
ps_st <- function(carto=NULL,
                   data=NULL,
                   ID.area=NULL,
                   ID.year=NULL,
                   crimes=NULL,
                   Expcrimes=NULL,
                   prior.spatial=1,      ## Selected prior distributions for the coefficients of the  spatial P-splines (1: RW1 ; 2: RW2)
                   prior.temporal=1,     ## Selected prior distributions for the coefficients of the  temporal P-splines (1: RW1 ; 2: RW2)
                   prior.interaction=2,  ## Selected prior distributions for the espatio-temporal random effects (1: Type I; 2: Type II; 3: Type III; 4: Type IV)
                   order.B=3,            ## order B-splines
                   k.long=NULL,          ## Number of internal intervals (longitude)
                   k.lati=NULL,           ## Number of internal intervals (latitude)
                   k.time=NULL,          ## Number of internal intervals (time)
                   centered=TRUE,        ## TRUE: to center the smooth functions; FALSE: to center the coefficients
                   strategy="simplified.laplace",
                   rerun=FALSE
                   ){
  
  ## Check for errors
  ## ----------------------------------------
  if(is.null(carto))
    stop("the carto argument is missing")
  if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
    stop("the carto argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
  if(is.null(data))
    stop("the data argument is missing")
  if(is.null(ID.area))
    stop("the ID.area argument is missing")
  if(is.null(crimes))
    stop("the crimes argument is missing")
  if(is.null(Expcrimes))
    stop("the Expcrimes argument is missing")
  
  if(!(prior.spatial %in% c(1:2)))
    stop("invalid prior.spatial argument")
  if(!(prior.temporal %in% c(1:2)))
    stop("invalid prior.temporal argument")
  if(!(prior.interaction %in% c(1:4)))
    stop("invalid prior.interaction argument")
  
  if(!(strategy %in% c("gaussian","simplified.laplace","laplace","adaptative")))
    stop("invalid strategy argument")

  
  ## ----------------------------------------
  ## Pre-processing data
  ## ----------------------------------------
  cat("STEP 1: Pre-processing data\n")
  
  ## Number of areas and number of time periods
  n <- length(unique(data[,ID.area]))
  t <- length(unique(data[,ID.year]))
  
  ## Parameter of the model (number of crimes)
  k <- length(crimes)
  
  ## Row-wisw Kronecker product
  Rten <- function(X1,X2){
    one1 <- matrix(1,1,ncol(X1))
    one2 <- matrix(1,1,ncol(X2))
    kronecker(X1,one2)*kronecker(one1,X2)
  }
  
  
  ## Construct the spatial B-spline basis (Bs) 
  ## ----------------------------------------
  
  ## Longitudes scale into the [0,1] interval
  x1 <- sp::coordinates(carto)[,1]
  x1 <- (x1-min(x1))/(max(x1)-min(x1))
  if(is.null(k.long)){ndx1 <- ceiling(min(length(unique(x1))/4, 40))} else {ndx1 <- k.long}
  dis1 <- (max(x1)-min(x1))/ndx1
  x1l <- min(x1)-dis1*0.05
  x1r <- max(x1)+dis1*0.05
  dx1 <- (x1r-x1l)/ndx1
  knots1 <- seq(x1l-order.B*dx1, x1r+order.B*dx1, by=dx1)
  
  ## The horizontal B-spline basis
  B1 <- splines::spline.des(knots1,x1,order.B+1, 0*x1)$design
  k1 <- dim(B1)[2]
  
  ## Latitudes scale into the [0,1] interval
  x2 <- sp::coordinates(carto)[,2]
  x2 <- (x2-min(x2))/(max(x2)-min(x2))
  if(is.null(k.lati)){ndx2 <- ceiling(min(length(unique(x2))/4, 40))} else {ndx2 <- k.lati}
  ndx2 <- ceiling(min(length(unique(x2))/4, 40))
  dis2 <- (max(x2)-min(x2))/ndx2
  x2l <- min(x2)-dis2*0.05
  x2r <- max(x2)+dis2*0.05
  dx2 <- (x2r-x2l)/ndx2
  knots2 <- seq(x2l-order.B*dx2, x2r+order.B*dx2, by=dx2)
  
  ## The vertical B-spline basis
  B2 <- splines::spline.des(knots2,x2,order.B+1, 0*x2)$design
  k2 <- dim(B2)[2]
  
  ## The spatial B-spline basis
  Bs <- Rten(B2,B1)
  ks<- dim(Bs)[2]
  
  
  ## Construct the temporal B-spline basis (Bt) 
  ## ----------------------------------------
  ## Times scale into the [0,1] interval
  xt <- unique(data[,ID.year])
  xt <- (xt-min(xt))/(max(xt)-min(xt))
  if(is.null(k.time)){ndxt <- ceiling(min(length(unique(xt))/4, 40))} else {ndxt <- k.time}
  dis.t <- (max(xt)-min(xt))/ndxt
  xtl <- min(xt)- dis.t*0.05
  xtr <- max(xt)+ dis.t*0.05
  dxt <- (xtr-xtl)/ndxt
  knots.t <- seq(xtl- order.B*dxt, xtr+order.B*dxt, by=dxt)
  
  ## The temporal B-spline basis
  Bt <- splines::spline.des(knots.t,xt,order.B+1, 0*xt)$design
  kt <- dim(Bt)[2]
  
  
  
  ## Adjacency and precision matrices
  ## ----------------------------------------

  ## Spatial structure matrices
  P1 <- crossprod(diff(diag(k1), differences=prior.spatial))
  P2 <- crossprod(diff(diag(k2), differences=prior.spatial))

  I.P1 <- kronecker(Matrix::Diagonal(nrow(P2),1), P1)
  P2.I <- kronecker(P2, Matrix::Diagonal(nrow(P1),1))
  
  Cmat.s <- list(inla.as.sparse(I.P1),inla.as.sparse(P2.I))
  
  
  ## Temporal structure matrix
  Pt <- crossprod(diff(diag(kt),differences=prior.temporal))
  Cmat.t <- list(inla.as.sparse(Pt))
  
  
  ## R.st: Precision matrices for spatio-temporal interaction
  spdep::nb2INLA("carto_nb.graph", spdep::poly2nb(carto))
  g <- INLA::inla.read.graph("carto_nb.graph")
  Qs <- matrix(0, g$n, g$n)
  for(i in 1:g$n){
    Qs[i,i]=g$nnbs[[i]]
    Qs[i,g$nbs[[i]]]=-1
  }
  Qt <- crossprod(diff(diag(t), differences=prior.temporal))


  if(prior.interaction %in% c(1)){
    R.st <- diag(n*t)
    r.def.st <- 0
  }
  if(prior.interaction %in% c(2)){
    R.st <- kronecker(Qt,diag(n))
    r.def.st <- prior.temporal*n
  }
  if(prior.interaction %in% c(3)){
    R.st <- kronecker(diag(t),Qs)
    r.def.st <- t
  }
  if(prior.interaction %in% c(4)){
    R.st <- kronecker(Qt, Qs)
    r.def.st <- prior.temporal*n+t-1
  }
  
  
  ## Define appropriate constraints matrices (A.constr)
  ## ----------------------------------------
  if(centered){ 
    A.constr.s<- matrix(apply(Bs,2,sum), nrow=1, ncol=ks)
    A.constr.t<- matrix(apply(Bt,2,sum), nrow=1, ncol=kt)
  } else {
    A.constr.s <- matrix(1, nrow=1, ncol=ks)
    A.constr.t <- matrix(1, nrow=1, ncol=kt)
  }
  
  
  if(prior.interaction %in% c(1)){
    A.constr.st <- matrix(1, nrow=1, ncol=n*t)
  }
  if(prior.interaction %in% c(2)){
    A.constr.st <- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
    if(prior.spatial==2 | prior.temporal==2){
      A.constr.st <- rbind(A.constr.st, kronecker(matrix(1:t,nrow=1, ncol=t), matrix(c(1, rep(0,n-1)),nrow=1, ncol=n)))
    }
  }
  if(prior.interaction %in% c(3)){
    A.constr.st <- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
  }
  if(prior.interaction %in% c(4)){
    A1 <- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
    if(prior.spatial==2 | prior.temporal==2){
      A1 <- rbind(A1, kronecker(matrix(1:t,nrow=1, ncol=t), matrix(c(1, rep(0,n-1)),nrow=1, ncol=n)))
    }
    A2 <- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
    A.constr.st <- rbind(A1[-1,], A2)
  }
  
  ## Design matrices for random effects (D.matrix)
  ## ----------------------------------------
  
  ## D.inter
  D.inter <- matrix(1, nrow = n*t, ncol=1)
  
  ## D.spat = I_{J} (x) 1_{t} (x) Bs
  D.spat <- kronecker(matrix(1, nrow = t, ncol=1), Bs)
  
  ## D.temp = I_{J} (x) Bt (x) 1_{n}
  D.temp <- kronecker(Bt, matrix(1, nrow=n, ncol=1))
  
  ## D.spat.temp
  D.spat.temp <- diag(1,t*n)
  
  ## Design matrix
  D.matrix <- cbind(D.inter, D.spat, D.temp,D.spat.temp)
  
  
  ## Define appropriate hyperprior distributions (spatio-temporal)
  ## ----------------------------------------
  sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"
  
  
  ## ----------------------------------------
  ## Fitting model
  ## ----------------------------------------
  
  cat("STEP 2: Fitting model with INLA (this may take a while...)\n")
  
  
  ## To compute patterns (as linear combinations of the log-risks)
  ## ----------------------------------------
  lc.spat.univ <- INLA::inla.make.lincombs(idx = Bs)
  names(lc.spat.univ) <- paste("spatial.",as.character(seq(1,n,1)),sep="")
  
  lc.temp.univ <- INLA::inla.make.lincombs(idy = Bt)
  names(lc.temp.univ) <- paste("temporal.",as.character(seq(1,t,1)),sep="")
  
  all.lc <- c(lc.spat.univ,lc.temp.univ)
  
  
  
  ## Implementation of the spatio-temporal P-spline models
  ## ----------------------------------------
  ## Define the 'formula'
  form <- "O ~ -1 +"
  form <- paste(form, paste("I", collapse="+"), sep=" ")
  form <- paste(form, "+ f(idx, model = 'generic3', Cmatrix=Cmat.s, rankdef=prior.spatial, constr=FALSE, extraconstr=list(A=A.constr.s, e=rep(0,dim(A.constr.s)[1])),
    hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)) )", sep=" ")
  form <- paste(form, "+ f(idy, model = 'generic3', Cmatrix=Cmat.t, rankdef=prior.temporal, constr=FALSE, extraconstr=list(A=A.constr.t, e=rep(0,dim(A.constr.t)[1])),
    hyper=list(prec1=list(prior=sdunif) ))", sep=" ")
  form <- paste(form, "+ f(idxy, model = 'generic0', Cmatrix=R.st, rankdef=r.def.st, constr=FALSE, extraconstr=list(A=A.constr.st, e=rep(0,dim(A.constr.st)[1])),
    hyper=list(prec=list(prior=sdunif) ))", sep=" ")
  formula.inla <- stats::as.formula(form)
  
  
  
  ## Fit the INLA model
  ## ----------------------------------------
  Models <- vector("list",k)
  names(Models) <- crimes
  # l<-4
  for(l in 1:k){
    ## Data for INLA
    df.inla <- data[,c(ID.area,ID.year,crimes[l],Expcrimes[l])]
    names(df.inla) <- c("ID.s","ID.t","O","E")
    df.inla <- df.inla[order(df.inla$ID.t, df.inla$ID.s),]
    Data.INLA <- list(O = df.inla$O, E = df.inla$E,
                      I= c(1, rep(NA, ks + kt + n*t)),
                      idx = c(rep(NA,1), 1:ks, rep(NA, kt + n*t)),
                      idy = c(rep(NA, 1 + ks), 1:kt, rep(NA, n*t)),
                      idxy = c(rep(NA, 1 + ks + kt), 1:(n*t))
                      )
    
    ## Fit the INLA model
    Models[[l]] <- INLA::inla(formula.inla, data = Data.INLA, E = E, family = "poisson",
                             control.predictor= list(compute = TRUE, A= D.matrix, link=1, cdf=c(log(1))),
                             control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                             lincomb=all.lc,
                             control.inla = list(strategy=strategy),
                             debug = F,
                             verbose = F)
    if(rerun){ Models[[l]]<- INLA::inla.rerun(Models[[l]]) }
    
    Models[[l]]$model <- list(prior.spatial=paste0("RW",prior.spatial),
                             prior.temporal=paste0("RW",prior.temporal),
                             prior.interaction=paste0("Type ",prior.interaction),
                             centered=centered,
                             rerun=rerun,
                             crime=crimes[[l]]
                             )
    rm(list = c("df.inla","Data.INLA"))
  }
  return(Models)
}
########################################################################
########################################################################