########################################################################
## M-models (INLA)
########################################################################
'mult.temp.Psplines' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL){
    ########################################################################
    ## theta
    ########################################################################
    interpret.theta = function(){
      sigma.j<- sapply(theta[as.integer(1:k)], function(x) { exp(-0.5*x) })
      corre.j<- sapply(theta[as.integer((k+1):(k*(k+1)/2))], function(x){
        (2* exp(x))/(1+exp(x)) -1})
      
      Rho<- diag(1,k) 
      Rho[lower.tri(Rho, diag = FALSE)] <- corre.j
      Rho[upper.tri(Rho, diag = FALSE)] <- t(Rho)[upper.tri(Rho)]
      
      sigma.mat<- matrix(sigma.j, ncol = 1) %*% matrix(sigma.j, nrow = 1)
      
      Covar<- Rho * sigma.mat
      Prec<- solve(Covar)
      
      sigma.lambda<- exp(-(1/2)*theta[as.integer(k*(k+1)/2+1)])
      
      return (list(sigma.j = sigma.j, corre.j=corre.j, Covar=Covar, Prec=Prec, sigma.lambda=sigma.lambda))
    }
    ########################################################################
    ## Graph of precision function; i.e., a 0/1 representation of precision matrix
    ########################################################################
    graph = function(){
      prec<- matrix(1, nrow = k, ncol=k)
      P<- Matrix::Diagonal(nrow(W), apply(W, 1, sum))+W
      G <- kronecker(prec, P)
      return (G)
    }
    ########################################################################
    ## Precision matrix
    ########################################################################
    Q = function(){
      param <- interpret.theta()
      P<- (param$sigma.lambda^(-2) ) * (Matrix::Diagonal(nrow(W), apply(W, 1, sum)) - W)
      Q<- kronecker(param$Prec, P)
      return (Q)
    }
    ########################################################################
    ## Mean of model
    ########################################################################
    mu = function(){
      return(numeric(0))
    }
    ########################################################################
    ## log.norm.const
    ########################################################################
    log.norm.const = function() {
      val <- numeric(0)
      return (val)
    }
    ########################################################################
    ## log.prior
    ########################################################################
    log.prior = function() {
      ## return the log-prior for the hyperparameters.
      param = interpret.theta()
      ##############################
      ## Wishart prior for Covar
      ##############################
      sigma2 <- 1 # Wishart parameter
      val = log(MCMCpack::dwish(W = param$Covar, v = 2*k+1, S = diag(rep(sigma2, k))))
      ##############################
      ## sigma (1:J)
      ##############################
      val = val + sum(-(5/2)* theta[as.integer(1:k)])
      ##############################
      ## rho
      ##############################
      val = val + sum( log(2) + theta[as.integer( (k+1):(k*(k+1)/2) )] - 2*log(1 + exp(theta[as.integer( (k+1):(k*(k+1)/2) )])) )
      ##############################
      ## sigma lambda ~ Unif(0,100) - smoothing parameter
      ##############################
      val = val + log(0.005)- (0.5)* theta[as.integer(k*(k+1)/2+1)]
      return (val)
    }
    ########################################################################
    ## initial
    ########################################################################
    initial = function() {
      ## return initial values
      return ( c(rep(1,k), rep(0, k*(k-1)/2), 0.5) )
    }
    ########################################################################
    ########################################################################
    quit = function() {
      return (invisible())
    }
    val = do.call(match.arg(cmd), args = list())
    return (val)
  }
########################################################################
########################################################################