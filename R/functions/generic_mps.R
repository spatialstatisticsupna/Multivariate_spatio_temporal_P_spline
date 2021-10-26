########################################################################
## Multivariate spatial Psplines (Whishart prior)
########################################################################
'generic_mps_spat' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL){
    ###############################
    ## cache
    ###############################
    envir = parent.env(environment())
    if (!exists("cache.done", envir = envir)) {
      M.W2 <- ( Matrix::Diagonal(nrow(W2), apply(W2, 1, sum)) - W2)
      M.W1 <- ( Matrix::Diagonal(nrow(W1), apply(W1, 1, sum)) - W1)
      assign("M.W2", M.W2, envir = envir)
      assign("M.W1", M.W1, envir = envir)
      assign("cache.done", TRUE, envir = envir)
    }
    ###############################
    ## theta
    ###############################
    interpret.theta = function(){
      sigma.j <- sapply(theta[as.integer(1:k)], function(x) { exp(-0.5*x) })
      corre.j <- sapply(theta[as.integer(k+1:(k*(k-1)/2))], function(x){ (2* exp(x))/(1+exp(x)) -1})
      
      Rho <- diag(1,k) 
      Rho[lower.tri(Rho, diag = FALSE)] <- corre.j
      Rho[upper.tri(Rho, diag = FALSE)] <- t(Rho)[upper.tri(Rho)]
      
      sigma.mat <- matrix(sigma.j, ncol = 1) %*% matrix(sigma.j, nrow = 1)
      Covar <- Rho * sigma.mat
      Prec <- solve(Covar)
      
      sigma.lambda <- sapply(theta[as.integer(k*(k+1)/2+1)], function(x) { exp(-0.5*x) })
      
      return (list(sigma.j = sigma.j, corre.j=corre.j, Covar=Covar, Prec=Prec, sigma.lambda=sigma.lambda))
    }
    ###############################
    ## Graph of precision function; i.e., a 0/1 representation of precision matrix
    ###############################
    graph = function(){
      return (Q())
    }
    ###############################
    ## Precision matrix
    ###############################
    Q = function(){
      param <- interpret.theta()
      P <- (param$sigma.lambda^(-2))* M.W2 + M.W1
      Q <- kronecker(param$Prec, P)
      Q <- inla.as.sparse(Q)
      return (Q)
    }
    ###############################
    ## Mean of model
    ###############################
    mu = function(){
      return(numeric(0))
    }
    ###############################
    ## log.norm.const
    ###############################
    log.norm.const = function() {
      val <- numeric(0)
      return (val)
    }
    ###############################
    ## log.prior
    ###############################
    log.prior = function() {
      ## return the log-prior for the hyperparameters
      param = interpret.theta()
      ##############
      ## Wishart prior for Covar
      ##############
      sigma2 <- 1 # Wishart parameter
      val = log(MCMCpack::dwish(W = param$Covar, v = 2*k+1, S = diag(rep(sigma2, k))))
      ##############
      ## sigma (1:J)
      ##############
      val = val + sum(-(5/2)* theta[as.integer(1:k)])
      ##############
      ## rho
      ##############
      val = val + sum( log(2) + theta[as.integer( k+1:(k*(k-1)/2) )] - 
                         2*log(1 + exp(theta[as.integer( k+1:(k*(k-1)/2) )])) )
      ##############
      ## sigma lambda ~ Unif(0,1000) - smoothing parameter
      ##############
      val = val + sum(log(0.0005)- (0.5)* theta[as.integer(k*(k+1)/2+1)])
      return (val)
    }
    ###############################
    ## initial
    ###############################
    initial = function() {
      ## return initial values
      return( as.vector(initial.values) )
    }
    ###############################
    ###############################
    quit = function() {
      return (invisible())
    }
    if (!length(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    
    return (val)
  }
########################################################################
########################################################################

########################################################################
## Multivariate temporal Psplines (Whishart prior)
########################################################################
'generic_mps_temp' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL){
    ###############################
    ## cache
    ###############################
    envir = parent.env(environment())
    if (!exists("cache.done", envir = envir)) {
      M.W <- Matrix::Diagonal(nrow(W), apply(W, 1, sum)) - W
      assign("M.W", M.W, envir = envir)
      assign("cache.done", TRUE, envir = envir)
    }
    ###############################
    ## theta
    ###############################
    interpret.theta = function(){
      sigma.j <- sapply(theta[as.integer(1:k)], function(x) { exp(-0.5*x) })
      corre.j <- sapply(theta[as.integer(k+1:(k*(k-1)/2))], function(x){ (2* exp(x))/(1+exp(x)) -1})
      
      Rho <- diag(1,k) 
      Rho[lower.tri(Rho, diag = FALSE)] <- corre.j
      Rho[upper.tri(Rho, diag = FALSE)] <- t(Rho)[upper.tri(Rho)]
      
      sigma.mat <- matrix(sigma.j, ncol = 1) %*% matrix(sigma.j, nrow = 1)
      Covar <- Rho * sigma.mat
      Prec <- solve(Covar)
      
      return (list(sigma.j = sigma.j, corre.j=corre.j, Covar=Covar, Prec=Prec))
    }
    ###############################
    ## Graph of precision function; i.e., a 0/1 representation of precision matrix
    ###############################
    graph = function(){
      return (Q())
    }
    ###############################
    ## Precision matrix
    ###############################
    Q = function(){
      param <- interpret.theta()
      P <- M.W
      Q <- kronecker(param$Prec, P)
      Q <- inla.as.sparse(Q)
      return (Q)
    }
    ###############################
    ## Mean of model
    ###############################
    mu = function(){
      return(numeric(0))
    }
    ###############################
    ## log.norm.const
    ###############################
    log.norm.const = function() {
      val <- numeric(0)
      return (val)
    }
    ###############################
    ## log.prior
    ###############################
    log.prior = function() {
      ## return the log-prior for the hyperparameters
      param = interpret.theta()
      ##############
      ## Wishart prior for Covar
      ##############
      sigma2 <- 0.1 # Wishart parameter
      val = log(MCMCpack::dwish(W = param$Covar, v = 2*k+1, S = diag(rep(sigma2, k))))
      ##############
      ## sigma (1:J)
      ##############
      val = val + sum(-(5/2)* theta[as.integer(1:k)])
      ##############
      ## rho
      ##############
      val = val + sum( log(2) + theta[as.integer( k+1:(k*(k-1)/2) )] - 
                         2*log(1 + exp(theta[as.integer( k+1:(k*(k-1)/2) )])) )
      return (val)
    }
    ###############################
    ## initial
    ###############################
    initial = function() {
      ## return initial values
      return( as.vector(initial.values) )
    }
    ###############################
    ###############################
    quit = function() {
      return (invisible())
    }
    if (!length(theta)) theta = initial()
    val = do.call(match.arg(cmd), args = list())
    return (val)
  }
########################################################################
########################################################################