########################################################################
## M-models (INLA)
########################################################################
'inla.rgeneric.Mmodel.icar.fe' <-
  function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                   "log.prior", "quit"), theta = NULL){
    ########################################################################
    ## theta
    ########################################################################
    interpret.theta = function(){
      ## The k*k parameters are the entries in the M matrix, by cols.
      M <- matrix(theta, ncol = k)
      return (list(M = M))
    }
    ########################################################################
    ## Graph of precision function; i.e., a 0/1 representation of precision matrix
    ########################################################################
    graph = function(){
      MI <- kronecker(Matrix(1, ncol = k, nrow = k), Diagonal(nrow(W), 1))
      IW <- Diagonal(nrow(W), 1) + W
      BlockIW <- bdiag(replicate(k, IW, simplify = FALSE))
      G <- (MI %*% BlockIW) %*% MI
      return (G)
    }
    ########################################################################
    ## Precision matrix
    ########################################################################
    Q = function(){
      param <- interpret.theta()
      M.inv <- solve(param$M)
      MI <- kronecker(M.inv, Diagonal(nrow(W), 1))
      D <- as.vector(apply(W, 1, sum))
      BlockIW <- bdiag(lapply(1:k, function(i){
        Diagonal(x = D) -  W
      }))
      Q <- (MI %*% BlockIW) %*% kronecker(t(M.inv), Diagonal(nrow(W), 1))
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
      ## Wishart prior for t(M)*M
      ##############################
      sigma2 <- 1000 # Wishart parameter
      val = log(MCMCpack::dwish(W = crossprod(param$M), v = k, S = diag(rep(sigma2, k))))
      return (val)
    }
    ########################################################################
    ## initial
    ########################################################################
    initial = function() {
      ## return initial values
      return ( c(as.vector(diag(rep(1, k)))))
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