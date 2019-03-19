#'Proposal distribution for GLM
#'
#'@param param vector of parameters for GLM
#'@param chain_cov covariance of chain up to now
#'@param adapt whether to use the chain covariance to adapt the proposal distribution is 0 or 1
#'
#'@export
#'@return param_prop. The proposed parameters.
GLMproposal = function(param, chain_cov, adapt) {
  
  no_param = length(param)
  
  if (adapt) {
    sigma = (2.38 ^ 2) * chain_cov / no_param #'optimal' scaling of chain covariance
    param_prop = rmvnorm(n = 1, mean = param, sigma = sigma) 
  } else {
    sigma = (1e-2) ^ 2 * diag(no_param) / no_param #this is an inital proposal covariance, see [Mckinley et al 2014]
    param_prop = rmvnorm(n = 1, mean = param, sigma = sigma)
  }
  
  return(param_prop[1,])
}



#'Prior distributions for the GLM parameters
#'
#'@param param vector of parameters for GLM
#'
#'@export
#'@return vector of prior probabilities
GLMprior = function(param) {
  
  Prior = rep(0,2)
  
  #GLM
  jj = grep("^log.adm05", names(param)) # select country parameters (parameter_type=2) the sd.prior=2 is from Kevin's original code create status
  sd.prior = 2
  
  Prior[1] =  - 0.5 * sum((param[jj] / sd.prior) ^ 2) # adjustment for reduced variation between countries?
  
  Prior[2] =  sum(dnorm(param[grepl("^log.adm05", names(param)) == FALSE],
                        mean = 0,
                        sd = 30,
                        log = TRUE
  ))
  # this term is for normally distributed non-country parameters : normal distrib with high sd (flat distrib)
  out = as.numeric( Prior )
  return( out )
}






#'likelihood for the GLM parameters
#'
#'@param beta coefficients 
#'@param x covariates
#'@param y explanatory variable
#'
#'@export
#'@return log likelihood
GLMlike = function(beta, x, y) {
  # beta are the model coefficients,
  # x the independent covariates (needs one column for the Intercept),
  # and y the binary outcome_
  # model predictions pi = p, 1-pi = q
  eta = as.numeric(x %*% beta)
  logq = -exp(eta) # -exp(X*beta) (= log(1-q) )
  logp = log(1-exp(logq)) # 
  #q = exp(logq)
  #p = 1-q
  logl = sum(logp[y==1]) + sum(logq[y==0])
  return (logl)
}





#'mcmc step for the GLM parameters
#'
#'@param param vector of parameters for GLM
#'@param x covariates
#'@param y explanatory variable
#'@param chain_cov covariance of chain up to now
#'@param adapt whether to use the chain covariance to adapt the proposal distribution is 0 or 1
#'@param accCurrent current posterior probability
#'
#'@export
#'@return new param, new posterior probability, whether step was accepted or not
GLM_MCMC_step = function(param, x, y, chain_cov, adapt,  accCurrent) {
  
  ### propose new param ###
  param_prop = GLMproposal(param, chain_cov, adapt)
 
  ### priors ###
  prior_prop = sum( GLMprior(param_prop) )
  
  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    
    like_prop = GLMlike(param_prop, x, y)
    
    ### accept prob ###
    accProp = like_prop + prior_prop

    p_accept= accProp - accCurrent
    if(is.na(p_accept) ){ 
      p_accept = -Inf
      }
  } else {
    p_accept = log(0)
  }
  
  ## accept/reject step:
  tmp = (runif(1))
  if( tmp< min( exp(p_accept), 1)  ) { # accept:
    param = param_prop
    accCurrent = accProp
    accept = 1
  } else { # reject:
    accept = 0
  }
  
  return(list(param = param,  accCurrent = accCurrent, accept = accept)  )
}

#'mcmc for the GLM parameters
#'
#'@param pars_ini vector of inital parameters for GLM
#'@param x covariates
#'@param y explanatory variable
#'@param Niter number of interations
#'@param plot_chain whether to plot the first parameter's chain (TRUE or FALSE)
#'@param name_dir where to save output
#'
#'@export
#'@return saves chain in a folder of name paste0(name_dir,"/","GLM_chain",fileIndex,"_output.csv")
GLM_MCMC = function(Niter, name_dir, pars_ini, x, y, plot_chain){
  
  burnin = min(2*length(pars_ini), Niter)
  
  ### find posterior probability at start ###
  out = GLM_MCMC_step(pars_ini, x, y, chain_cov=1, adapt=0, accCurrent=-Inf)
  
  chain = posteriorProb = acceptRate = NULL
  
  fileIndex = 0
  iter = 1
  
  for (iter in iter:Niter){
    #save current step
    param = out$param;  names(param) = names(pars_ini);  accCurrent = out$accCurrent;  accept = out$accept
    
    chain = rbind(chain, param);  
    posteriorProb = rbind(posteriorProb, accCurrent)
    acceptRate = rbind(acceptRate, accept)
    
    if(iter>100 & plot_chain) plot(chain[,1] , type= "l") 
    
    if (iter %% 1000 == 0){
      
      colnames(acceptRate) = "acceptRate"
      colnames(posteriorProb) = "posteriorProb"
      if (iter %% 10000 == 0){
        fileIndex  = iter/10000
      }
      
      out_state = cbind(chain,  posteriorProb, acceptRate)[min((fileIndex * 10000+1),iter):iter,]
      
      write.csv(out_state, paste0(name_dir,"/","GLM_chain",fileIndex,"_output.csv") ) 
    }
    
    #adapt?
    if (iter>burnin & runif(1)<0.9){ #adapt
      adapt = 1
      chain_cov  = cov(chain[max(nrow(chain)-10000, 1):nrow(chain),])
    } else {
      adapt = 0
      chain_cov = 1
    }
    
    #new step
    out = GLM_MCMC_step(param, x, y, chain_cov, adapt, accCurrent)
    
  }
  
}
