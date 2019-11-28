
#'MCMC for Sero
#'
#'@param pars_ini vector of initial parameters for seroprevalence
#'@param seroout processed serological data
#'@param foi_const_surv vector of additional foi for each survey (R0 model)
#'@param vc_agg3d aggregated vaccination array
#'@param pop_agg3d aggregated population array
#'@param pop_moments_agg population moments aggregated
#'@param dim_year years of interest
#'@param dim_age ages of interest
#'@param p_at_survey population proportion at time of serological survey
#'@param P_tot_survey total population at time of serological survey
#'@param inc_v3d_agg the incidence of vaccination (aggregated)
#'@param model_type either "R0" or "Foi"
#'@param parameter_type indexing for the different types of parameters
#'@param prob_Foi log of the model prior for Foi model
#'@param posterior_distributions distribution parameters for each model pseudo priors
#'@param ign sero surveys to ignore (good for bug hunting)
#'@param name_dir where to save the output
#'@param Niter number of mcmc iterations
#'@param run_id for parallel running- differentiates the different chains. Default = 1
#'
#'@export
#'@return saves the output in name_dir
Sero_Gibbs_MCMC = function(pars_ini,
                           seroout,
                           foi_const_surv,
                           vc_agg3d,
                           pop_agg3d,
                           pop_moments_agg,
                           dim_year,
                           dim_age,
                           p_at_survey,
                           P_tot_survey,
                           inc_v3d_agg,
                           model_type,
                           parameter_type,
                           prob_Foi,
                           posterior_distributions,
                           ign,
                           name_dir,
                           Niter,
                           run_id = 1) {

  burnin = min(Niter, 2*length(pars_ini))


  ### find posterior probability at start ###
  out = Sero_Gibbs_MCMC_step(pars_ini,
                             seroout,
                             foi_const_surv,
                             vc_agg3d,
                             pop_agg3d,
                             pop_moments_agg,
                             dim_year,
                             dim_age,
                             p_at_survey,
                             P_tot_survey,
                             inc_v3d_agg,
                             chain_cov = 1,
                             adapt = 0,
                             model_type,
                             parameter_type,
                             accCurrent=-Inf,
                             prob_Foi,
                             posterior_distributions,
                             ign)
  fileIndex = 0
  iter = 1
  chain = model_chain = posteriorProb = acceptRate = R0modelPrior = NULL

  for (iter in iter:Niter){
    #save current step
    param = out$param
    names(param) = names(pars_ini)
    model_type = out$model_type
    accCurrent = out$accCurrent
    accept = out$accept

    chain = rbind(chain, param)
    model_chain = rbind(model_chain, model_type == "R0")
    posteriorProb = rbind(posteriorProb, accCurrent)
    acceptRate = rbind(acceptRate, accept)


    if(iter>10 & run_id<100) plot(model_chain , type= "l") #plot(na.omit(posteriorProb), type="l") #

    if (iter %% 100 == 0){

      colnames(acceptRate) = "acceptRate"
      colnames(posteriorProb) = "posteriorProb"
      colnames(model_chain) = "model_chain"
      if (iter %% 10000 == 0){
        fileIndex  = iter/10000
      }

      out_state = cbind(chain, model_chain, posteriorProb, acceptRate, prob_Foi)[min((fileIndex * 10000+1),iter):iter,]

      write.csv(out_state,
                paste0(name_dir,"/","sero_chain_",fileIndex,"_output",run_id,".csv") )
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
    out = Sero_Gibbs_MCMC_step(param = param,
                               seroout,
                               foi_const_surv = foi_const_surv,
                               vc_agg3d = vc_agg3d,
                               pop_agg3d = pop_agg3d,
                               pop_moments_agg = pop_moments_agg,
                               dim_year = dim_year,
                               dim_age = dim_age,
                               p_at_survey,
                               P_tot_survey,
                               inc_v3d_agg,
                               chain_cov,
                               adapt,
                               model_type,
                               parameter_type,
                               accCurrent,
                               prob_Foi,
                               posterior_distributions,
                               ign)


  } #for loop end

} #function



#'MCMC step for Sero
#'
#'@param param vector of parameters for seroprevalence
#'@param seroout processed serology data
#'@param foi_const_surv vector of additional foi for each survey (R0 model)
#'@param vc_agg3d aggregated vaccination array
#'@param pop_agg3d aggregated population array
#'@param pop_moments_agg population moments aggregated
#'@param dim_year years of interest
#'@param dim_age ages of interest
#'@param p_at_survey population proportion at time of serological survey
#'@param P_tot_survey total population at time of serological survey
#'@param inc_v3d_agg the incidence of vaccination (aggregated)
#'@param chain_cov covariance of chain up to now
#'@param adapt whether to use the chain covariance to adapt the proposal distribution is 0 or 1
#'@param model_type either "R0" or "Foi"
#'@param parameter_type indexing for the different types of parameters
#'@param accCurrent current posterior probability
#'@param prob_Foi log of the model prior for Foi model
#'@param posterior_distributions distribution parameters for each model pseudo priors
#'@param ign sero surveys to ignore (good for bug hunting)
#'
#'@export
#'@return next step in Metropolis-within-gibbs (modout)
Sero_Gibbs_MCMC_step = function(param,
                                seroout,
                                foi_const_surv,
                                vc_agg3d,
                                pop_agg3d,
                                pop_moments_agg,
                                dim_year,
                                dim_age,
                                p_at_survey,
                                P_tot_survey,
                                inc_v3d_agg,
                                chain_cov,
                                adapt,
                                model_type,
                                parameter_type,
                                accCurrent,
                                prob_Foi,
                                posterior_distributions,
                                ign) {

  ############################# PARAMETER STEP ############################################################
  parout = Sero_Gibbs_MCMC_parameter_step (param = param,
                                           seroout = seroout,
                                           foi_const_surv = foi_const_surv,
                                           vc_agg3d = vc_agg3d,
                                           pop_agg3d = pop_agg3d,
                                           pop_moments_agg = pop_moments_agg,
                                           dim_year = dim_year,
                                           dim_age = dim_age,
                                           p_at_survey = p_at_survey,
                                           P_tot_survey = P_tot_survey,
                                           inc_v3d_agg = inc_v3d_agg,
                                           chain_cov = chain_cov,
                                           adapt = adapt,
                                           model_type = model_type,
                                           parameter_type = parameter_type,
                                           accCurrent = accCurrent,
                                           prob_Foi = prob_Foi,
                                           posterior_distributions = posterior_distributions,
                                           ign = ign)
  param = parout$param
  accept = parout$accept
  like = parout$like
  prior = parout$prior
  accCurrent = parout$accCurrent

  ############################# MODEL STEP ########################################################

  modout = Sero_Gibbs_MCMC_model_step (param = param,
                                       seroout = seroout,
                                       foi_const_surv = foi_const_surv,
                                       vc_agg3d = vc_agg3d,
                                       pop_agg3d = pop_agg3d,
                                       pop_moments_agg = pop_moments_agg,
                                       dim_year = dim_year,
                                       dim_age = dim_age,
                                       p_at_survey = p_at_survey,
                                       P_tot_survey = P_tot_survey,
                                       inc_v3d_agg = inc_v3d_agg,
                                       chain_cov = chain_cov,
                                       adapt = adapt,
                                       model_type = model_type,
                                       parameter_type = parameter_type,
                                       accCurrent = accCurrent,
                                       prob_Foi = prob_Foi,
                                       posterior_distributions = posterior_distributions,
                                       ign = ign,
                                       like = like,
                                       prior = prior)


  return(modout)
}


#'MCMC parameter step for Sero
#'
#'@param param vector of parameters for seroprevalence
#'@param seroout processed serology data
#'@param foi_const_surv vector of additional foi for each survey (R0 model)
#'@param vc_agg3d aggregated vaccination array
#'@param pop_agg3d aggregated population array
#'@param pop_moments_agg population moments aggregated
#'@param dim_year years of interest
#'@param dim_age ages of interest
#'@param p_at_survey population proportion at time of serological survey
#'@param P_tot_survey total population at time of serological survey
#'@param inc_v3d_agg the incidence of vaccination (aggregated)
#'@param chain_cov covariance of chain up to now
#'@param adapt whether to use the chain covariance to adapt the proposal distribution is 0 or 1
#'@param model_type either "R0" or "Foi"
#'@param parameter_type indexing for the different types of parameters
#'@param accCurrent current posterior probability
#'@param prob_Foi log of the model prior for Foi model
#'@param posterior_distributions distribution parameters for each model pseudo priors
#'@param ign sero surveys to ignore (good for bug hunting)
#'
#'@export
#'@return next step in the metropolis-within-gibbs for the parameters (param, like, prior, accept)
Sero_Gibbs_MCMC_parameter_step = function(param,
                                          seroout,
                                          foi_const_surv,
                                          vc_agg3d,
                                          pop_agg3d,
                                          pop_moments_agg,
                                          dim_year,
                                          dim_age,
                                          p_at_survey,
                                          P_tot_survey,
                                          inc_v3d_agg,
                                          chain_cov,
                                          adapt,
                                          model_type,
                                          parameter_type,
                                          accCurrent,
                                          prob_Foi,
                                          posterior_distributions,
                                          ign) {

  model_type_converse = switch(model_type,
                               "Foi" = "R0",
                               "R0" = "Foi")

  ### sample `off' parameters only if prob_Foi is not 0 or 1
  if (is.finite(prob_Foi) & prob_Foi != 0 ){

    pos_dis = posterior_distributions %>% dplyr::filter(model_type == model_type_converse)

    param[grep(model_type_converse,names(param))] = switch(model_type_converse,
                                                           "R0" = log( rtrunc(1,
                                                                               spec = "norm",
                                                                               a=1, b=Inf,
                                                                               mean = pos_dis[,1],
                                                                               sd = pos_dis[,2])),
                                                           "Foi" = log( rtrunc(1,
                                                                              spec = "norm",
                                                                              a=0, b=Inf,
                                                                              mean =pos_dis[,1],
                                                                              sd=pos_dis[,2] ) ))

  } #end if

  ### Accept/ reject `on' parameters ###
  param_prop = param
  out = SEROproposal(param, chain_cov, adapt, model_type)

  param_prop[grep(model_type_converse,names(param), invert = TRUE)] = out$param_prop[grep(model_type_converse,
                                                                                          names(param),
                                                                                          invert = TRUE)]

  ### priors ###
  prior_prop = sum( SEROprior(param_prop,  parameter_type, model_type) )
  prior = sum( SEROprior(param,  parameter_type, model_type) )


  ### if prior finite, evaluate likelihood ###
  if (is.finite(prior_prop)) {
    like_prop = sum(SEROlike( param_prop,
                              seroout,
                              foi_const_surv,
                              vc_agg3d,
                              pop_agg3d,
                              pop_moments_agg,
                              dim_year,
                              dim_age,
                              p_at_survey,
                              P_tot_survey,
                              inc_v3d_agg,
                              model_type,
                              ign), na.rm =TRUE)

    like = sum( SEROlike( param,
                          seroout,
                          foi_const_surv,
                          vc_agg3d,
                          pop_agg3d,
                          pop_moments_agg,
                          dim_year,
                          dim_age,
                          p_at_survey,
                          P_tot_survey,
                          inc_v3d_agg,
                          model_type,
                          ign), na.rm =TRUE)

    ### accept/ reject ###
    accProp = like_prop + prior_prop

    if(is.finite(like) & is.finite(prior)){
      accCurrent = like + prior
    } else {accCurrent = -Inf}

    p_accept= accProp - accCurrent
    # rm(accCurrent)
  } else {
    p_accept = log(0)
    like = -Inf
  }

  ## accept/reject step:
  accept = 0
  tmp = log(runif(1))
  if( tmp< min( (p_accept), log(1) ) ) { # accept:
    param = param_prop
    accept = 1
  }

  return(list(param = param, accept = accept, like = like, prior = prior, accCurrent = accCurrent))
}



#'MCMC model step for Sero
#'
#'@param param vector of parameters for seroprevalence
#'@param seroout processed serology data
#'@param foi_const_surv vector of additional foi for each survey (R0 model)
#'@param vc_agg3d aggregated vaccination array
#'@param pop_agg3d aggregated population array
#'@param pop_moments_agg population moments aggregated
#'@param dim_year years of interest
#'@param dim_age ages of interest
#'@param p_at_survey population proportion at time of serological survey
#'@param P_tot_survey total population at time of serological survey
#'@param inc_v3d_agg the incidence of vaccination (aggregated)
#'@param chain_cov covariance of chain up to now
#'@param adapt whether to use the chain covariance to adapt the proposal distribution is 0 or 1
#'@param model_type either "R0" or "Foi"
#'@param parameter_type indexing for the different types of parameters
#'@param accCurrent current posterior probability
#'@param prob_Foi log of the model prior for Foi model
#'@param posterior_distributions distribution parameters for each model pseudo priors
#'@param ign sero surveys to ignore (good for bug hunting)
#'@param like value of the likelihood from the Sero_Gibbs_MCMC_parameter_step
#'@param prior value of the prior from the Sero_Gibbs_MCMC_parameter_step
#'
#'@export
#'@return next step in the metropolis-within-gibbs for the model type (param, model_type, accCurrent, accept)
Sero_Gibbs_MCMC_model_step = function(param,
                                      seroout,
                                      foi_const_surv,
                                      vc_agg3d,
                                      pop_agg3d,
                                      pop_moments_agg,
                                      dim_year,
                                      dim_age,
                                      p_at_survey,
                                      P_tot_survey,
                                      inc_v3d_agg,
                                      chain_cov,
                                      adapt,
                                      model_type,
                                      parameter_type,
                                      accCurrent,
                                      prob_Foi,
                                      posterior_distributions,
                                      ign,
                                      like,
                                      prior) {

  if (is.finite(prob_Foi) & prob_Foi != 0 ){

    ## Propose new model ##
    out = SEROproposal(param, chain_cov, adapt, model_type)
    model_type_prop = out$model_type_prop


    ### priors ###
    prior_prop = sum( SEROprior(param, parameter_type, model_type_prop) )

    ### pseudopriors ###
    pseudoout = SEROpseudoprior(model_type,
                                param,
                                posterior_distributions,
                                prob_Foi)

    pseudoprior = pseudoout$pseudoprior
    modelprior = pseudoout$modelprior
    rm(pseudoout)

    pseudoout_prop = SEROpseudoprior(model_type_prop,
                                     param,
                                     posterior_distributions,
                                     prob_Foi)

    pseudoprior_prop = pseudoout_prop$pseudoprior
    modelprior_prop = pseudoout_prop$modelprior
    rm(pseudoout_prop)

    ### if prior finite, evaluate likelihood ###
    if (is.finite(prior_prop) & is.finite(modelprior_prop)) {

      like_prop = sum( SEROlike( param,
                                 seroout,
                                 foi_const_surv,
                                 vc_agg3d,
                                 pop_agg3d,
                                 pop_moments_agg,
                                 dim_year,
                                 dim_age,
                                 p_at_survey,
                                 P_tot_survey,
                                 inc_v3d_agg,
                                 model_type_prop,
                                 ign), na.rm =TRUE)

      ### accept/ reject ###
      accProp = like_prop + prior_prop + modelprior_prop + pseudoprior_prop


      if( is.finite(like) & is.finite(prior) & is.finite(modelprior)){
        accCurrent = like + prior + modelprior + pseudoprior
      } else {
        accCurrent = -Inf
      }

      p_accept= accProp - accCurrent
    } else {
      p_accept = log(0)

    }

    if(is.na(p_accept)){p_accept = log(0)} #if both infinite do not move

    ## accept/reject step:
    tmp = (runif(1))
    if( tmp< min( exp(p_accept), 1) ) { # accept:
      model_type = model_type_prop
      accCurrent = accProp
      accept = 1
    } else { # reject:

      accept = 0
    }
  } else if(is.infinite(prob_Foi)){
    model_type = "R0"
    accept = 1
  } else if(prob_Foi==0){
    model_type = "Foi"
    accept = 1
  }

  return(list(param = param, model_type = model_type, accCurrent = accCurrent, accept = accept))
}








#'Pseudo prior probability
#'
#'@param param vector of parameters for seroprevalence
#'@param model_type either "R0" or "Foi"
#'@param prob_Foi log of the model prior for Foi model
#'@param posterior_distributions distribution parameters for each model pseudo priors
#'
#'@export
#'@return pseudo prior probabilities (pseudoprior) and model prior probability (modelprior)
SEROpseudoprior = function(model_type,
                           param,
                           posterior_distributions,
                           prob_Foi){

  model_type_converse = switch(model_type,
                               "Foi" = "R0",
                               "R0" = "Foi")

  pos_dis = posterior_distributions %>% filter(model_type == model_type_converse)

  pseudoprior = switch(model_type,
                       "Foi" = sum(dtrunc(exp( param[grep(model_type_converse, names(param))] ),
                                          spec = "norm",
                                          a=1, b=Inf,
                                          mean =pos_dis[,1],
                                          sd=pos_dis[,2] ,
                                          log=TRUE)),
                       "R0" = sum(dtrunc(exp( param[grep(model_type_converse, names(param))] ),
                                         spec = "norm",
                                         a=0, b=Inf,
                                         mean =pos_dis[,1],
                                         sd=pos_dis[,2],
                                         log=TRUE)) )

  modelprior = switch(model_type,
                      "Foi" = prob_Foi,
                      "R0" = log(1-exp(prob_Foi) ) )


  return(list(pseudoprior = pseudoprior, modelprior = modelprior))
}


#'Proposal distribution for Sero
#'
#'@param param vector of parameters for seroprevalence
#'@param chain_cov covariance of chain up to now
#'@param adapt whether to use the chain covariance to adapt the proposal distribution is 0 or 1
#'@param model_type either "R0" or "Foi"
#'
#'@export
#'@return proposed parameters (param_prop) and proposed model type (model_type_prop)
SEROproposal = function(param, chain_cov, adapt, model_type) {
  #this adapts the proposal distribution covariance when adapt = 1
  no_param = length(param)

  if (adapt) {
    param_prop = rmvnorm(n = 1,
                         mean = param,
                         sigma = (2.38 ^ 2) * chain_cov / no_param
    ) #'optimal' scaling of chain covariance

  } else {
    sigma = (1e-2) ^ 2 * diag(no_param) / no_param #this is an inital proposal covariance, see [Mckinley et al 2014]
    param_prop = rmvnorm(n = 1,
                         mean = param,
                         sigma = sigma)

  }

  names(param_prop) = names(param)

  if (runif(1)<0.5){
    model_type_prop = "Foi"
  } else {
    model_type_prop = "R0"
  }
  return(list(param_prop = param_prop, model_type_prop = model_type_prop) )
}





#'Prior distribution for Sero
#'
#'@param param vector of parameters for seroprevalence
#'@param parameter_type indexing for the different types of parameters
#'@param model_type either "R0" or "Foi"
#'
#'@export
#'@return vector of prior probabilities
SEROprior = function(param,  parameter_type, model_type ) {

  #ADJUST FOR LOG TRANSFORM
  correction =  sum(param)
  params = exp(param)

  #declare
  Prior = rep(0,4)

  #vac eff
  Prior[1]= log( dtrunc(params[ which(parameter_type == 1)],"norm",a=0, b=1, mean = 0.975, sd = 0.05) )  #KEVINS PAPER

  #transmission parameters
  Prior[2] = switch(model_type,
                    "Foi" = sum(dexp(params[grep("Foi", names(params))] , rate = 0.001, log =TRUE)),
                    "R0" = sum(dexp(params[grep("R0", names(params))] - 1, rate = 0.001, log =TRUE)))

  #vc.factor
  Prior[3] =  dunif(params[which(parameter_type == 4)],
                    min = 0,
                    max = 1,
                    log = TRUE) #  vc.factor

  Prior[4] = correction

  out = as.numeric( Prior )
  #print(out)
  return( out )
}





#'Likelihood distribution for Sero
#'
#'@param param vector of parameters for seroprevalence
#'@param seroout processed serology data
#'@param foi_const_surv vector of additional foi for each survey (R0 model)
#'@param vc_agg3d aggregated vaccination array
#'@param pop_agg3d aggregated population array
#'@param pop_moments_agg population moments aggregated
#'@param dim_year years of interest
#'@param dim_age ages of interest
#'@param p_at_survey population proportion at time of serological survey
#'@param P_tot_survey total population at time of serological survey
#'@param inc_v3d_agg the incidence of vaccination (aggregated)
#'@param model_type either "R0" or "Foi"
#'@param ign indices of surveys to ignore in the likelihood- useful for bug hunting
#'
#'@export
#'@return log likelihood
SEROlike = function(param,
                    seroout,
                    foi_const_surv,
                    vc_agg3d,
                    pop_agg3d,
                    pop_moments_agg,
                    dim_year,
                    dim_age,
                    p_at_survey,
                    P_tot_survey,
                    inc_v3d_agg,
                    model_type,
                    ign) {

  no_sero_surveys = length(seroout$survey_dat)

  #sort out parameters
  vac_eff = exp(param[grep("vac_eff", names(param))]) # ADJUST FOR LOG TRANSFORM
  vcfac_CMRs = exp(param[grep("vc", names(param))]) # ADJUST FOR LOG TRANSFORM

  vcfac = seroout$vc_factor
  vcfac[seroout$sero_studies == "CMRs"] = vcfac_CMRs

  ### declare lnL
  lnL = rep(NA, no_sero_surveys)
  for (surveyIndex in 1:no_sero_surveys) {

    vc_agg_ag = NULL
    ag = findInterval(0:100, dplyr::pull(seroout$survey_dat[[surveyIndex]], age_min) )
    vc_agg_ag = aggregate(c(as.matrix(vc_agg3d[surveyIndex, paste0(seroout$study_years[surveyIndex]),]
                                      * pop_agg3d[surveyIndex, paste0(seroout$study_years[surveyIndex]),])),
                          by = list(ag = ag), sum)$x /
      aggregate(c(as.matrix(pop_agg3d[surveyIndex, paste0(seroout$study_years[surveyIndex]),])),
                by = list(ag = ag), sum)$x


    lnL[surveyIndex] = SEROlike_onesurvey(model_type,
                                          transmission_intensity =
                                            switch(model_type,
                                                   "R0" = exp(param[1 + no_sero_surveys + surveyIndex]),
                                                   "Foi" = exp(param[1 + surveyIndex])),
                                          seroout,
                                          pop = pop_agg3d[surveyIndex, paste0(seroout$study_years[surveyIndex]),],
                                          vc = vcfac[surveyIndex] * vac_eff * vc_agg_ag,
                                          surveyIndex = surveyIndex,
                                          foi_const = foi_const_surv[surveyIndex],
                                          dim_year,
                                          dim_age,
                                          p_at_survey,
                                          P_tot_survey,
                                          inc_v3d_agg,
                                          pop_moments_agg,
                                          vcfac_arg = vcfac[surveyIndex],
                                          vac_eff_arg = vac_eff)


  }


  #### IGNORE  ####
  lnL[c(ign)] = NaN


  return(lnL)
}



#' likelihood for one serological survey in either model
#'
#' @param model_type whether in R0 or Foi model
#' @param transmission_intensity value of foi or R0 for survey
#' @param seroout processed serology data
#' @param pop population in year and place of survey
#' @param age_groups1 age groups in one survey
#' @param vc vaccine coverage, defaults to 0
#' @param index_survey index of serological survey
#' @param foi_const constant force of infection to add to R0 model
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_at_survey_3d population proportion at survey as an array
#' @param P_tot_survey_2d population total at survey
#' @param inc_v3d_agg incidence of vacciantion at survey
#' @param pop_moments_agg aggregated population moments at survey
#' @param vcfac_arg vaccine coverage factor
#' @param vac_eff_arg vaccine efficacy
#'
#' @return log likelihood for one serological survey from either model
#' @export
SEROlike_onesurvey = function(model_type,
                              transmission_intensity,
                              seroout,
                              pop,
                              vc = 0,
                              surveyIndex,
                              foi_const = 0,
                              dim_year,
                              dim_age,
                              p_at_survey,
                              P_tot_survey,
                              inc_v3d_agg,
                              pop_moments_agg,
                              vcfac_arg = 1,
                              vac_eff_arg) {
  #sort seroout
  samples1 = dplyr::pull(seroout$survey_dat[[surveyIndex]], samples)
  positives1 = dplyr::pull(seroout$survey_dat[[surveyIndex]], positives)
  age_groups1 =  dplyr::pull(seroout$survey_dat[[surveyIndex]], age_min)
  t0_vac = seroout$t0_vac
  study_years = seroout$study_years


  if(length(positives1) !=length(samples1) |
     length(positives1) != length(age_groups1)){
    stop("SEROlike_onesurvey: incompatible dimensions_\n")
  }

  if(model_type == "Foi"){
    serop = foi_calculate_seroprevalence(foi = transmission_intensity,
                                         age_groups1,
                                         pop,
                                         vc)
  } else if(model_type == "R0"){
    serop = R0_calculate_seroprevalence(index_survey=surveyIndex,
                                        R0=transmission_intensity,
                                        foi_const,
                                        t0_vac,
                                        dim_year,
                                        dim_age,
                                        study_years,
                                        p_at_survey,
                                        P_tot_survey,
                                        inc_v3d_agg,
                                        age_groups=age_groups1,
                                        vcfac_arg=vcfac_arg,
                                        pop_moments_agg = pop_moments_agg,
                                        vac_eff_arg=vac_eff_arg)
  }
  #avoiding infinite Loglike
  log_serop = ifelse(serop==0, -50, log(serop))
  log_1_minus_serop = ifelse(serop==1, -50, log(1-serop))

  #rearranging
  Loglike = sum(lgamma(samples1+1)-lgamma(positives1+1)-lgamma(samples1-positives1+1) +
                  positives1*log_serop + (samples1-positives1)*log_1_minus_serop)


  return(Loglike)
}
