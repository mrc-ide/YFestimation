#' function estimate a value for the model prior in order to have even mixing.
#' Note, this is based on initial parameter values so may need adjusting.
#' 
#' @param pars_ini vector of initial parameters for seroprevalence 
#' @param seroout processed serology data
#' @param foi_const_surv vector of additional foi for each survey (R0 model)
#' @param vc_agg3d aggregated vaccination array
#' @param pop_agg3d aggregated population array
#' @param pop_moments_agg population moments aggregated
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_at_survey population proportion at time of serological survey
#' @param P_tot_survey total population at time of serological survey
#' @param inc_v3d_agg the incidence of vaccination (aggregated)
#' @param parameter_type indexing for the different types of parameters
#' @param ign indices of surveys to ignore in the likelihood- useful for bug hunting
#' 
#' @export
#' @return Estimate of log model prior for Foi model

setup_modelprior = function(pars_ini,
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
                            parameter_type,
                            ign){
  
  
  # calculate likelihoods at start
  likeR0 =  sum(SEROlike(pars_ini,
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
                         model_type = "R0",
                         ign), na.rm = TRUE)
  
  likeFoi =  sum(SEROlike(pars_ini,
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
                          model_type = "Foi",
                          ign), na.rm = TRUE)
  
  #calculate priors at start
  priorR0 = sum( SEROprior(pars_ini,  parameter_type, "R0") )
  
  priorfoi = sum( SEROprior(pars_ini,  parameter_type, "Foi") )
  
  #calculate posterior for each
  postR0 = likeR0 + priorR0
  
  postFoi = likeFoi + priorfoi
  
  
  #solve for the model prior for foi model (log of)
  #this is assuming that we want the posteriors to be 
  #approximately equal. The initial parameters need to 
  #be fairly good for this to be useful ie. starting
  #in an area of high probability for both models
  prob_Foi = log( 1/(exp(postFoi - postR0)-1))
  
  
  return(prob_Foi)
}


#' function estimate a value for the pseudo prior in order to have even mixing.
#' Note, this is based on initial parameter values so may need adjusting.
#' 
#' @param pars_ini vector of initial parameters for seroprevalence 
#' @param posterior_distributions distribution parameters for each model pseudo priors
#' @param scale_model_type the model whose posterior distribution sd is to be scaled
#' @param tolerance how close the difference between the pseudo prior probabilities 
#'                  shoud be at the initial parameters. Defaults to 0.01
#' 
#' 
#' @export
#' @return Estimate of scaling factor of scale_model_type posterior distributions sd- 
#'         assuming they are normally distributed

setup_pseudoprior = function(pars_ini,
                             posterior_distributions,
                             scale_model_type,
                             tolerance = 0.01){
  
  prob_Foi = log(0.5) #does not matter what this is
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
  #first iteration- two ends of the interval
  scaling_factor_long = seq(from = 1e-10, to = 5, length.out = 100)
  difference = Inf
  iter = 0
  
  #get starting interval so as to not miss anything
  difference_vec = NULL
  for(sf in scaling_factor_long){
    #set new scaling factor
    pos_dis = posterior_distributions
    pos_dis$sd[pos_dis$model_type == scale_model_type] = 
      sf*pos_dis$sd[pos_dis$model_type == scale_model_type]
    
    #calculate new pseudo prior
    pseudoR0 = sum(SEROpseudoprior("R0", pars_ini, pos_dis, prob_Foi)$pseudoprior)
    pseudoFoi = sum(SEROpseudoprior("Foi", pars_ini, pos_dis, prob_Foi)$pseudoprior)
    
    difference_vec = rbind(difference_vec, pseudoFoi - pseudoR0)
  }
  
  #sign change 
  sign_change = which(diff(sign(difference_vec)) != 0)
  
  #set up interval
  scaling_factor = c(scaling_factor_long[sign_change],
                     (scaling_factor_long[sign_change]+
                        scaling_factor_long[sign_change+1])/2,
                     scaling_factor_long[sign_change+1])
  ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
  
  while(iter<100 & abs(difference)>tolerance){
    
    difference_vec = NULL
    for(sf in scaling_factor){
      #set new scaling factor
      pos_dis = posterior_distributions
      pos_dis$sd[pos_dis$model_type == scale_model_type] = 
        sf*pos_dis$sd[pos_dis$model_type == scale_model_type]
      
      #calculate new pseudo prior
      pseudoR0 = sum(SEROpseudoprior("R0", pars_ini, pos_dis, prob_Foi)$pseudoprior)
      pseudoFoi = sum(SEROpseudoprior("Foi", pars_ini, pos_dis, prob_Foi)$pseudoprior)
      
      difference_vec = rbind(difference_vec, pseudoFoi - pseudoR0)
    }
    
    #assign new boundaries of the interval, first check boundaries are good
    if(all(sign(difference_vec)>0) | all(sign(difference_vec)<0)){
      stop("interval not wide enough")
    }
    
    difference = difference_vec[2]
    out_scaling = scaling_factor[2]
    iter = iter+1
    
    if(sign(difference_vec[3]) != sign(difference_vec[2])){
      scaling_factor[1] = scaling_factor[2] 
      scaling_factor[2] = (scaling_factor[1] + scaling_factor[3])/2
    } else {
      scaling_factor[3] = scaling_factor[2]
      scaling_factor[2] = (scaling_factor[1] + scaling_factor[3])/2
    }
    
    if(iter == 100){stop("maximum number of iterations reached")}
  }
  return(out_scaling)
}