

#'calculate the probability of detection for both models
#'
#' @param x glm covariates
#' @inheritParams ii
#' @param seroout processed serology data
#' @param params estimated parameters for both serology and GLM from either R0 or Foi model
#' @param dat environmental and occurrence data
#' @param t0_vac_africa time of first vaccination for admins in Africa endemic zone
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_prop_3d prpoprtion of population in array form
#' @param P_tot_2d total population in 2d form
#' @param inc_vc3d incidence of vaccination in array form
#' @param pop_moments_whole population moments for all admins
#' @param varsin_nc covariates in use
#' @param vc30_agg vaccine coverage aggregated over time and space for observation period
#' @param pop30_agg population aggregated over time and space for observation period
#' @param model_type either "R0" or "Foi"
#'
#' @return the probability of detection
#' @export
fun_calc_pdetect_multi_both = function(x,
                                       ii,
                                       seroout,
                                       params,
                                       dat,
                                       t0_vac_africa,
                                       dim_year,
                                       dim_age,
                                       p_prop_3d,
                                       P_tot_2d,
                                       inc_v3d,
                                       pop_moments_whole,
                                       varsin_nc,
                                       vc30_agg,
                                       pop30_agg,
                                       model_type){

  adm1s = seroout$adm1s
  sero_studies = seroout$sero_studies
  no_sero_surveys = seroout$no_sero_surveys

  # calculate the b term from the left part of equation 6(plos Med paper)
  adm = match(unlist(adm1s), vc30_agg$adm0_adm1)
  t0_vac_adm = t0_vac_africa[adm]

  # now I have to calculate the Number of infections over observation period among each PROVINCE covered by a survey
  vac_eff = params[1]

  if(model_type == "R0"){

    R0_vec = NULL
    for(k in 1:no_sero_surveys) { # repeat the numbers of time of each adm1s
      R0_vec = c(R0_vec,
                 rep(params[which(names(params)==paste0("R0_", sero_studies[k]))],
                     length(adm1s[[k]])) )
    }

    res = R0_recurrence_seroprev_whole(adm = adm,
                                       dat = dat,
                                       R0 = R0_vec,
                                       t0_vac_adm,
                                       dim_year,
                                       dim_age,
                                       p_prop_3d,
                                       P_tot_2d,
                                       inc_v3d,
                                       pop_moments_whole,
                                       vac_eff_arg = vac_eff)

    # res returns the total Nb of infection on the whole period, I only need 1984 to 2019
    Ninf = rowSums( res$Ninf_t_province[,which(dim_year==1984):which(dim_year==2019)] )
    Ninf = ifelse(Ninf<1, 1, Ninf)

  } else if(model_type == "Foi"){

    vec_foi = vec_vc_factor = NULL
    for(i in 1:no_sero_surveys) {
      vec_foi = c(vec_foi,
                  rep(params[which(names(params)==paste("Foi_",sero_studies[i],sep=""))],
                      length(adm1s[[i]])))
      if(!is.na(seroout$vc_factor[i])) {
        vec_vc_factor = c(vec_vc_factor, rep(seroout$vc_factor[i],length(adm1s[[i]])))
      } else {
        vec_vc_factor = c(vec_vc_factor,
                          rep(params[which(names(params)==paste("vc_factor_",sero_studies[i],sep=""))],
                              length(adm1s[[i]])))
      }
    }

    vec_foi = vec_foi[match(vc30_agg[adm,1],unlist(adm1s))]
    vec_vc_factor = vec_vc_factor[match(vc30_agg[adm,1], unlist(adm1s))]

    popvc = (1-vac_eff*vec_vc_factor*vc30_agg[adm,-(1)])*pop30_agg[adm,-(1)]
    expfoi = exp(outer(-vec_foi, 0:100))
    Ninf = vec_foi*rowSums(expfoi*popvc)
  }

  #calcualte glm predictions
  mypreds_nc = fun_calcPred(coefs = params[ii],
                            newdata=x,
                            type="link",
                            varsin=varsin_nc)

  #calculate b
  b = mypreds_nc[adm]-log(Ninf)
  names(b) = unlist(adm1s)

  return( b )
}





#'calculate the probability of detection for both models
#'
#'@param x glm covariates
#'@param ii indices of glm parameters
#'@param seroout processed serology data
#'@param params estimated parameters for both serology and GLM from either R0 or Foi model
#'@param dat environmental and occurrence data
#'@param t0_vac_africa time of first vaccination for admins in Africa endemic zone
#'@param dim_year years of interest
#'@param dim_age ages of interest
#'@param p_prop_3d prpoprtion of population in array form
#'@param P_tot_2d total population in 2d form
#'@param inc_vc3d incidence of vaccination in array form
#'@param pop1 population data
#'@param vc2d vaccination coverages in 2d form
#'@param varsin_nc covariates in use
#'@param polydeg degree of polynomial, defaults to 5
#'@param R0_lookup table of R0 and infection numbers for each admin
#'@param model_type either "R0" or "Foi"
#'
#'@return transmission_whole, transmission intensty across Africa endemic zone
#'@export
fun_calc_transmission_Africa = function(x,
                                        ii,
                                        seroout,
                                        params,
                                        dat,
                                        t0_vac_africa,
                                        dim_year,
                                        dim_age,
                                        p_prop_3d,
                                        P_tot_2d,
                                        inc_v3d,
                                        pop1,
                                        vc2d,
                                        varsin_nc,
                                        polydeg=5,
                                        R0_lookup,
                                        model_type) {

  #get aggregated vc and pop over observation period
  aggout=create_pop30_agg_vc30_agg(pop1, vc2d)

  #glm predictions
  mypreds_nc  = fun_calcPred(coefs = as.numeric(params)[ii],
                             newdata=x,
                             type="link",
                             varsin=varsin_nc)

  #probability of detection
  p_detect =  fun_calc_pdetect_multi_both(x,
                                          ii,
                                          seroout,
                                          params,
                                          dat,
                                          t0_vac_africa,
                                          dim_year,
                                          dim_age,
                                          p_prop_3d,
                                          P_tot_2d,
                                          inc_v3d,
                                          pop_moments_whole,
                                          varsin_nc,
                                          aggout$vc30_agg,
                                          aggout$pop30_agg,
                                          model_type)


  p_detect_link = mean(p_detect)

  #calculating number of infections over the observation period for the whole region
  Ninf_whole = exp( mypreds_nc - p_detect_link)

  #translate this number into a transmission intensity for each model type
  if(model_type == "Foi"){
    pop_vc_moments = aggout$pop_vc_moments

    if(polydeg>ncol(pop_vc_moments)) error("fun_calc_transmission_Africa: invalid value for polydeg.\n")

    z = -Ninf_whole

    if(polydeg>0) for(i in 1:polydeg) {
      z = cbind(z,(-1)^(i+1)*pop_vc_moments[,i+1]/factorial(i-1))
    }

    transmission_whole = sapply(1:nrow(x), function(i) polyroot(z[i,]))
    transmission_whole[abs(Arg(transmission_whole))<=1e-10] = Re(transmission_whole)[abs(Arg(transmission_whole))<=1e-10]
    transmission_whole[abs(Arg(transmission_whole))>1e-10] = NA

    dt = dim(transmission_whole)
    transmission_whole = as.numeric(transmission_whole)
    dim(transmission_whole) = dt
    transmission_whole = apply(transmission_whole,2,min,na.rm=T)

  } else if(model_type == "R0"){

    transmission_whole = rep(NA, length(Ninf_whole))

    for (i in 1: length(t0_vac_africa) ){
      inf_bound = max( findInterval(floor(Ninf_whole[i]), R0_lookup[i,], left.open = TRUE) ,
                       1)
      sup_bound = inf_bound + 1

      if( sup_bound <= ncol(R0_lookup)){
        x_R0_1 = R0_lookup[i,inf_bound]# find inf bound
        x_R0_2 = R0_lookup[i,sup_bound] # sup bound

        y_R0_1 = as.numeric(colnames(R0_lookup)[inf_bound]) # find inf bound
        y_R0_2 = as.numeric(colnames(R0_lookup)[sup_bound]) # sup bound

        ## manually calculate the linear interpolation
        a_lin_inter = (y_R0_2 - y_R0_1)/(x_R0_2-x_R0_1)
        b_lin_inter = y_R0_1 - a_lin_inter*x_R0_1

        transmission_whole[i] = a_lin_inter*Ninf_whole[i] + b_lin_inter

      } else if (sup_bound  > ncol(R0_lookup)){
        transmission_whole[i] = as.numeric(colnames(R0_lookup)[ncol(R0_lookup)])
      }
    }

  }

  return(transmission_whole)
}
