#' function to calculate seroprevalence by recurrence relation in R0 model,
#' forms the basis of R0_recurrence_seroprev_survey and R0_recurrence_seroprev_whole
#'
#'
#' @param R0 value of R0 for survey
#' @param foi foi before vaccination
#' @param foi_const constant force of infection to add to R0 model
#' @param time_limit time of vaccination
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_at population proportion
#' @param P_tot population total
#' @param inc_v incidence of vaccination
#' @param vac_eff_arg vaccine efficacy
#'
#' @return all statuses in R0 model
#' @export
R0_recurrence_seroprev = function(R0,
                                  foi,
                                  foi_const = 0,
                                  time_limit,
                                  dim_year,
                                  dim_age,
                                  p_at,
                                  P_tot,
                                  inc_v,
                                  vac_eff_arg) {

  # foi_const is a constant FOI for environmental transmissions
  ageVec=c(0:100)
  length_year=length(dim_year)
  length_age=length(dim_age)

  # first create the table to keep the results
  s_atau1=rep(NA, length_year*length_age)
  dim(s_atau1)= c(length_year, length_age)
  Stot = Ninf_t = Nenv_t = rep(NA, length_year)
  n_inf = s_atau3 = s_atau2=  s_atau1bis= n_env = s_atau1

  # create the 4 compartiments accounting for natural infection i (+/-) and vaccination v (+/-)
  i_neg_v_neg_tau1 = i_neg_v_neg_tau2 = i_neg_v_neg_tau3 =  # i- v-
    i_pos_v_neg_tau1 = i_pos_v_neg_tau2 = i_pos_v_neg_tau3 = # i+ v-
    i_neg_v_pos_tau1 = i_neg_v_pos_tau2 = i_neg_v_pos_tau3 = # i- v+
    i_pos_v_pos_tau1 = i_pos_v_pos_tau2 = i_pos_v_pos_tau3 = matrix( rep(NA, length_year*length_age), nrow=length_year) # i+ v+

  i_neg_v_neg_tau1_bis = i_neg_v_pos_tau1_bis = i_pos_v_neg_tau1_bis = i_pos_v_pos_tau1_bis =
    matrix( rep(NA, length(dim_year)*length_age), nrow=length_year) # this is for environmental infection

  # then calculate the rest of the prevalence

  for (y in  1:length(dim_year) ){

    if ( dim_year[y] < time_limit ){
      # equation (4) at t+tau1
      i_neg_v_neg_tau1[y,] = exp( - (foi+foi_const)*ageVec )
      i_pos_v_neg_tau1[y,] = 1 - exp( - (foi+foi_const)*ageVec )
      i_neg_v_pos_tau1[y,] = 0 # no vaccination
      i_pos_v_pos_tau1[y,] = 0 # no vaccination

      # aging
      i_neg_v_neg_tau2[y,] = i_neg_v_neg_tau3[y,] = exp( - (foi+foi_const)*(ageVec+1) )
      i_pos_v_neg_tau2[y,] = i_pos_v_neg_tau3[y,] = 1- exp( - (foi+foi_const)*(ageVec+1) )
      i_neg_v_pos_tau2[y,] = i_neg_v_pos_tau3[y,] = 0
      i_pos_v_pos_tau2[y,] = i_pos_v_pos_tau3[y,] = 0

      #record infections
      n_inf[y,]=  (i_pos_v_neg_tau2[y,] - i_pos_v_neg_tau1[y,]) * p_at[y,] * P_tot[y]
      Ninf_t[y] = sum(n_inf[y,], na.rm = TRUE)
      n_env[y,] = (1-exp(-foi_const)) *  i_neg_v_neg_tau1[y,] * p_at[y,] * P_tot[y]
      Nenv_t[y] = sum(n_env[y,], na.rm = TRUE)
      Stot[y] =   sum(i_neg_v_neg_tau1[y,] * p_at[y,], na.rm = TRUE)

    } else if (dim_year[y] == time_limit){

      i_neg_v_neg_tau1[y,] = exp( - (foi+foi_const)*ageVec )
      i_pos_v_neg_tau1[y,] = 1 - exp( - (foi+foi_const)*ageVec )
      i_neg_v_pos_tau1[y,] = 0 # no vaccination
      i_pos_v_pos_tau1[y,] = 0 # no vaccination

      # aging
      i_neg_v_neg_tau2[y,] = exp( - (foi+foi_const)*(ageVec+1) )
      i_pos_v_neg_tau2[y,] = 1 - exp( - (foi+foi_const)*(ageVec+1) )
      i_neg_v_pos_tau2[y,] = 0
      i_pos_v_pos_tau2[y,] = 0

      #record infections
      n_inf[y,]=  (i_pos_v_neg_tau2[y,] - i_pos_v_neg_tau1[y,]) * p_at[y,] * P_tot[y] #same as before
      Ninf_t[y] = sum(n_inf[y,], na.rm = TRUE)#same as before
      n_env[y,] = (1-exp(-foi_const)) *  i_neg_v_neg_tau1[y,] * p_at[y,] * P_tot[y] #same as before
      Nenv_t[y] = sum(n_env[y,], na.rm = TRUE) #same as before
      Stot[y] =   sum(i_neg_v_neg_tau1[y,] * p_at[y,], na.rm = TRUE)#same as before

      #tau3: vaccination arise
      i_neg_v_neg_tau3[y,] = i_neg_v_neg_tau2[y,] * (1 - vac_eff_arg * inc_v[y,]) # eq (8)
      i_pos_v_neg_tau3[y,] = i_pos_v_neg_tau2[y,] * (1 - vac_eff_arg * inc_v[y,]) # eq (8)
      i_neg_v_pos_tau3[y,] = i_neg_v_pos_tau2[y,] + i_neg_v_neg_tau2[y,] * (vac_eff_arg * inc_v[y,])
      i_pos_v_pos_tau3[y,] = i_pos_v_pos_tau2[y,] + i_pos_v_neg_tau2[y,] * (vac_eff_arg * inc_v[y,])

    } else if (dim_year[y] > time_limit){

      # children born susceptible and unvaccinated
      i_neg_v_neg_tau1[y,1] = 1
      i_pos_v_neg_tau1[y,1] = 0
      i_neg_v_pos_tau1[y,1] = 0
      i_pos_v_pos_tau1[y,1] = 0

      for (a in 2:length_age){
        # tau1: aging   eq(4)
        i_neg_v_neg_tau1[y,a] = i_neg_v_neg_tau3[y-1,a-1]
        i_pos_v_neg_tau1[y,a] = i_pos_v_neg_tau3[y-1,a-1]
        i_neg_v_pos_tau1[y,a] = i_neg_v_pos_tau3[y-1,a-1]
        i_pos_v_pos_tau1[y,a] = i_pos_v_pos_tau3[y-1,a-1]
      }
      # tau1_bis: environmental infection
      i_neg_v_neg_tau1_bis[y,] = exp( - (foi_const)*ageVec )*i_neg_v_neg_tau1[y,]
      i_pos_v_neg_tau1_bis[y,] = i_pos_v_neg_tau1[y,] + (i_neg_v_neg_tau1[y,] - i_neg_v_neg_tau1_bis[y,])
      i_neg_v_pos_tau1_bis[y,] = i_neg_v_pos_tau1[y,] # no vaccination
      i_pos_v_pos_tau1_bis[y,] = i_pos_v_pos_tau1[y,]  # no vaccination

      n_env[y,] = (1-exp(-foi_const))*i_neg_v_neg_tau1[y,]*p_at[y,]*P_tot[y]
      Nenv_t[y] = sum(n_env[y,], na.rm = TRUE)

      if(dim_year[y] < 1990) {
        Stot[y] = sum(i_neg_v_neg_tau1_bis[y,] * p_at[y,], na.rm=TRUE) # total susceptibles are only those uninfected and unvaccinated
      } else if (dim_year[y] >= 1990){ # from 1990, we have ages btw 77 and 100y, what we don't have before
        Stot[y] = sum(i_neg_v_neg_tau1_bis[y, !is.na(i_neg_v_neg_tau1_bis[y,])] * p_at[y, !is.na(i_neg_v_neg_tau1_bis[y,])])/
                                                  sum(p_at[y, !is.na(i_neg_v_neg_tau1_bis[y,])])
      }
      Ninf_t[y] = max (0, P_tot[y]*(Stot[y] - 1/R0) ) # eq (5)
      n_inf[y,] = Ninf_t[y] * i_neg_v_neg_tau1_bis[y,] * p_at[y,]/ Stot[y] # eq (6)

      #tau2: new infections arise
      i_neg_v_neg_tau2[y,] = i_neg_v_neg_tau1_bis[y,] - n_inf[y,]/ (p_at[y,] * P_tot[y] ) # eq (7)
      i_pos_v_neg_tau2[y,] = i_pos_v_neg_tau1_bis[y,] + n_inf[y,]/ (p_at[y,] * P_tot[y] )
      i_neg_v_pos_tau2[y,] = i_neg_v_pos_tau1_bis[y,]
      i_pos_v_pos_tau2[y,] = i_pos_v_pos_tau1_bis[y,]

      #tau3: vaccination arise
      i_neg_v_neg_tau3[y,] = i_neg_v_neg_tau2[y,] * (1-vac_eff_arg * inc_v[y,]) # eq (8)
      i_pos_v_neg_tau3[y,] = i_pos_v_neg_tau2[y,] * (1-vac_eff_arg * inc_v[y,]) # eq (8)
      i_neg_v_pos_tau3[y,] = i_neg_v_pos_tau2[y,] + i_neg_v_neg_tau2[y,] * (vac_eff_arg * inc_v[y,])
      i_pos_v_pos_tau3[y,] = i_pos_v_pos_tau2[y,] + i_pos_v_neg_tau2[y,] * (vac_eff_arg * inc_v[y,])
    }

  } # end for (y in  1:length_year )


  return(list(i_neg_v_neg_tau3 = i_neg_v_neg_tau3,
              i_pos_v_neg_tau3 = i_pos_v_neg_tau3,
              i_neg_v_pos_tau3 = i_neg_v_pos_tau3,
              i_pos_v_pos_tau3 = i_pos_v_pos_tau3,
              Ninf_t = Ninf_t,
              Stot = Stot,
              Nenv_t = Nenv_t,
              n_inf = n_inf))
}




#' function to calculate seroprevalence everywhere by recurrence relation in R0 model
#'
#' @param adm admin 1 unit of interest
#' @param dat environmental and occurrence data
#' @param R0 value of R0 for survey
#' @param t0_vac_adm time of vaccination for admin, scalar
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_prop_3d population proportion
#' @param P_tot_2d population total
#' @param inc_v3d incidence of vaccination
#' @param pop_moments_whole aggregated population moments
#' @param vac_eff_arg vaccine efficacy
#'
#' @return all statuses in R0 model for admin
#' @export
R0_recurrence_seroprev_whole = function(adm,
                                        dat,
                                        R0,
                                        t0_vac_adm,
                                        dim_year,
                                        dim_age,
                                        p_prop_3d,
                                        P_tot_2d,
                                        inc_v3d,
                                        pop_moments_whole,
                                        vac_eff_arg) {
  #
  #THIS NOW ACTS AS A WRAPPER FOR fun_recurrence_seroprev WHICH CONTAINS THE one COPY OF THE EQUATIONS
  #COMBINING THE FUNCTIONS IN THIS WAY SHOULD REDUCE INCONSISTENCIES BETWEEN DIFFERENT MODEL VERSIONs
  #SUCH AS: THE INCLUSION OF ENVIRONMENTAL TRANSMISSION

  # Declare outputs
  n_inf_at = s_atau3 = s_atau2= s_atau1= Stot_at = Ninf_t = NULL

  # create the 4 compartments accounting for natural infection i (+/-) and vaccination v (+/-)
  i_neg_v_neg_tau3 = i_pos_v_neg_tau3 = i_neg_v_pos_tau3 = i_pos_v_pos_tau3 = NULL

  # I firstly need to calculate for each adm the FOI at t0_vac
  pop_mom=NULL
  for(j in 1:length(adm)){
    pop_mom = rbind (pop_mom, pop_moments_whole[adm[j], ])
  }
  pop_mom = cbind(dat$adm0_adm1[adm], pop_mom)

  #get foi before vaccination
  foi = YFburden::foi_prevac(adm = adm, R0=R0, pop_moments = pop_mom, polydeg=6)

  for (i in 1:length(adm) ){

    time_limit = min(1981, t0_vac_adm[i])  #THIS IS CONSISTENT?
    p_at= p_prop_3d[adm[i], ,]
    P_tot=P_tot_2d[adm[i], ]
    inc_v=inc_v3d[adm[i], ,]

    seroprev_out = R0_recurrence_seroprev(R0 = R0[i],
                                          foi = foi[i],
                                          foi_const = 0,
                                          time_limit = time_limit,
                                          dim_year = dim_year,
                                          dim_age = dim_age,
                                          p_at = p_at,
                                          P_tot = P_tot,
                                          inc_v = inc_v,
                                          vac_eff_arg = vac_eff_arg
    )

    if(i==1){

      i_neg_v_neg_tau3 = seroprev_out$i_neg_v_neg_tau3 #FIRST JUST ASSIGN MATRIX
      i_pos_v_neg_tau3 = seroprev_out$i_pos_v_neg_tau3
      i_neg_v_pos_tau3 = seroprev_out$i_neg_v_pos_tau3
      i_pos_v_pos_tau3 = seroprev_out$i_pos_v_pos_tau3
      Ninf_t = seroprev_out$Ninf_t
      Stot_at = seroprev_out$Stot
      n_inf_at = seroprev_out$n_inf

    } else if(i==2){

      i_neg_v_neg_tau3 = abind::abind(i_neg_v_neg_tau3, seroprev_out$i_neg_v_neg_tau3, along=0) #THEN BIND ALONG NEW FIRST DIMENSION, RESULT IS AN ARRAY
      i_pos_v_neg_tau3 = abind::abind(i_pos_v_neg_tau3, seroprev_out$i_pos_v_neg_tau3, along=0)
      i_neg_v_pos_tau3 = abind::abind(i_neg_v_pos_tau3, seroprev_out$i_neg_v_pos_tau3, along=0)
      i_pos_v_pos_tau3 = abind::abind(i_pos_v_pos_tau3, seroprev_out$i_pos_v_pos_tau3, along=0)
      Ninf_t = abind::abind(Ninf_t, seroprev_out$Ninf_t,along=0)
      Stot_at = abind::abind(Stot_at, seroprev_out$Stot,along=0)
      n_inf_at = abind::abind(n_inf_at, seroprev_out$n_inf,along=0)

    } else {

      i_neg_v_neg_tau3 = abind::abind(i_neg_v_neg_tau3, seroprev_out$i_neg_v_neg_tau3, along=1) #NOW FIRST DIMENSION EXISTS, BIND ALONG THAT
      i_pos_v_neg_tau3 = abind::abind(i_pos_v_neg_tau3, seroprev_out$i_pos_v_neg_tau3, along=1)
      i_neg_v_pos_tau3 = abind::abind(i_neg_v_pos_tau3, seroprev_out$i_neg_v_pos_tau3, along=1)
      i_pos_v_pos_tau3 = abind::abind(i_pos_v_pos_tau3, seroprev_out$i_pos_v_pos_tau3, along=1)
      Ninf_t = abind::abind(Ninf_t, seroprev_out$Ninf_t,along=1)
      Stot_at = abind::abind(Stot_at, seroprev_out$Stot,along=1)
      n_inf_at = abind::abind(n_inf_at, seroprev_out$n_inf,along=1)

    }

  }

  return(list(i_neg_v_neg_tau3 = i_neg_v_neg_tau3,
              i_pos_v_neg_tau3 = i_pos_v_neg_tau3,
              i_neg_v_pos_tau3 = i_neg_v_pos_tau3,
              i_pos_v_pos_tau3 = i_pos_v_pos_tau3,
              Stot_at = Stot_at,
              Ninf_t_province = Ninf_t,
              n_inf_at_province = n_inf_at))
}




#' function to calculate seroprevalence at survey by recurrence relation in R0 model
#'
#' @param R0 value of R0 for survey
#' @param foi_const constant force of infection to add to R0 model
#' @param t0_vac time of vaccination, scalar
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_at_survey population proportion at one survey
#' @param P_tot_survey population total at one survey
#' @param inc_v3d incidence of vaccination at one survey (from inv_v3d_agg)
#' @param pop_moments aggregated population moments at one survey (from pop_moments_agg)
#' @param vac_eff_arg vaccine efficacy
#'
#' @return all statuses in R0 model for one survey
#' @export
R0_recurrence_seroprev_survey = function(R0,
                                         foi_const,
                                         t0_vac,
                                         dim_year,
                                         dim_age,
                                         p_at_survey,
                                         P_tot_survey,
                                         inc_v3d,
                                         pop_moments,
                                         vac_eff_arg) {


  #THIS NOW ACTS AS A WRAPPER FOR fun_recurrence_seroprev WHICH CONTAINS THE one COPY OF THE EQUATIONS
  #COMBINING THE FUNCTIONS IN THIS WAY SHOULD REDUCE INCONSISTENCIES BETWEEN DIFFERENT MODEL VERSIONs
  #SUCH AS: THE INCLUSION OF ENVIRONMENTAL TRANSMISSION

  #calculate foi
  foi = YFburden::foi_prevac(adm=NA, R0=R0, pop_moments=pop_moments, polydeg=6)

  #time_limit
  time_limit=t0_vac

  #run recurrence
  seroprev_out = R0_recurrence_seroprev( R0=R0,
                                         foi=foi,
                                         foi_const = foi_const,
                                         time_limit=time_limit,
                                         dim_year=dim_year,
                                         dim_age=dim_age,
                                         p_at=p_at_survey,
                                         P_tot=P_tot_survey,
                                         inc_v=inc_v3d,
                                         vac_eff_arg=vac_eff_arg)


  return(list(i_neg_v_neg_tau3 = seroprev_out$i_neg_v_neg_tau3,
              i_pos_v_neg_tau3 = seroprev_out$i_pos_v_neg_tau3,
              i_neg_v_pos_tau3 = seroprev_out$i_neg_v_pos_tau3,
              i_pos_v_pos_tau3 = seroprev_out$i_pos_v_pos_tau3,
              Ninf_t_survey = seroprev_out$Ninf_t,
              Stot_at_survey = seroprev_out$Stot,
              Nenv_t_survey = seroprev_out$Nenv_t))
}




#' function to calculate aggregated seroprevalence at survey in R0 model
#'
#' @param index_survey index of serological survey
#' @param R0 value of R0 for survey
#' @param foi_const constant force of infection to add to R0 model, defaults to 0
#' @param t0_vac time of vaccination, vector
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param study_years years of sero surveys, vector
#' @param p_at_survey_3d population proportion at survey as an array
#' @param P_tot_survey_2d population total at surveys
#' @param inc_v3d_agg incidence of vacciantion at surveys
#' @param pop_moments_agg aggregated population moments at surveys
#' @param vcfac_arg vaccine coverage factor, defaults to 1
#' @param vac_eff_arg vaccine efficacy
#' @param age_groups age groups in survey
#'
#' @return seroprevalence at age groups
#' @export
R0_calculate_seroprevalence = function(index_survey,
                                       R0,
                                       foi_const = 0,
                                       t0_vac,
                                       dim_year,
                                       dim_age,
                                       study_years,
                                       p_at_survey_3d,
                                       P_tot_survey_2d,
                                       inc_v3d_agg,
                                       age_groups,
                                       pop_moments_agg,
                                       vcfac_arg = 1,
                                       vac_eff_arg) {

  #calculate the seroprevalence at survey
  res = R0_recurrence_seroprev_survey(R0=R0,
                                      foi_const=foi_const,
                                      t0_vac=t0_vac[index_survey],
                                      dim_year=dim_year,
                                      dim_age=dim_age,
                                      p_at_survey= p_at_survey_3d[index_survey,,],
                                      P_tot_survey= P_tot_survey_2d[index_survey,],
                                      inc_v3d=inc_v3d_agg[index_survey,,],
                                      pop_moments=pop_moments_agg[index_survey,],
                                      vac_eff_arg=vac_eff_arg)

  if(vcfac_arg > 1 | vcfac_arg < 0){
    stop(paste0("vcfac is wrong, value: ",vcfac_arg))
  }

  res_out = vcfac_arg*(1 - res$i_neg_v_neg_tau3[study_years[[index_survey]]-1939,]) +
    (1-vcfac_arg)*(1 - res$i_neg_v_neg_tau3[study_years[[index_survey]]-1939,]/
                     (res$i_pos_v_neg_tau3[study_years[[index_survey]]-1939,] +
                        res$i_neg_v_neg_tau3[study_years[[index_survey]]-1939,]))



  ag = findInterval(0:100, age_groups)
  preval1 = res_out * p_at_survey_3d[index_survey, study_years[[index_survey]]-1939,]

  sero_out = aggregate(preval1, by=list(age_group=ag), sum, na.rm=TRUE)$x/
    aggregate(p_at_survey_3d[index_survey, study_years[[index_survey]]-1939,],
              by=list(age_group=ag), sum)$x

  # here if needed, make a condition for vc=1
  sero_out = ifelse( sero_out < 1e-10, 0, sero_out) # rounding may lead to very low (1e-16) negative values


  return(sero_out)
}

#' function to calculate aggregated seroprevalence at survey in R0 model
#'
#' @param dim_adm admin 1 units
#' @param dat data on environment and occurrence
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param p_at_survey_3d population proportion  as an array
#' @param P_tot_survey_2d population total
#' @param inc_v3d incidence of vacciantion
#' @param pop_moments_agg aggregated population moments
#' @param pop_moments_whole popualtion moments
#' @param vcfac_arg vaccine coverage factor, defaults to 1
#' @param vac_eff_arg vaccine efficacy
#' @param granularity how many intervals in R0 values, defaults to "coarse"
#'
#' @return file of R0 related to number of infections over the who time period
#' @export
create_R0_lookup = function(dim_adm,
                            dat,
                            dim_year,
                            dim_age,
                            p_prop_3d,
                            P_tot_2d,
                            inc_v3d,
                            pop_moments_whole,
                            pop_moments_agg,
                            vac_eff_arg = 1,
                            granularity = "coarse") {

  dn_R0 = switch(granularity,
                 "fine" = c(seq(1, 1.099 , by=0.001), seq(1.1, 1.198, b=0.002), seq(1.2, 1.495, by=0.005), seq(1.5, 1.99, by=0.01),
                            seq(2, 2.9, by= 0.05 ), seq(3,3.9, by = 0.1), seq(4,10, by=0.2)),
                 "medium" = c(seq(1,2,by=0.01), seq(2.1,4.9,by=0.1), seq(5,10, by=1) ),
                 "coarse" = c(seq(1,1.9,by=0.1), seq(2,4.5,by=0.5), seq(5,10, by=1) ) )


  R0_lookup = rep(NA, length(dim_adm)*length(dn_R0))
  dim(R0_lookup) = c(length(dim_adm),length(dn_R0))
  colnames(R0_lookup) = dn_R0
  rownames(R0_lookup) = dat$adm0_adm1
  adm_ind = 1:length(dim_adm)

  ptm = proc.time()
  for (r in 1:length(dn_R0)){
    print(paste(r,"/", length(dn_R0), ";   R0=", dn_R0[r]))
    R0_rep = rep(dn_R0[r], length(dim_adm))

    tmp = R0_recurrence_seroprev_whole(adm = adm_ind,
                                       dat = dat,
                                       R0 = R0_rep,
                                       t0_vac_adm = t0_vac_africa,
                                       dim_year = dim_year,
                                       dim_age = dim_age,
                                       p_prop_3d = p_prop_3d,
                                       P_tot_2d = P_tot_2d,
                                       inc_v3d = inc_v3d,
                                       pop_moments_whole = pop_moments_whole,
                                       vac_eff_arg = vac_eff_arg
    )$Ninf_t_province # this is for each year from 1940 to 2050

    R0_lookup[,r] = rowSums(tmp[, which(dim_year==1984):which(dim_year==2019)])
  }
  save(R0_lookup, file="R0_lookup_table.Rdata")

}



### use fun_calc_pdetect_Multi_both and fun_calc_transmission_africa
