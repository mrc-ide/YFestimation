#' Import and bind chains for estimation
#'
#' This functions takes a file path, finds and imports all
#' the .csv files in that directory and then tries to bind them.
#' If specified, it will also remove a burnin period and thin.
#'
#' @param path file path to estimation directory
#' @param burnin burnin iterations to remove
#' @param thin factor to thin chain by
#'
#' @export
#' @return mcmc_out - data frame of estimation
get_chains = function(path,
                      burnin,
                      thin){


  #get files
  temp = list.files(path, pattern = "\\.csv$")
  temp = paste0(path, "/", temp)
  temp = temp[order(file.info(temp)$mtime)]
  temp = temp[file.info(temp)$size>5000]


  #read
  l=lapply(temp,read.csv)
  for(i in 1:length(l)){l[[i]] = l[[i]][,-1]}

  #bind
  mcmc_out=data.table::rbindlist( l )

  #remove burnin
  if(burnin>nrow(mcmc_out)){stop("burnin exceeds chain length")}
  mcmc_out = mcmc_out[burnin:nrow(mcmc_out),]

  #remove NA rows
  mcmc_out = mcmc_out %>% tidyr::drop_na()

  #thin
  mcmc_out=mcmc_out[seq(1,nrow(mcmc_out),thin),]

  #make data frame
  mcmc_out = as.data.frame(mcmc_out)

  return(mcmc_out)
}





#' Plot the prior and posterior deinsities
#'
#' Function to plot the prior and posterior densities for parameters.
#' Note that the prior distributions are hard coded here.
#'
#' @param mcmc_out chain
#' @param prior_col colour of prior distributions, defaults to "red"
#' @param glm_or_sero whether plotting glm parameters or shared serology, either "GLM" or "SERO".
#'
#' @export
#' @return plot of prior and posterior distributions
plot_prior_post = function(mcmc_out,
                           prior_col = "red",
                           glm_or_sero ){

  if(glm_or_sero == "GLM"){
    #set up output plot
    par(mfrow=c(3,2), mar = c(4,4,2,1) + 0.1)

    #get the variables that are country factors
    jj = grep("^log.adm05", names(mcmc_out))
    sd.prior = 2 #hardcoded- check matches GLMprior

    p=as.data.frame(mcmc_out)[,names(mcmc_out)[jj]]
    t = as.vector(as.matrix(p) )

    #plot for country factors
    plot(density(t),
         ylim=c(0,1), main="", xlab = "Country factor")
    vec = seq(min(t)-100,max(t)+100, by = (max(t)-min(t))/1000)
    polygon(vec, dnorm(vec, mean = 0, sd = sd.prior), col = prior_col, border = prior_col)
    polygon(density(t), col = rgb(0,0,0,0.5))


    #plot for everything else
    kk = grep("^log.adm05", names(mcmc_out), invert = TRUE)[1:4]
    for (i in kk){
      q = as.numeric(as.data.frame(mcmc_out)[,names(mcmc_out)[i]])

      plot(density(q), main = "", xlab=names(mcmc_out)[i])
      vec = seq(min(q)-100, max(q)+100, by = (max(q)-min(q))/1000)
      polygon(vec,dnorm(vec, mean = 0, sd = 30), col = prior_col, border = prior_col)
      polygon(density(q), col = rgb(0,0,0,0.5))
    }

  }

  if(glm_or_sero == "SERO"){

    par(mfrow=c(1,2))

    #plot vaccine efficacy
    plot(density(exp(mcmc_out$vac_eff)), xlim = c(min(exp(mcmc_out$vac_eff)),1),
         main = "", xlab="Vaccine efficacy")

    vec = seq(min(exp(mcmc_out$vac_eff)), 1.2, length.out = 1000)
    polygon(vec, dnorm(vec, mean = 0.975, sd = 0.05), col = prior_col, border = prior_col)

    polygon(density(exp(mcmc_out$vac_eff)), col = rgb(0,0,0,0.5))


    #plot vaccine coverage factor for CMRs
    plot(density(exp(mcmc_out$vc_factor_CMRs)), xlim = c(0,1),
         main = "", xlab="Vaccine factor CMRs")

    vec = seq(-0.01, 1.01, length.out = 1000)
    polygon(vec,dunif(vec), col = prior_col, border = prior_col)

    polygon(density(exp(mcmc_out$vc_factor_CMRs)), col = rgb(0,0,0,0.5))
  }
}









#' Internal mean function
#'
#' @param dat vector to be summarised
#' @param idx index
#'
#' @return mean
#'
mean_fun = function(dat, idx) {
  mean(dat[idx], na.rm = TRUE)
}



#' Calculate Bayes factors and bootstrap for uncertainty
#'
#' @param mcmc_out chains from serology
#' @param col colour for histogram, defaults to "orchid"
#' @param plotout whether to plot or not- defaults to TRUE
#'
#' @export
#' @return histogram of model evidence, proportion of iterations
#'         spent on R0 model and log_Bayes_R0_Foi
#'
calculate_bootstrap_Bayes = function(mcmc_out,
                                     col = "orchid",
                                     plotout = TRUE){

  bootobj = boot::boot(data = mcmc_out$model_chain,
                       statistic = mean_fun, R = 1e3)

  hist(bootobj$t, col = col, xlab = "Proportion of time spent on R0 model", main = NULL)

  #Bayes factors
  prob_Foi=  mcmc_out$prob_Foi[1]
  prob_R0 =  log(1-exp(prob_Foi))

  proportion_R0 = quantile(bootobj$t, probs = c(0.025,0.5,0.975))
  log_Bayes_R0_FOI = log(proportion_R0) + prob_Foi - log(1-proportion_R0) -prob_R0

  return(list(proportion_R0 = proportion_R0, log_Bayes_R0_FOI = log_Bayes_R0_FOI))
}









#' Internal function to get hpd
#'
#' @param mcmc_out_chunk section of mcmc_out
#' @export
#' @return hpd of chunk
hpd_out = function(mcmc_out_chunk,
                   log_scale = TRUE){

  mcmc_out1 = mcmcplots::convert.mcmc.list(mcmc_out_chunk)
  hpd_tmp = coda::HPDinterval(mcmc_out1)[[1]]

  if(is.vector(mcmc_out_chunk)){
    hpd_tmp = cbind(hpd_tmp, median(mcmc_out_chunk, na.rm = TRUE))
  } else {
    hpd_tmp = cbind(hpd_tmp, apply(mcmc_out_chunk, 2, median, na.rm = TRUE))
  }

  if(log_scale){
    hpd_tmp = exp(hpd_tmp[, c(1,3,2)])
  } else {
    hpd_tmp = hpd_tmp[, c(1,3,2)]
  }

  return(hpd_tmp)
}


#' Calculating the High probability density regions for all parameters
#'
#' @param mcmc_out chains from production space serology estimation
#' @param glm_mcmc_out chains from glm estimation
#'
#' @return hpd for R0 and Foi models with glm parameters
#' @export
calculate_hpd = function(mcmc_out,
                         glm_mcmc_out){

  if(!anyNA(mcmc_out)){
    ### hpd for serology ###
    #full uncertainty for shared parameters
    hpd_sero = hpd_out(mcmc_out[,grep("^vac", names(mcmc_out))])

    #filter by model type
    mcmc_out_r = mcmc_out %>% filter(model_chain ==1)
    mcmc_out_f = mcmc_out %>% filter(model_chain ==0)

    #hpd by model type
    hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out_f[,grep("^Foi", names(mcmc_out))]) )
    hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out_r[,grep("^R0", names(mcmc_out))])  )

    #last shared parameter
    hpd_sero = rbind(hpd_sero, hpd_out(mcmc_out[,grep("^vc", names(mcmc_out))]))
  }

  if(!anyNA(glm_mcmc_out)){
    ### hpd for glm ###
    hpd_glm = hpd_out( subset( glm_mcmc_out, select = -c(posteriorProb, acceptRate)),
                       log_scale = FALSE)
  }

  if(!anyNA(mcmc_out) & !anyNA(glm_mcmc_out)){

    rownames(hpd_sero)[c(1, nrow(hpd_sero))] = c("vac_eff", "vc_factor_CMRs")

    ### bind all together for each model ###
    R0_param = rbind(hpd_sero[grep("^vac", rownames(hpd_sero)),],
                     hpd_glm,
                     hpd_sero[ grep("^R0", rownames(hpd_sero) ),],
                     hpd_sero[grep("^vc", rownames(hpd_sero)),])

    Foi_param = rbind(hpd_sero[grep("^vac", rownames(hpd_sero)),],
                      hpd_glm,
                      hpd_sero[ grep("^Foi", rownames(hpd_sero) ),],
                      hpd_sero[grep("^vc", rownames(hpd_sero)),])

    rownames(R0_param)[c(1,nrow(R0_param))] =
      rownames(Foi_param)[c(1,nrow(Foi_param))] = c("vac_eff", "vc_factor_CMRs")



    out = list(Foi_param = Foi_param, R0_param = R0_param, hpd_sero = hpd_sero, hpd_glm)
  } else if(anyNA(glm_mcmc_out)){
    rownames(hpd_sero)[c(1, nrow(hpd_sero))] = c("vac_eff", "vc_factor_CMRs")
    out = list(hpd_sero = hpd_sero)
  } else {
    out = list(hpd_glm = hpd_glm)
  }

  return(out)
}





#'calculate predicitions from glm
#'
#'@param coefs cofficients from glm
#'@param newdata environmental data in
#'@param type of link defaults to "response"
#'@param varsin which parameters to include
#'
#'@export
#'@return predictions
fun_calcPred = function(coefs,
                        newdata,
                        type = "response",
                        varsin = NA) {

  if(!(type %in% c("link","response"))) {
    stop("fun_calcPred: invalid type of predictions specified_\n")
  }
  if(is.na(varsin[1])) {
    varsin = 1:length(coefs)
  }
  eta = newdata[,varsin] %*% coefs[varsin]
  if(type=="link") {
    preds = eta #X*beta
  } else if(type=="response") {
    preds = 1-exp(-exp(eta)) # q = 1 - exp (- exp(X*beta))
  }
  return(preds)
}





#' plot a map of glm output alongside map of data
#'
#' @param shp0 shapefile for Africa countries
#' @param shp1 shapefile for Africa admin 1 units
#' @param c34 list of countries
#' @param dat occurrence and environmental data
#' @param Est_beta estimated coefficients
#' @param x covariates
#' @param colours what colours to use- length = 1000
#'
#' @export
plot_glm_map = function(shp0,
                        shp1,
                        c34,
                        dat,
                        Est_beta,
                        x,
                        colours){

  par(mfrow=c(1,2))

  ### data ###
  plot(shp0, xlim=c(-15,45),ylim=c(-20,35))
  mm0 = match(shp0$ISO,c34) #
  plot(shp0[is.na(mm0),],col="black",add=TRUE)

  pres= dat$adm0_adm1[dat$cas.or.out>0]

  mm1<-match(pres, shp1$adm0_adm1)


  plot(shp1[mm1,], col=colours[750], add=TRUE)
  plot(shp0,lwd=2, add=TRUE)
  plot(shp1,lwd=1, add=TRUE)

  ### model ###
  glmpreds_tmp = fun_calcPred( Est_beta,x,type="response")

  mybreaks = seq(0, 1.0001, length.out=1001)
  mycols =  colours
  mm = match(shp1$adm0_adm1, dat$adm0_adm1)
  vcols = findInterval(glmpreds_tmp, mybreaks)


  plot(shp0, xlim=c(-15,45), ylim=c(-20,35))
  mm0 = match(shp0$ISO,c34) #
  plot(shp0[is.na(mm0),], col="black",add=TRUE)
  plot(shp1[!is.na(mm),], col=mycols[vcols], xlim=c(-15,45), ylim=c(-20,30) , lty=0, add=TRUE)
  plot(shp0, xlim=c(-15,45),ylim=c(-20,35), add=TRUE)
  fields::image.plot(legend.only=TRUE, breaks=mybreaks, col=mycols, zlim=c(0,1), horizontal = TRUE)

}



#' plot the seroprevalence predictions with data and credible intervals
#'
#' @param seroout processed serology data
#' @param hpd_sero high probability density region of parameter estimates for product space estimation of seroprevalence models
#' @param pop_agg3d population aggregated in array form
#' @param dim_year years of interest
#' @param dim_age ages of interest
#' @param inc_v3d incidence of vaccination in array form
#' @param pop_moments_agg aggregated population moments
#' @param vc_agg3d vaccination coverage aggregated in array form
#'
#' @export
#' @return figure of all serological surveys and predictions
plot_sero_predictions = function(seroout,
                                 hpd_sero,
                                 pop_agg3d,
                                 dim_year,
                                 dim_age,
                                 inc_v3d,
                                 pop_moments_agg,
                                 vc_agg3d){



  p = create_pop_at_survey(pop_agg3d,
                           seroout$sero_studies,
                           dim_year)

  p_at_survey_3d=p$p_at_survey_3d
  P_tot_survey_2d=p$P_tot_survey_2d


  for (i in 1:3){

    #set parameters
    param =  hpd_sero[,i]
    v_e = hpd_sero[1, i]

    seroprev_predict_survey_r = seroprev_predict_survey_f =  NULL

    for (index_survey in 1:seroout$no_sero_surveys){

      vcfac = ifelse(!is.na(seroout$vc_factor[index_survey]),
                     seroout$vc_factor[index_survey],
                     param[grep("vc_fac", names( param))])

      # SEROPREVALENCE R0 #
      R0surv = param[grep("R0", names( param))][index_survey]
      res = R0_recurrence_seroprev_survey(R0=R0surv,
                                          foi_const = 0,
                                          t0_vac=seroout$t0_vac[index_survey],
                                          dim_year=dim_year,
                                          dim_age=dim_age,
                                          p_at_survey= p_at_survey_3d[index_survey,,],
                                          P_tot_survey= P_tot_survey_2d[index_survey,],
                                          inc_v3d=inc_v3d_agg[index_survey,,],
                                          pop_moments=pop_moments_agg[index_survey,],
                                          vac_eff_arg = v_e)


      res_out = vcfac*(1 - res$i_neg_v_neg_tau3[seroout$study_years[[index_survey]]-1939,]) +
        (1-vcfac)*(1 - res$i_neg_v_neg_tau3[seroout$study_years[[index_survey]]-1939,]/
                     (res$i_pos_v_neg_tau3[seroout$study_years[[index_survey]]-1939,] +
                        res$i_neg_v_neg_tau3[seroout$study_years[[index_survey]]-1939,]))

      seroprev_predict_survey_r = rbind(seroprev_predict_survey_r,  res_out)

      # SEROPREVALENCE FOI #
      Foisurv = param[grep("^Foi", names( param))][index_survey]

      vc_agg_ag = vc_agg3d[index_survey, paste0(seroout$study_years[index_survey]),]
      vc=vcfac*v_e*vc_agg_ag

      pop = pop_agg3d[index_survey,paste0(seroout$study_years[index_survey]),]

      sero1 = c(as.matrix((1-exp(-Foisurv*0:100))*pop))
      sero_ag = aggregate(sero1, by=list(age_group=0:100), sum)
      sero_ag$x = sero_ag$x/aggregate(pop, by=list(age_group=0:100), sum)$x

      res_out = 1 - (1-sero_ag$x[sero_ag$age_group>0])*(1-vc)

      seroprev_predict_survey_f = rbind(seroprev_predict_survey_f,  res_out)

    } # end of for(index_survey)

    seroprev_predict_survey_print_r = data.frame(sero_studies=rep(seroout$sero_studies, times = 1),
                                                 seroprev_predict_survey_r)
    seroprev_predict_survey_print_f = data.frame(sero_studies=rep(seroout$sero_studies, times = 1),
                                                 seroprev_predict_survey_f)


    if (i==1) {
      r_sero_predictions_lo=seroprev_predict_survey_print_r[1:seroout$no_sero_surveys,2:102]
      f_sero_predictions_lo=seroprev_predict_survey_print_f[1:seroout$no_sero_surveys,2:102]
    }
    if (i==2) {
      r_sero_predictions=seroprev_predict_survey_print_r[1:seroout$no_sero_surveys,2:102]
      f_sero_predictions=seroprev_predict_survey_print_f[1:seroout$no_sero_surveys,2:102]
    }
    if (i==3) {
      r_sero_predictions_hi=seroprev_predict_survey_print_r[1:seroout$no_sero_surveys,2:102]
      f_sero_predictions_hi=seroprev_predict_survey_print_f[1:seroout$no_sero_surveys,2:102]
    }


  } #end of for(i)


  ### PLOT ###
  par(mfrow=c(4,10),  par(mar = rep(2, 4)))
  for (i in 1:seroout$no_sero_surveys){
    #get upper bound for plot axes
    maxsero = max(as.numeric(f_sero_predictions_hi[i,1:86]),
                  na.rm = TRUE)

    #pull out key parts of data for plotting
    age_points = ( seroout$survey_dat[[i]]$age_min+seroout$survey_dat[[i]]$age_max )/2
    seroprevalence = seroout$survey_dat[[i]]$positives / seroout$survey_dat[[i]]$samples

    #plot the data
    plot( age_points,
          seroprevalence,
          main = paste0((seroout$sero_studies[i])),
          cex = 1,
          ylab = "Seroprevalence", pch=19, col = alpha("black",0.5),
          xlim=c(0,85), ylim=c(0,maxsero+0.1 ), xlab = "Age")

    #add the HPD of predictions
    polygon(c(85:0,0:85), c(rev(as.numeric(f_sero_predictions_hi[i,1:86])),
                            (as.numeric(f_sero_predictions_lo[i,1:86]))) ,
            border=NA, col = rgb(0, 0, 1,0.2) )
    polygon(c(85:0,0:85), c(rev(as.numeric(r_sero_predictions_hi[i,1:86])),
                            (as.numeric(r_sero_predictions_lo[i,1:86]))) ,
            border=NA, col = rgb(1, 0, 0,0.2) )

    #add the median predictions
    lines(0:85,as.numeric(r_sero_predictions[i,1:86]),
          type="l",col = "red")
    lines(0:85,as.numeric(f_sero_predictions[i,1:86]),
          type="l", main = paste0((seroout$sero_studies[i])),
          ylab = "Seroprevalence", col = "blue" )

    ## add binomial confidence for the data
    conf_int = Hmisc::binconf(seroout$survey_dat[[i]]$positives,
                              seroout$survey_dat[[i]]$samples)

    arrows(age_points,
           conf_int[,2],
           age_points,
           conf_int[,3], length=0.05, angle=90, code=3)

    points(age_points, conf_int[,1])
  }
}


#' calculate and plot the transmission intensity across the African endemic region
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
#'@param shp1  shapefiles at admin 1
#'@param shp0 shape files at country level
#'@param Countries list of countries
#'@param colours vector of length 100
#'
#'@export
#'@return map of transmission intensity
plot_transmission_intensity = function(x,
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
                                       polydeg,
                                       R0_lookup,
                                       model_type,
                                       shp1,
                                       shp0,
                                       Countries,
                                       colours){


  runs = fun_calc_transmission_Africa(x,
                                      ii ,
                                      seroout ,
                                      params ,
                                      dat ,
                                      t0_vac_africa ,
                                      dim_year ,
                                      dim_age ,
                                      p_prop_3d ,
                                      P_tot_2d ,
                                      inc_v3d ,
                                      pop1,
                                      vc2d,
                                      varsin_nc,
                                      polydeg,
                                      R0_lookup,
                                      model_type)

  runs = switch(model_type,
                "R0" = runs -1,
                "Foi" = runs)


  mybreaks= seq(min(log10(runs)), max(log10(runs))+0.01, length.out=101)
  mm = match(shp1$adm0_adm1,dat$adm0_adm1)
  vcols = findInterval(log10(runs),mybreaks)

  plot(shp0, xlim=c(-15,45),ylim=c(-20,30))
  mm0 = match(shp0$ISO,Countries$c34) #
  plot(shp0[is.na(mm0),],col="grey70",add=T)
  plot(shp1[!is.na(mm),],col=colours[vcols], xlim=c(-15,45),ylim=c(-20,30) , lty=0, add=T)

  plot(shp0, lwd=2, add=T)


  if(model_type == "R0"){
    image.plot(legend.only=TRUE, breaks = mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
               axis.args = list(at = c(-5:0, log10(9)), labels =c("1.00001", "1.0001", "1.001", "1.01", "1.1","2","10+"), las =2),
               legend.mar = 3.5)
  } else {
    image.plot(legend.only=TRUE, breaks=mybreaks, col=colours, zlim=c(0,1), horizontal = TRUE,
               axis.args = list(at = c(-5:1), labels =c("1e-05","1e-04", "0.001", "0.01", "0.1","1","10"), las =2),
               legend.mar = 3.5)
  }
}
