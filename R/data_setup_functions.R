# function file to work with model_reestimation

#' load environmental data from path
#'
#' @param filepath file path to data csv
#' @param c34 list of countries
#' @param dat_full data if already imported, defaults to NA
#' @param delete_surv_AGO remove AGO from surveillance quality, defaults to TRUE
#'
#'
#' @return depi and dat, dependant variable index and data
#' @export
launch_env_dat = function(filepath, c34, dat_full = NA, delete_surv_AGO = TRUE) {

  if(!is.na(filepath)){
    dat_full = read.csv(filepath, stringsAsFactors=FALSE)
  } else if(is.na(dat_full)){
    stop("no valid filepath or data entered")
  }
  #which presence/abence outcome to consider. Alternatives : "cases" # "outbreaks" #
  depvar = "cas.or.out"
  depi = match(depvar, names(dat_full)) # column number for the chosen outcome

  # excluding LC1,3,5,15 as they only ever cover <5% -
  ex_i = match(c("LC1","LC3","LC5","LC15"), names(dat_full))
  dat_full = dat_full[,-ex_i]

  # adding a categorical "dominant land cover" variable:
  LC_i = grep("LC",names(dat_full))
  LC_dom = names(dat_full)[LC_i][ apply(dat_full[,LC_i],1, which.max )]
  dat_full = cbind(dat_full, LC_dom) # 37 potential covariates/levels

  #setting NA to 0
  dat_full$surv.qual.adm0[is.na(dat_full$surv.qual.adm0)] =
    dat_full$surv.qual.adm1[is.na(dat_full$surv.qual.adm1)] = 0

  # setting the surv.qual.adm0 for AGO to 0 (very few cases reported):
  if(delete_surv_AGO){
    dat_full$surv.qual.adm0[dat_full$adm0=="AGO"] = 0
  }

  #adding log(surv_qual_adm0)
  dat_full = dplyr::mutate(dat_full,
                           log.surv.qual.adm0 = ifelse(is.infinite(log(surv.qual.adm0)),
                                                       0, log(surv.qual.adm0)))

  # a same categorical variable for all countries within the YFSD
  dat_full = dplyr::mutate(dat_full,
                           adm05 = ifelse(surv.qual.adm0>0,
                                          "AFR", adm0))

  dat =  dplyr::filter(dat_full, adm0 %in% c34)


  v1 = apply(dat,2,var, na.rm = TRUE)
  for(i in 8:(ncol(dat))) {

    if(!is.factor(dat[,i]) & !is.character(dat[,i])) {

      dat[,i] = dat[,i]/sqrt(v1[i])
      dat_full[,i] = dat_full[,i]/sqrt(v1[i])

    }
    # Explanation:
    # If we fit the model on dat and if we want to project estimates on dat_full, variable from dat_full
    # need to be expressed on the same scale than those from dat , thus we normalize dat_full relatively to dat

  }

  dat = dat[,names(dat)!="surv_qual_adm0"]
  depi = match(depvar,names(dat))

  # I do that because pop is ordered that way and we need to match both pop and dat
  dat = dat[ order(dat$adm0_adm1), ]

  return(list(depi = depi, dat = dat))
}




#' fit the glm
#'
#' @param dat data with wenvironmental and occurrence information
#' @param depi dependant variables
#' @param models list of models- only use first one
#'
#'
#' @return coefficients of glm
#' @export
fit_glm = function(dat, depi, models) {

  fm_best = models[1] #take first model

  bm = glm(as.formula(fm_best), data=dat, family=binomial(link="cloglog"))

  # setting up the evaluation of the likelihood:
  beta = coefficients(bm)

  vl=NULL
  for(i in 1:ncol(dat)) {
    if(length(grep(names(dat)[i], names(beta)))>0) {
      vl = c(vl,i) # select the variables used in the GLM
    }
  }


  x = cbind(Intercept=1,dat[,vl])
  j_expand = !sapply(1:ncol(x), function(i) is.numeric(x[,i]) & !is.factor(x[,i]))

  x_num = NULL
  for(j in 1:ncol(x)) {   # create indicative variables for adm05

    if(j_expand[j]) {
      tab = table(1:nrow(x),x[,j])
      colnames(tab) = paste(colnames(x)[j],colnames(tab),sep="")
      x_num = cbind(x_num,tab[,-1])
    } else {
      x_num = cbind(x_num,x[,j])
      colnames(x_num)[ncol(x_num)] = colnames(x)[j]
    }
  }
  x = x_num # covariate matrix
  rm(x_num)
  y = dat[,depi]  # dependant variable vector

  beta0 = coefficients(bm)
  names(beta0)[names(beta0)=="(Intercept)"] = "Intercept"

  mm = match(names(beta0),colnames(x))
  x = x[,mm]

  beta0 = beta0[match(colnames(x),names(beta0))]

  return(list(beta0=beta0, x=x, y=y))

}






#' old popualtion data import function
#'
#' @param path path to data
#' @param c_country names of countries
#' @param adm admin level
#'
#'
#' @return population by age and year
#' @export
import_pop_data_LS2014 = function(path, c_country, adm) {

  pop_year_age = NULL

  for(adm0 in c_country){
    tmp = read.csv(paste0(path,"Population", "/from_LandScan2014/", adm, "/pop_year_age_",
                          adm, "_", adm0, "_LandScan2014.csv"),
                   header=TRUE)
    pop_year_age = rbind(pop_year_age, tmp)
  }

  pop_year_age[,adm] = paste(pop_year_age[,"adm0"], pop_year_age[,adm], sep="_")
  colnames(pop_year_age)[colnames(pop_year_age)==adm] = paste("adm0", adm, sep = "_")

  return(pop_year_age)
}





#' population data import function
#'
#' @param path path to data
#' @param c_country names of countries
#' @param adm admin level
#'
#'
#' @return population by age and year
#' @export
import_pop_data_LS2017 = function(path, c_country, adm) {

  pop_year_age = read.csv(paste0(path,"Population/", "from_Landscan2014/",
                                 "population_by_adm1_year_age_1950-2100_LandScan2014_gadm2_countries.csv"),
                          header=TRUE)

  pop_year_age = pop_year_age[pop_year_age$adm0 %in% c_country,]

  pop_year_age[,adm] = paste(pop_year_age[,"adm0"], pop_year_age[,adm], sep="_")
  colnames(pop_year_age)[colnames(pop_year_age)==adm] = paste("adm0", adm, sep = "_")

  return(pop_year_age)
}





#' internal fucntion to add early years to population data
#'
#' @param pop2d 2-dimensional population data
#'
#'
#' @return 2-dimensional population data with early years
#' @export
#'
add_1940_1950 = function(pop2d){
  if(min(pop2d$year) != 1950)  stop("pop2d should start in 1950")

  for(y in 1949:1940) {
    pop1_early = pop2d[pop2d$year==y+1,]
    pop1_early$year = y
    pop2d = rbind(pop1_early,pop2d)
  }

  return(pop2d)
}






#' function to reorganise population data
#'
#' @param pop2d 2-dimensional population data
#' @param admin level
#'
#'
#' @return pop3d, P_tot_2d, p_prop_3d ie. population as an array,
#' the total population in a matrix and the population proportion in an array
Make_pop3d_Ptot2d_Pprop3d = function(pop2d, adm){

  adm0_adm = paste0("adm0_", adm)

  dim_adm  = unique(pop2d[,adm0_adm])
  dim_year = as.numeric(unique(pop2d$year))
  dim_age  = names(pop2d)[4:length(names(pop2d))]

  length_adm=length(dim_adm);   length_year=length(dim_year);   length_age=length(dim_age)

  #pop3d =3d array- dim [adm, year, age]
  pop3d=rep(NA, nrow(pop2d)*length_age)
  dim(pop3d)=c(length_adm, length_year, length_age)

  for(ageIndex in 1:length_age) {
    for(yearIndex in min(dim_year):max(dim_year)) {
      mm = match( pop2d[,adm0_adm][pop2d$year == yearIndex], dim_adm )
      pop3d[mm, yearIndex-min(dim_year)+1, ageIndex] = pop2d[pop2d$year == yearIndex, ageIndex+3]
    }
  }

  #P_tot2d = sum(pop3d(adm,year))
  P_tot_2d = rep(NA, length_adm*length_year);   dim(P_tot_2d)=c(length_adm,length_year)

  #p_prop_3d = pop3d/P_tot_2d
  p_prop_3d=rep(NA, length_adm*length_year*length_age);   dim(p_prop_3d)=dim(pop3d)

  for(admIndex in 1:length_adm){
    for (yearIndex in 1:length_year){

      P_tot_2d[admIndex, yearIndex] = sum(pop3d[admIndex, yearIndex,], na.rm=TRUE)

      for (ageIndex in 1:length_age){
        p_prop_3d[admIndex,yearIndex,ageIndex] = pop3d[admIndex, yearIndex, ageIndex]/P_tot_2d[admIndex, yearIndex]
      }
    }
  }

  #names
  dimnames(p_prop_3d)[[1]] = dimnames(pop3d)[[1]] = rownames(P_tot_2d) = dim_adm
  dimnames(p_prop_3d)[[2]] = dimnames(pop3d)[[2]] = colnames(P_tot_2d) = dim_year
  dimnames(p_prop_3d)[[3]] = dimnames(pop3d)[[3]] = dim_age

  return(list(pop3d = pop3d, P_tot_2d = P_tot_2d, p_prop_3d = p_prop_3d))

}


#' function to import and get pop data in array
#'
#' @param path path to data
#' @param c_country list of countries
#' @param dat environmental and occurrence data
#'
#'
#' @export
#' @return pop1, pop3d, P_tot_2d, p_prop_3d
get_pop_data_3d = function(path, c_country, dat) {

  adm="adm1"
  adm0_adm = paste0("adm0_", adm)

  pop1 = import_pop_data_LS2017(path, c_country, adm)

  #restrict to lines in dat THIS IS CHANGED FOR 2017!
  pop1 = pop1[pop1[,adm0_adm] %in% dat[,adm0_adm],]

  #make sure data starts from 1940
  pop1 = add_1940_1950(pop1)

  out = Make_pop3d_Ptot2d_Pprop3d(pop1, adm)

  return(list(pop1 = pop1, pop3d = out$pop3d, P_tot_2d = out$P_tot_2d, p_prop_3d = out$p_prop_3d))
}



#' get vaccine data into an array
#'
#' @param vc2d 2-dimensional vaccination data
#' @param adm defunct
#'
#'
#' @return vc3d
#' @export
transform_into_vc3d = function(vc2d, adm = NA){

  vc_tmp = tidyr::gather(vc2d, age, coverage, -c(year, adm0_adm1))

  vc3d = array(data = vc_tmp$coverage,
               dim=c(length(unique(vc_tmp$adm0_adm1)),
                     length(unique(vc_tmp$year)),
                     length(unique(vc_tmp$age))),
               dimnames=list(unique(vc_tmp$adm0_adm1),
                             unique(vc_tmp$year),
                             unique(vc_tmp$age))
  )

  return(vc3d)
}



#' calculate the first instance of vaccine coverage
#'
#' @param vc3d vaccination data in array form
#'
#'
#'
#' @return inc_v3d
#' @export
calc_incidence_vac_general = function(vc3d){

  inc_v3d = vc3d*NA

  for (admIndex in 1:dim(vc3d)[1] ) { #

    inc_v3d[admIndex,1,]= 0 # no vaccination in 1940

    for(yearIndex in 2:(dim(vc3d)[2]-1)){

      inc_v3d[admIndex,yearIndex,1] = vc3d[admIndex, yearIndex+1,1] # for a0, incidence = coverage in the a1 age group the next year

      for(age in 2:(dim(vc3d)[3]-1)) {

        if(!is.na(vc3d[admIndex, yearIndex, age]) & vc3d[admIndex, yearIndex, age]==1){ #if vc is already =1, incidence =0

          inc_v3d[admIndex,yearIndex,age]=0

        } else {

          inc_v3d[admIndex, yearIndex, age] =  ( vc3d[admIndex, yearIndex+1, age+1] - vc3d[admIndex, yearIndex, age] ) /
            (1 - vc3d[admIndex,yearIndex,age] ) # among those susceptible at y-1
          inc_v3d[admIndex, yearIndex, age] = ifelse(inc_v3d[admIndex,yearIndex,age]<10e-15,
                                                     0, inc_v3d[admIndex,yearIndex,age])
        }# with rounding, some values become negative
      }
      inc_v3d[admIndex, yearIndex, dim(vc3d)[3]]=0 # incidence vaccination = 0 for age=100
    }

    inc_v3d[admIndex,dim(vc3d)[2],]  = 0# incidence vaccination = 0 for year=2050
    inc_v3d[admIndex,dim(vc3d)[2],1] = vc3d[admIndex,dim(vc3d)[2],1] # except among new borns
  }

  return(inc_v3d)
}



#' get time of first vaccination
#'
#' @param vc3d vaccination data in array form
#'
#'
#' @return t0_vac_africa
#' @export
calc_t0_vac_africa = function(vc3d){

  t0_vac_africa = rep(NA, dim(vc3d)[1])

  for (i in 1:dim(vc3d)[1]){

    year_i  = 1940
    sum_vac = 0

    while(sum_vac == 0 & year_i <= 2050) {

      sum_vac = sum( vc3d[i, year_i-1940+1,], na.rm=TRUE)
      year_i = year_i + 1

    }

    t0_vac_africa[i]=year_i-2

  }

  names(t0_vac_africa) = dimnames(vc3d)[[1]]

  return(t0_vac_africa)
}




#' calculate population moments
#'
#' @param p_prop proportion of population
#' @param t0_vac_ intial vaccination time in area of interest
#' @param dim_adm admins of interest
#' @param dim_year years of interest
#' @param dim_age ages of interest
#'
#' @return pop_moments_whole
#' @export
calc_pop_moments = function(p_prop,
                            t0_vac_,
                            dim_adm,
                            dim_year,
                            dim_age) {

  ageVec=c(0:100)

  minYear=1983 #min year of data
  maxYear=2019 #last year of data
  print("Remember, last year of data is hardcoded into calc_pop_moments")

  maxDuration=maxYear-minYear + 1 #duration of data for second for loop

  length_adm=length(dim_adm)

  pop_moments_whole = rep(NA, length_adm*6)
  dim(pop_moments_whole)=c(length_adm,6)

  for (admIndex in 1:length_adm){
    if (t0_vac_[admIndex]<maxYear){

      vec = which(dim_year==1940):max(which(dim_year==t0_vac_[admIndex]), 1)

      pop_moments_whole_tmp = rep(NA, length(vec)*6)
      dim(pop_moments_whole_tmp) = c(length(vec), 6)

      for(yearIndex in vec){
        for(i in 1:6){
          pop_moments_whole_tmp[yearIndex-min(vec)+1, i] = sum(p_prop[admIndex, yearIndex,]*ageVec^(i-1),
                                                               na.rm = TRUE)
        }
      }

    } else if (t0_vac_[admIndex]>=maxYear){

      pop_moments_whole_tmp = rep(NA, maxDuration*6)
      dim(pop_moments_whole_tmp) = c(maxDuration, 6)
      vec = which(dim_year==minYear):which(dim_year==maxYear)

      for(yearIndex in vec){
        for(i in 1:6){
          pop_moments_whole_tmp[yearIndex-min(vec)+1, i] = sum(p_prop[admIndex, yearIndex,]*ageVec^(i-1),
                                                               na.rm = TRUE)
        }
      }
    }
    pop_moments_whole[admIndex,] = apply(pop_moments_whole_tmp,2,mean)
  }
  return(pop_moments_whole)
}




#' calculate the population moments for aggregated data
#'
#' @param pop_agg3d population aggregated in array
#' @param t0_vac_ initial vaccination time for area of interest
#' @param dim_year years of interest
#' @param study_years years of serological studies
#'
#'
#' @return pop_moments_agg
#' @export
calc_pop_moments_agg = function(pop_agg3d,
                                t0_vac_,
                                dim_year,
                                study_years) {

  ageVec=c(0:100)# a matrix of age with the same dimension that pop_agg for 1 fixed year

  minYear=1983 #min year of data
  maxYear=2019 #last year of data
  print("Remember, last year of data is hardcoded into calc_pop_moments_agg")

  maxDuration=maxYear-minYear + 1 #duration of data for second for loop
  n_serosurveys=length(study_years)

  pop_moments_agg = rep(NA, n_serosurveys*6)
  dim(pop_moments_agg)=c(n_serosurveys,6)

  for (index_survey in 1:n_serosurveys){

    if (t0_vac_[index_survey]<study_years[[index_survey]]){

      vec = which(dim_year==1940):which(dim_year==t0_vac_[index_survey])
      pop_moments_agg_tmp = rep(NA, length(vec)*6)
      dim(pop_moments_agg_tmp) = c(length(vec), 6)

      for(year in vec){
        p_prop_agg = pop_agg3d[index_survey, year,]/sum(pop_agg3d[index_survey, year,],
                                                        na.rm = TRUE)

        for(i in 1:6){
          pop_moments_agg_tmp[year, i] = sum(p_prop_agg*ageVec^(i-1),
                                             na.rm=TRUE)
        }
      }

    } else if (t0_vac_[index_survey]>=maxYear){

      vec = which(dim_year==minYear):which(dim_year==maxYear)
      pop_moments_agg_tmp = rep(NA, length(vec)*6)
      dim(pop_moments_agg_tmp) = c(length(vec), 6)

      for(year in vec){
        p_prop_agg = pop_agg3d[index_survey, year,]/sum(pop_agg3d[index_survey, year,],
                                                        na.rm = TRUE)
        for(i in 1:6){
          pop_moments_agg_tmp[year-min(vec)+1, i] = sum(p_prop_agg*ageVec^(i-1),
                                                        na.rm = TRUE)
        }
      }
    }
    pop_moments_agg[index_survey,] = apply(pop_moments_agg_tmp,2,mean)
  }
  return(pop_moments_agg)
}





#' aggregating population and vaccination data for serological surveys
#'
#' @param pop1 population data
#' @param vc2d vaccine coverage in 2-d
#' @param sero_studies list of serological studies
#' @param adm1s list of admins for each survey
#'
#'
#' @return pop_agg, vc_agg
#' @export
Make_aggregate_pop_vc = function(pop1,
                                 vc2d,
                                 sero_studies,
                                 adm1s) {

  no_sero_surveys=length(sero_studies)

  vc_tmp = pop_tmp = pop_agg = vc_agg = NULL

  for(yearIndex in 1940:max(pop1$year)) { # if needed replace 2050 by max(unlist(study_years))
    for(surveyIndex in 1:no_sero_surveys) {

      pop_tmp = pop1 %>% filter(year == yearIndex &
                                  adm0_adm1 %in% adm1s[[surveyIndex]])

      vc_tmp = vc2d %>% filter(year == yearIndex &
                                 adm0_adm1 %in% adm1s[[surveyIndex]])

      if (nrow(pop_tmp)==1) {

        pop_tmp %<>% mutate(adm0 = sero_studies[surveyIndex])
        vc_tmp %<>% mutate(adm0 = sero_studies[surveyIndex])

      } else { # need to do some summing/averaging...

        if(anyNA(vc_tmp)){
          pop_tmp[is.na(vc_tmp)] = NA
        }



        vc_tmp2 = colSums(select(pop_tmp, -c(year, adm0_adm1)) * select(vc_tmp, -c(year, adm0_adm1)),na.rm=TRUE)/
          colSums(select(pop_tmp, -c(year, adm0_adm1)),na.rm=TRUE)


        names(vc_tmp2) = paste0("a",0:100)

        vc_tmp = data.frame(adm0 = sero_studies[surveyIndex],
                            adm0_adm1 = NA,
                            year = yearIndex,
                            t(vc_tmp2) )

        pop_tmp = colSums(select(pop_tmp, -c(year, adm0_adm1)))

        names(pop_tmp) = paste0("a",0:100)

        pop_tmp = data.frame(adm0 = sero_studies[surveyIndex],
                             adm0_adm1 = NA,
                             year = yearIndex,
                             t(pop_tmp))

      }
      pop_agg = rbind(pop_agg,pop_tmp)
      vc_agg = rbind(vc_agg,vc_tmp)

    }
  }

  pop_agg[,-(1:3)][is.na(pop_agg[,-(1:3)])] = 0
  vc_agg[,-(1:3)][is.na(vc_agg[,-(1:3)])] = 0

  return(list(pop_agg = pop_agg, vc_agg = vc_agg))
}






#' Aggregate the population and vaccination data into an array
#'
#' @param pop1 population data
#' @param vc2d vaccine coverage in 2-d
#' @param sero_studies list of serological studies
#' @param adm1s list of admins for each survey
#'
#'
#' @return pop_agg3d, vc_agg3d
#' @export
Make_aggregate_pop_vc_3d = function(pop1,
                                    vc2d,
                                    sero_studies,
                                    adm1s){

  agg=Make_aggregate_pop_vc(pop1, vc2d, sero_studies, adm1s)

  pop_agg=agg$pop_agg
  vc_agg=agg$vc_agg

  dim_survey = sero_studies
  dim_year = as.numeric(names(table(pop_agg$year)))
  dim_age = names(pop1 %>% select(-c(adm0_adm1, year)))

  length_survey=length(dim_survey);   length_year=length(dim_year);   length_age=length(dim_age)

  ## pass in 3d
  pop_agg3d=rep(NA, nrow(pop_agg)*length_age)
  dim(pop_agg3d)=c(length_survey, length_year, length_age)
  vc_agg3d=rep(NA, nrow(vc_agg)*length_age)
  dim(vc_agg3d)=c(length_survey, length_year,length_age)

  for(ageIndex in 1:length_age) { # dim_age
    for(yearIndex in min(dim_year):max(dim_year)) {

      mm = match(pop_agg$adm0[pop_agg$year == yearIndex], dim_survey)
      pop_agg3d[mm,yearIndex-min(dim_year)+1,ageIndex] = pop_agg[pop_agg$year==yearIndex, ageIndex+3]
      vc_agg3d[mm,yearIndex-min(dim_year)+1,ageIndex] = vc_agg[vc_agg$year==yearIndex, ageIndex+3]
    }
  }

  dimnames(pop_agg3d)[[1]] = dimnames(vc_agg3d)[[1]] = dim_survey
  dimnames(pop_agg3d)[[2]] = dimnames(vc_agg3d)[[2]] = dim_year
  dimnames(pop_agg3d)[[3]] = dimnames(vc_agg3d)[[3]] = dim_age

  return(list(pop_agg3d = pop_agg3d, vc_agg3d = vc_agg3d))

}





#'get the population at the survey times and places
#'
#'@param pop_agg3d aggregated population
#'@param dim_survey surveys of interest
#'@param dim_year years of interest
#'
#'
#'@return p_at_survey_3d, P_tot_survey_2d
#'@export
create_pop_at_survey = function(pop_agg3d, dim_survey, dim_year) {

  # Create p_at_survey and P_tot, population structure and number, aggregated at survey provinces
  p_at_survey_3d=rep(NA, dim(pop_agg3d)[1]*dim(pop_agg3d)[2]*dim(pop_agg3d)[3])
  dim(p_at_survey_3d)=dim(pop_agg3d)
  dimnames(p_at_survey_3d) =dimnames(pop_agg3d)

  P_tot_survey_2d=  matrix(rep(NA, dim(pop_agg3d)[1]*dim(pop_agg3d)[2]), nrow=length(dim_survey))
  rownames(P_tot_survey_2d) = dim_survey; colnames(P_tot_survey_2d) = dim_year

  for (i in 1:length(dim_survey)){
    for (Y in  1:length(dim_year) ){
      P_tot_survey_2d[i,Y] = sum(pop_agg3d[i,Y,],
                                 na.rm=TRUE)
      p_at_survey_3d[i,Y,] = pop_agg3d[i,Y,]/ P_tot_survey_2d[i,Y]
    }
  }

  return(list(p_at_survey_3d = p_at_survey_3d, P_tot_survey_2d = P_tot_survey_2d))
}




#' get aggregated population and vaccination across observation period
#'
#' @param pop1 population data
#' @param vc2d vaccine coverage data
#'
#' @export
#' @return pop30_agg, pop_vc_moments and vc30_agg, population and vaccination aggregated over observation period
create_pop30_agg_vc30_agg = function(pop1, vc2d){
  #print("Remember, final year of data is hardcoded into create_pop30_agg_vc30_agg")

  ### subset by year ###
  pop30 = dplyr::filter(pop1, year>1983 & year<2019);  pop30 = pop30[order(pop30$adm0_adm1),]
  vc30 = dplyr::filter(vc2d, year>1983 & year<2019)

  pop30[is.na(pop30) | is.na(vc30)] = 0
  vc30[is.na(vc30)] = 0

  ### aggregate first age group ###
  pop30_agg = aggregate(pop30[,"a0"],
                        by=list(adm0_adm1=pop30$adm0_adm1),
                        sum) # a0 is the 4th column,

  vc30_agg = aggregate(vc30[,"a0"]*pop30[,"a0"],
                       by=list(adm0_adm1=vc30$adm0_adm1),
                       sum) # number of vaccinated people by age, aggregated over 30y

  vc30_agg$x = vc30_agg$x/pop30_agg$x

  names(pop30_agg)[2] = names(vc30_agg)[2] = names(pop30)["a0"]

  ### then other age groups ###
  for(i in grep("a1", names(pop30))[1]:ncol(pop30)) {

    pop30_agg = cbind(pop30_agg, aggregate(pop30[,i],
                                           by=list(adm0_adm1=pop30$adm0_adm1),
                                           sum)$x)

    vc30_agg = cbind(vc30_agg,
                     aggregate(vc30[,i]*pop30[,i],
                               by=list(adm0_adm1=vc30$adm0_adm1),
                               sum)$x/
                                      aggregate(pop30[,i],
                                                by=list(adm0_adm1=pop30$adm0_adm1),
                                                sum)$x)

    names(pop30_agg)[ncol(pop30_agg)] = names(vc30_agg)[ncol(vc30_agg)] = names(pop30)[i]
  }

  vc30_agg[is.na(vc30_agg)] = 0


  amat = matrix(0:100, nrow=nrow(pop30_agg), ncol=101, byrow=TRUE)
  pop_vc_moments = data.frame(adm0_adm1 = pop30_agg$adm0_adm1,
                              m0 = rowSums(pop30_agg[,-1]*(1-vc30_agg[,-1]),
                                           na.rm=TRUE),
                              m1 = rowSums(pop30_agg[,-1]*(1-vc30_agg[,-1])*amat,
                                           na.rm=TRUE),
                              m2 = rowSums(pop30_agg[,-1]*(1-vc30_agg[,-1])*amat*amat,
                                           na.rm=TRUE),
                              m3 = rowSums(pop30_agg[,-1]*(1-vc30_agg[,-1])*amat*amat*amat,
                                           na.rm=TRUE),
                              m4 = rowSums(pop30_agg[,-1]*(1-vc30_agg[,-1])*amat*amat*amat*amat,
                                           na.rm=TRUE),
                              m5 = rowSums(pop30_agg[,-1]*(1-vc30_agg[,-1])*amat*amat*amat*amat*amat,
                                           na.rm=TRUE))

  return(list(pop30_agg=pop30_agg, vc30_agg=vc30_agg, pop_vc_moments=pop_vc_moments))
}



#' extract key information from serological data
#'
#' @param Serology dataframe of serological data
#' @param adm which gamd- defaults to "gadm2"
#'
#' @return survey_dat, adm1s, sero_studies, study_years, vc_factor, t0_vac, no_sero_studies
#' @export
process_serology = function(Serology, adm = "gadm2"){

  #extracting values as the file is read as a tibble
  sero_studies =  unique(Serology$country_zone)                                      #serology locations
  study_years = unique(Serology[,c("country_zone","year") ] )[ ,2]    #years each survey occurred
  no_sero_surveys = length(sero_studies)                                             #Number of sero surveys in database

  vc_factor = unique(Serology[,c("country_zone","vc_factor") ] )[ ,2]  #whether to account for vaccination in serology or not
  t0_vac = unique(Serology[,c("country_zone","t0_vac") ] )[ ,2]         #first incidence of vaccination, pulled from tibble
  names(t0_vac)=unique(Serology[,c("country_zone","t0_vac") ] )[ ,1]

  adm1s_tibble = unique(Serology[ , c("country_zone", adm)])[ ,2]                   #shape file admin locations in tibble form
  adm1s=list()
  for (surveyIndex in 1:no_sero_surveys) {
    adm1s[surveyIndex] = strsplit(as.character(adm1s_tibble[surveyIndex]),",")
  }

  survey_dat=list()
  for (surveyIndex in 1:no_sero_surveys) {
    survey_dat[[surveyIndex]] = filter(Serology,
                                       country_zone==sero_studies[surveyIndex])[,c("samples","positives","age_min", "age_max")]
  }

  return(list(survey_dat = survey_dat,
              adm1s = adm1s,
              sero_studies = sero_studies,
              study_years = study_years,
              vc_factor = vc_factor,
              t0_vac = t0_vac,
              no_sero_surveys = no_sero_surveys))
}
