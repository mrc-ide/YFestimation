

#' function to calculate aggregated seroprevalence at survey in foi model
#' 
#' @param foi value of force of infection for survey
#' @param pop population in year and place of survey
#' @param age_groups age groups in survey
#' @param vc vaccine coverage, defaults to 0
#' 
#' @return seroprevalence at age groups
#' @export
foi_calculate_seroprevalence = function(foi,
                                        age_groups,
                                        pop,
                                        vc) {
  
  
  #check lengths
  if(length(vc)==1) {
    vc = rep(vc,length(age_groups))
  } else {
    if(length(vc)!=length(age_groups)) {
      stop("foi_calculate_seroprevalence: incompatible length of vaccination coverage vc_\n")
    }
  }
  
  #check vaccine coverage makes sense
  if(!all(!is.na(vc)) | min(vc)<0 | max(vc)>1){
    stop("foi_calculate_seroprevalence: invalid values for vaccination coverage vc_\n")
  }
  
  #calculate average seroprevalence across age groups
  sero1 = c(as.matrix((1-exp(-foi*0:100))*pop)) 
  
  ag = findInterval(0:100,age_groups)
  
  sero_ag = aggregate(sero1, by=list(age_group=ag), sum)
  
  sero_ag$x = sero_ag$x/aggregate(pop,by=list(age_group=ag),sum)$x
  
  # proportion immune due to either vaccination or previous infection:
  sero_out = 1 - (1-sero_ag$x[sero_ag$age_group>0])*(1-vc)
  
  return(sero_out)
}


