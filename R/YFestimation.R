#' YFestimation: Estimate ecological niche and two transmission models for yellow fever in Africa.
#'
#' The YFestimation package estimates ecological niche model (GLM) and two transmission models (SERO)
#' for yellow fever. The GLM is estimated using conventional MCMC and the two transmission models are
#' estimated with product space MCMC allowing the user to find the model evidence for one
#' relative to the other. Models are based on Garske et al. (2014) and Jean et al. (2019). The
#' estimation approaches are detailed in Gaythorpe et al. (2019).
#'
#' @section YFestimation functions:
#'GLM_MCMC
#'GLM_MCMC_step
#'GLMlike
#'GLMprior
#'GLMproposal
#'Make_aggregate_pop_vc
#'Make_aggregate_pop_vc_3d
#'R0_calculate_seroprevalence
#'R0_recurrence_seroprev
#'R0_recurrence_seroprev_survey
#'R0_recurrence_seroprev_whole
#'SEROlike
#'SEROlike_onesurvey
#'SEROprior
#'SEROproposal
#'SEROpseudoprior
#'Sero_Gibbs_MCMC
#'Sero_Gibbs_MCMC_model_step
#'Sero_Gibbs_MCMC_parameter_step
#'Sero_Gibbs_MCMC_step
#'add_1940_1950
#'calc_incidence_vac_general
#'calc_pop_moments
#'calc_pop_moments_agg
#'calc_t0_vac_africa
#'calculate_bootstrap_Bayes
#'calculate_hpd
#'create_R0_lookup
#'create_pop30_agg_vc30_agg
#'create_pop_at_survey
#'fit_glm
#'foi_calculate_seroprevalence
#'fun_calcPred
#'fun_calc_pdetect_multi_both
#'fun_calc_transmission_Africa
#'get_chains
#'get_pop_data_3d
#'hpd_out
#'import_pop_data_LS2014
#'import_pop_data_LS2017
#'launch_env_dat
#'plot_glm_map
#'plot_prior_post
#'plot_sero_predictions
#'plot_transmission_intensity
#'process_serology
#'setup_modelprior
#'setup_pseudoprior
#'transform_into_vc3d
#' @docType package
#' @name YFestimation
NULL
