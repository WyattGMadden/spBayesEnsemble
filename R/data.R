#' USEPA Air Quality System (AQS) PM2.5 Data For Atlanta Metro Area, 2003-2005 and CMAQ-GRM, MODIS-GRM, and BART, Model Outputs
#'
#' Daily AQS PM2.5 readings matched with CMAQ CMT output, MODIS AOD reading, covariates of interest, and CMAQ-GRM, MODIS-GRM, and BART model outputs in Atlanta Metro Area from 2003-2005
#' Each row of the data is a unique AQS monitor location/day. 
#' Data included only if CMAQ CMT and MODIS AOD outputs are available for that day and monitoring location.
#' (ie. if the monitor is within a 12km x 12km grid square for which CMT data is available).
#'
#' @format ## `cmaq_aqs_matched`
#' A data frame with 11,404 rows and 16 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{pm25}{PM2.5 reading}
#'   \item{aod}{Aerosol optical depth value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{ctm}{Chemical transport model value}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   \item{cmaqgrm_model_estimate}{PM2.5 estimate from CMAQ-GRM model.}
#'   \item{cmaqgrm_model_sd}{PM2.5 standard deviation from CMAQ-GRM model.}
#'   \item{modisgrm_model_estimate}{PM2.5 estimate from MODIS-GRM model.}
#'   \item{modisgrm_model_sd}{PM2.5 standard deviation from MODIS-GRM model.}
#'   \item{bart_model_estimate}{PM2.5 estimate from BART model.}
#'   \item{bart_model_sd}{PM2.5 standard deviation from BART model.}
#'   ...
#' }
#' @source <https://www.nature.com/articles/jes201390>
"la_aqs"


#' Community Multiscale Air Quality (CMAQ) Chemical Transport Model (CMT) Data For Atlanta Metro Area, June 2004
#'
#' Daily CMAQ CMT output, MODIS AOD readings, covariates of interest, and CMAQ-GRM, MODIS-GRM, and BART model prediction outputs, in Atlanta Metro Area for June 2004. 
#' Data included for every 12km x 12km grid square in study area. 
#' Each row of the data is a unique location/day. 
#'
#' @format ## `cmaq_full`
#' A data frame with 72,000 rows and 15 columns:
#' \describe{
#'   \item{time_id}{Week identifier}
#'   \item{space_id}{Location identifier}
#'   \item{spacetime_id}{Year identifier}
#'   \item{aod}{Aerosol optical depth value}
#'   \item{elevation}{Elevation spatial covariate (US Geological Survey)}
#'   \item{forestcover}{Percentage of forest cover spatial covariate (2001 National Land Cover database)}
#'   \item{hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of major roadway length spatial covariate (2001 National Land Cover database)}
#'   \item{lim_hwy_length}{Sum of local roadway lengths spatial covariate (2001 National Land Cover database)}
#'   \item{point_emi_any}{Indicator of PM2.5 primary emission point source spatial covariate (2002 USEPA National Emissions Inventory)}
#'   \item{tmp}{Temperature spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{wind}{Wind spatio-temporal covariate (North American Land Data Assimilation Systems)}
#'   \item{ctm}{Chemical transport model value}
#'   \item{date}{Date of measurement}
#'   \item{x}{Location x-coordinate}
#'   \item{y}{Location y-coordinate}
#'   \item{cmaqgrm_model_estimate}{PM2.5 estimate from CMAQ-GRM model.}
#'   \item{cmaqgrm_model_sd}{PM2.5 standard deviation from CMAQ-GRM model.}
#'   \item{modisgrm_model_estimate}{PM2.5 estimate from MODIS-GRM model.}
#'   \item{modisgrm_model_sd}{PM2.5 standard deviation from MODIS-GRM model.}
#'   \item{bart_model_estimate}{PM2.5 estimate from BART model.}
#'   \item{bart_model_sd}{PM2.5 standard deviation from BART model.}
#'   ...
#' }
#' @source <https://www.nature.com/articles/jes201390>
"la_grid"

