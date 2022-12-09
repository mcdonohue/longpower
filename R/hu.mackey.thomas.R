#' @Title Function to Produce RCRM MSDR & Needed Statistics
#'
#' @param eff_data Efficacy Dataset, pre-filtered to only population and/or parameter of interest. 
#'          Dataset should include needed covariates for the model. Include Baseline & Post Baseline Assessments. 
#' @param outcome_var A string of the name of the outcome variable of interest (e.g. "AVAL")
#' @param sbj_var A string (if applicable) of the name of the Stuyd ID variable of interest (e.g. "SUBJID")
#' @param covar_vars A vector (if applicable) of strings of the names of the variable(s) of interest
#' @param time_var A string of the name of the variable that contains the relative time component of the efficacy variable. 
#'         (For example relative day compared to randomization date or treatment start date. i.e. "ADY")
#' @param time_unit A string of the time unit of interest. Possible options are "WEEKS", "DAYS", "YEARS", "MONTHS". 
#' @param assess_time A numeric vector of the scheduled assessments defined by the time-unit of interest (e.g c(0, .5, 1, 1.5, 2))
#'
#' @return Outputs RCRM Input Model & RCRM MSDR statistics
#' @export
#'
#' @examples
#' \dontrun{
#' library(rice)
#' library(dplyr)
#' library(nlme)
#' source("~/RCRM_MSDR/RCRM_MSDR_fun.R")
#'
#' #Part 1 Use CREAD 1/2 Pooled Data & Compute RCRM MSDR -----
#' #Read in CREAD Data
#' vads_path <- "root/clinical_studies/RO5490245/CDPT7412/share/pool_PostCSR_Main/prod/outdata_vad/"
#' acdr <- rice_read(paste0(vads_path, "acdr.sas7bdat")) # Efficacy Dataset with CDR Sum of Boxes Endpoint
#' rice_session_close()
#'
#' #Process Data if Required
#' acdr_fltrd <- acdr %>% filter(MITTFL == "Y" & ANL04FL == "Y"  &
#'                                 PARAM == 'CDR: Sum of Boxes' &
#'                                 TRT01P == "Placebo")
#'
#' #Call RCRM_MSDR Function
#' result_CREADPool <- RCRM_MSDR_fun( eff_data = acdr_fltrd
#'                                   ,outcome_var = "AVAL"
#'                                   ,sbj_var = "USUBJID"
#'                                   ,covar_vars = c("STUDYID", "STRATDS")
#'                                   ,time_var = "ADY"
#'                                   ,time_unit = "YEARS"
#'                                   ,assess_time = c(0.5, 1, 1.5, 2)
#'                                   )
#'
#' #CREAD Pooled MSDR Results
#' print(result_CREADPool$outputs_MSDR$MSDR.df[, c(1,4)])
#'
#' #CREAD Pooled RCRM Model Results
#' print(result_CREADPool$outputs_mod$RCRM.res.df)
#'
#' }
#'

RCRM_MSDR_fun <- function(
   eff_data
  ,outcome_var = "AVAL"
  ,sbj_var = "USUBJID"
  ,covar_vars = NULL
  ,time_var = "ADY"
  ,time_unit = "YEARS"
  ,assess_time
) {
  
  # Create RCRM Unit (Year/Weeks/Months/Days) Type of variable 
  # For Baseline Records Set Analysis Relative Day to 1
  
  if (max(names(eff_data) == "ABLFL") == 1) {  #Baseline Flag is present in the dataset. 
    eff_data[[time_var]] = ifelse(eff_data[["ABLFL"]] == "Y", 1, eff_data[[time_var]])
  }
  
  if (time_unit == "YEARS") {
    eff_data[["rcrm.var"]] = eff_data[[time_var]]/365.25
  } else if (time_unit == "MONTHS") {
    eff_data[["rcrm.var"]] = eff_data[[time_var]]/30.43
  } else if (time_unit == "WEEKS") {
    eff_data[["rcrm.var"]] = eff_data[[time_var]]/7
  } else if (time_unit == "DAYS") {
    eff_data[["rcrm.var"]] = eff_data[[time_var]]
  }
  
  # Covariates Processing
  covar_vars = paste0(covar_vars, collapse = " + ")
  
  # Independent Part of Model Processing
  RCRM_dep_formula = gsub( '^[+]', "", paste0(covar_vars, " + ", "rcrm.var"))
  
  ##Formula for Fixed Effect RCRM Model
  RCRM_fix_mod = as.formula(paste0(outcome_var, "~", RCRM_dep_formula))
  
  ##Formula for Random Effect RCRM Model
  RCRM_rndm_mod <- as.formula(paste0("~", "(", "rcrm.var", ")", "|", sbj_var))
  
  ##Run the Model Part 1 Pooled Placebo Data CREAD 1/2 ----
  RCRM_res <- nlme::lme(  RCRM_fix_mod
                         ,random = RCRM_rndm_mod
                         ,data = eff_data
                         ,na.action = na.omit
                         ,method = "REML"
                         ,control = lmeControl(opt="optim"))
  
  ##Get Needed Statistics from the RCRM Result ----
  beta <- as.numeric(fixed.effects(RCRM_res)["rcrm.var"])
  sigmasq.alpha <- as.numeric(VarCorr(RCRM_res)['(Intercept)', 'Variance'])
  sigmasq.beta <- as.numeric(VarCorr(RCRM_res)["rcrm.var", 'Variance'])
  sigmasq <- as.numeric(VarCorr(RCRM_res)['Residual', 'Variance'])
  row <- as.numeric(VarCorr(RCRM_res)["rcrm.var",'Corr'])
  
  ##Computation of MSDR
  t <- assess_time
  t.length <- length(assess_time)
  MSDR.num <- rep(NA, t.length)
  MSDR.deno <- rep(NA, t.length)
  MSDR <- rep(NA, t.length)
  MSDR.time <- rep(NA, t.length)
  
  for (i in 1:t.length) {
    
    MSDR.num[i] <- as.numeric(beta * t[i])
    MSDR.deno[i] <- as.numeric(sqrt(2*sigmasq + sigmasq.beta*(t[i])**2 - ((t[i]*row*sqrt(sigmasq.alpha)*sqrt(sigmasq.beta) - sigmasq)**2/(sigmasq.alpha + sigmasq))))
    
    MSDR[i] <- as.numeric(MSDR.num[i]/MSDR.deno[i])
    
    MSDR.time[i] <- paste0("Time: ", as.character(t[i]))
    
  }
  
  RCRM.res.df <- data.frame(beta,
                            sigmasq,
                            sigmasq.alpha,
                            sigmasq.beta,
                            row)
  
  MSDR.df <- data.frame(MSDR.time,
                        MSDR.num,
                        MSDR.deno,
                        MSDR,
                        stringsAsFactors = FALSE)
  
  result_lst <- list(
    inputs = list(RCRM_fixed_model = RCRM_fix_mod
                  ,RCRM_rndm_modl = RCRM_rndm_mod
                  ,time = assess_time),
    
    outputs_mod = list( RCRM_res = RCRM_res
                        ,RCRM.res.df = RCRM.res.df),
    
    outputs_MSDR = list(MSDR.df = MSDR.df)
  )
}
