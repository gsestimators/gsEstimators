#Binary Estimators input functions ---------------------------------------------

#CODE FOR THE OVERALL MLE FUNCTION ---------------------------------------------
b_ovr_mle_single <- function(succ, n_e, p_0, stage, max_stage){
  
#Vectorize arguments
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  
#Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }
  
#Calculate and return the overall MLE
  mle_ovr = (succ[stage]/n_e[stage]) - p_0
  return(mle_ovr)
}





#CODE FOR THE STAGE 1 MLE FUNCTION ---------------------------------------------
b_stg1_mle_single <- function(succ, n_e, p_0, stage, max_stage){
  
#Vectorize arguments 
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  
#Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  
#Calculate and return the stage 1 MLE 
  mle_stg1 = (succ[1]/n_e[1]) - p_0
  return(mle_stg1)
}






#CODE FOR THE MUE FUNCTION -----------------------------------------------------
b_mue_single <- function(succ, n_e, p_0, lower, upper,
                  search, stage, tolerance, max_iter, tail_type, n, max_stage){
  
#Vectorize arguments
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits
  errors <- c(errors, check_limits(lower, upper, stage = stage, 
                                   max_stage = max_stage))
  
#Check length of search region 
  if(length(search) != 2){
    errors <- c(errors, "Search region must be of length 2")
  }
  
#Check if search region lower bound is less than upper
  if(length(search) == 2){
    if(search[1]>=search[2]){
      errors <- c(errors, "Upper search limit must be greater than lower")
    }
  }
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  } 

  
#Transform data into correct format 
  data <- data_bnry_tran_single_arm(n = n_e, s = succ, p_0 = p_0, 
                                    lower = lower, upper = upper)
  
#Calculate and return the MUE 
  mue_val <- tryCatch(mue(stage = stage, 
                 max_stage = max_stage,
                 c_region = data$c_region,
                 info_vector = data$info_vector,
                 z = data$z[stage],
                 tolerance = tolerance,
                 theta_lower = search[1],
                 theta_upper = search[2],
                 tail_type = tail_type,
                 max_iter = max_iter,
                 n=n),
                 error = function(e) {
                   return(paste0("Error when computing estimator please check 
                                 input/control parameters ", e$message))
                 })
  return(mue_val)
}






#CODE FOR THE UMVUE FUNCTION ---------------------------------------------------
b_UMVUE_single <- function(succ, n_e, p_0, lower, upper, eps, abs_tol, rel_tol, stage, 
                    subdivisions, max_stage){
  
#Vectorize arguments 
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  
#Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits
  errors <- c(errors, check_limits(lower, upper, stage = stage, 
                                   max_stage = max_stage))
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

#Transform data into correct format 
  data <- data_bnry_tran_single_arm(n = n_e, s = succ, p_0 = p_0, 
                                    lower = lower, upper = upper)
  
#Calculate and return the UMVUE 
  umvue <- tryCatch(umvue(stage = stage, 
                 max_stage = max_stage, 
                 c_region = data$c_region,
                 info_vector = data$info_vector,
                 z = data$z[stage], 
                 eps = eps,
                 rel_tol = rel_tol,
                 abs_tol = abs_tol,
                 subdivisions = subdivisions),
                 error = function(e) {
                   return(paste0("Error when computing estimator please check
                                 input/control parameters ", e$message))
                 })
  
  return(umvue)
}






#CODE FOR THE ubc-MLE FUNCTION -------------------------------------------------
b_ubc_MLE_single <- function(succ, n_e, p_0, lower, upper, eps, abs_tol, rel_tol,
                      search, stage, max_stage){
  
#Vectorize arguments 
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits
  errors <- c(errors, check_limits(lower, upper, stage = stage, 
                                   max_stage = max_stage))
  
#Check length of search region 
  if(length(search) != 2){
    errors <- c(errors, "Search region must be of length 2")
  }
  
#Check if search region lower bound is less than upper
  if(length(search) == 2){
    if(search[1]>=search[2]){
      errors <- c(errors, "Upper search limit must be greater than lower")
    }
  }
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

#Transform data into correct format 
  data <- data_bnry_tran_single_arm(n = n_e, s = succ, p_0 = p_0, 
                                    lower = lower, upper = upper)
  
#Calculate overall MLE
  mle_ovr = (succ[stage]/n_e[stage]) - p_0
  
#calculate and return the ubc-MLE
  ubc_mle = tryCatch(ubc_MLE(mle = mle_ovr,
                    stage = stage,
                    max_stage = max_stage,
                    info_vector = data$info_vector,
                    c_region = data$c_region,
                    theta_lower = search[1],
                    theta_upper = search[2],
                    eps = eps,
                    rel_tol = rel_tol,
                    abs_tol = abs_tol),
                    error = function(e) {
                      return(paste0("Error when computing estimator please 
                                    check input/control parameters ", 
                                    e$message))
                    })
  
  return(ubc_mle)
}






#CODE FOR THE cMLE FUNCTION ----------------------------------------------------
b_cMLE_single <- function(succ, n_e, p_0, stage, stage_conditional, max_stage){
  
#Vectorize arguments 
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  
  #Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                p_0 = p_0, stage = stage, 
                                                max_stage = max_stage)
  
  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
  #Check whether trying to condition on stage 1 
  if(stage_conditional == 1){
    errors <- c(errors, "Conditioning on first stage is just the stage 1 MLE")
  }
  
  #Check whether conditioning on a stage larger than the maximum stage reached 
  if(stage_conditional > stage){
    errors <- c(errors, "Conditioning on a stage greater than maximum trial 
                stage")
  }
  
  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }
  
#Calculate and return the cMLE
  cMLE = (succ[stage_conditional] - succ[stage_conditional-1])/
    (n_e[stage_conditional] - n_e[stage_conditional-1]) - p_0
  return(cMLE)
}






#CODE FOR THE cMUE FUNCTION ----------------------------------------------------
b_cMUE_single <- function(succ, n_e, p_0, lower, upper,
                   search, stage, tolerance, max_iter, tail_type, n, max_stage){
  
#Vectorize arguments
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits
  errors <- c(errors, check_limits(lower, upper, stage = stage, 
                                   max_stage = max_stage))
  
#Check length of search region 
  if(length(search) != 2){
    errors <- c(errors, "Search region must be of length 2")
  }
  
#Check if search region lower bound is less than upper
  if(length(search) == 2){
    if(search[1]>=search[2]){
      errors <- c(errors, "Upper search limit must be greater than lower")
    }
  }
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

#Transform data into correct format
  data <- data_bnry_tran_single_arm(n = n_e, s = succ, p_0 = p_0, 
                                    lower = lower, upper = upper)
  
#Calculate and return the cMUE
  cMUE_val <- tryCatch(cmue(stage = stage, 
                   max_stage = max_stage,
                   c_region = data$c_region,
                   info_vector = data$info_vector,
                   z = data$z[stage],
                   tolerance = tolerance,
                   theta_lower = search[1],
                   theta_upper = search[2]),
                   error = function(e) {
                     return(paste0("Error when computing estimator please check 
                                   input/control parameters ", e$message))
                   })
  return(cMUE_val)
}






#CODE FOR THE UMVCUE FUNCTION --------------------------------------------------
b_umvcue_single <- function(succ, n_e, p_0, lower, upper, eps, abs_tol, rel_tol, 
                            stage, subdivisions, max_stage){
  
#Vectorize arguments 
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  
  #Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
  #Check for any errors in the limits
  errors <- c(errors, check_limits(lower, upper, stage = stage, 
                                   max_stage = max_stage))
  
  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

#Transform data into correct format
  data <- data_bnry_tran_single_arm(n = n_e, s = succ, p_0 = p_0, 
                                    lower = lower, upper = upper)
  
#Calculate UMVCUE 
  UMVCUE = tryCatch(umvcue(stage = stage, 
                  max_stage = max_stage, 
                  c_region = data$c_region,
                  info_vector = data$info_vector,
                  z = data$z[stage], 
                  eps = eps,
                  rel_tol = rel_tol,
                  abs_tol = abs_tol,
                  subdivisions = subdivisions),
                  error = function(e) {
                    return(paste0("Error when computing estimator please check 
                                  input/control parameters ", e$message))
                  })
  
#Return UMVCUE
  return(UMVCUE)
}






#CODE FOR THE cbc-MLE FUNCTION -------------------------------------------------
b_cbc_MLE_single <- function(succ, n_e, p_0, lower, upper, eps, abs_tol, rel_tol,
                      search, stage, stage_conditional, max_stage){
  
#Vectorize arguments 
  succ <- vector_input(succ)
  n_e <- vector_input(n_e)
  p_0 <- vector_input(p_0)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
  #Check for any errors in input data 
  errors <- binary_single_arm_check_sample_sizes(succ = succ, n_e = n_e, 
                                                 p_0 = p_0, stage = stage, 
                                                 max_stage = max_stage)
  
  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
  #Check for any errors in the limits
  errors <- c(errors, check_limits(lower, upper, stage = stage, 
                                   max_stage = max_stage))
  
  #Check length of search region 
  if(length(search) != 2){
    errors <- c(errors, "Search region must be of length 2")
  }
  
  #Check if search region lower bound is less than upper
  if(length(search) == 2){
    if(search[1]>=search[2]){
      errors <- c(errors, "Upper search limit must be greater than lower")
    }
  }
  
  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

#Transform data into correct format 
  data <- data_bnry_tran_single_arm(n = n_e, s = succ, p_0 = p_0, 
                                    lower = lower, upper = upper)
  
#Calculate overall MLE
  mle_ovr = (succ[stage]/n_e[stage]) - p_0
  
#Calculate and return the cbc-MLE
  cbc_mle_est = tryCatch(cbc_MLE(mle = mle_ovr,
                        stage = stage,
                        max_stage = max_stage,
                        info_vector = data$info_vector,
                        c_region = data$c_region,
                        theta_lower = search[1],
                        theta_upper = search[2],
                        eps = eps,
                        rel_tol = rel_tol,
                        abs_tol = abs_tol),
                        error = function(e) {
                          return(paste0("Error when computing estimator please 
                                        check input/control parameters ", 
                                        e$message))
                        })
  
  return(cbc_mle_est)
}
