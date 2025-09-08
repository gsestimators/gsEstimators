#Normal Estimators input functions ---------------------------------------------

#CODE FOR THE OVERALL MLE FUNCTION ---------------------------------------------
n_ovr_mle_paired <- function(mu_e, mu_c, var, n_p, stage, max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }
  
#Calculate and return the overall MLE 
  mle_ovr = (mu_e[stage] - mu_c[stage])
  return(mle_ovr)
}





#CODE FOR THE STAGE 1 MLE FUNCTION ---------------------------------------------
n_stage1_mle_paired <- function(mu_e, mu_c, var, n_p, stage, max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <- vector_input(var)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }
  
#Calculate and return the stage 1 MLE 
  mle_stg1 = (mu_e[1] - mu_c[1])
  return(mle_stg1) 
}





#CODE FOR THE MUE FUNCTION -----------------------------------------------------
n_mue_paired <- function(mu_e, mu_c, var, n_p, lower, upper, search, stage, 
                         tolerance, max_iter, tail_type, n, max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))
  
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
  data <- data_norm_tran_paired(n = n_p, variance = var, mu_e = mu_e, 
                                mu_c = mu_c, lower = lower, 
                                upper = upper)
  
#Calculate and return the MUE 
  mue_val = tryCatch(mue(stage = stage,
                max_stage = max_stage,
                c_region = data$c_region,
                info_vector = data$info_vector,
                z = data$z[stage],
                tolerance = tolerance,
                theta_lower = search[1],
                theta_upper = search[2],
                tail_type = tail_type,
                max_iter = max_iter,
                n = n),
                error = function(e) {
                  return(paste0("Error when computing estimator please check 
                                input/control parameters ", e$message))
                })
  return(mue_val) 
}





#CODE FOR THE UMVUE FUNCTION ---------------------------------------------------
n_umvue_paired <- function(mu_e, mu_c, var, n_p, lower, upper, 
                               eps, abs_tol, rel_tol, stage, subdivisions,
                               max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))
  
#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }  
  
#Transform data into correct format 
  data <- data_norm_tran_paired(n = n_p, variance = var, mu_e = mu_e, 
                                mu_c = mu_c, lower = lower, 
                                upper = upper)
  
#Calculate and return the MUE 
  umvue = tryCatch(umvue(stage = stage,
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
n_ubc_mle_paired <- function(mu_e, mu_c, var, n_p,  
                                 lower, upper, eps, abs_tol, rel_tol,
                                 search, stage, max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))
  
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
  data <- data_norm_tran_paired(n = n_p, variance = var, mu_e = mu_e, 
                                mu_c = mu_c, lower = lower, 
                                upper = upper)
  
#Calculate overall MLE 
  mle_ovr <- mu_e[stage] - mu_c[stage]
  
#Calculate and return the ubc-MLE   
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
                      return(paste0("Error when computing estimator please check 
                                    input/control parameters ", e$message))
                    })
  return(ubc_mle)
}





#CODE FOR THE cMLE FUNCTION ----------------------------------------------------
n_cmle_paired <- function(mu_e, mu_c, var, n_p, stage, 
                              stage_conditional, max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <- vector_input(var)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
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
  cMLE = (((mu_e[stage_conditional] * n_p[stage_conditional]) - 
            (mu_e[stage_conditional - 1] * n_p[stage_conditional - 1])) / 
            (n_p[stage_conditional] - n_p[stage_conditional - 1])) - 
            (((mu_c[stage_conditional] * n_p[stage_conditional]) - 
            (mu_c[stage_conditional - 1] * n_p[stage_conditional - 1])) / 
            (n_p[stage_conditional] - n_p[stage_conditional - 1])) 
  return(cMLE)
}





#CODE FOR THE cMUE FUNCTION ----------------------------------------------------
n_cmue_paired <- function(mu_e, mu_c, var, n_p, lower, upper, 
                              search, stage, tolerance, max_iter, tail_type, n, 
                              max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))
  
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
  data <- data_norm_tran_paired(n = n_p, variance = var, mu_e = mu_e, 
                                mu_c = mu_c, lower = lower, 
                                upper = upper)
  
#Calculate and return the MUE 
  cMUE_val = tryCatch(cmue(stage = stage,
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
n_umvcue_paired <- function(mu_e, mu_c, var, n_p, lower, upper,
                                eps, abs_tol, rel_tol, stage, subdivisions,
                                max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))

#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  } 
  
#Transform data into correct format 
  data <- data_norm_tran_paired(n = n_p, variance = var, mu_e = mu_e, 
                                mu_c = mu_c, lower = lower, 
                                upper = upper)
  
#Calculate and return the UMVCUE 
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
  return(UMVCUE)
}





#CODE FOR THE cbc-MLE FUNCTION -------------------------------------------------
n_cbc_mle_paired <- function(mu_e, mu_c, var, n_p, lower, 
                                 upper, eps, abs_tol, rel_tol, search, stage,
                                 max_stage){
  
#Vectorize arguments 
  mu_e <- vector_input(mu_e)
  mu_c<- vector_input(mu_c)
  n_p <- vector_input(n_p)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)
  
#Check for any errors in input data 
  errors <- normal_paired_data_check(mu_c = mu_c, mu_e = mu_e, var = var, 
                                     n_p = n_p, stage = stage, 
                                     max_stage = max_stage)
  
#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))
  
#Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))
  
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
  data <- data_norm_tran_paired(n = n_p, variance = var, mu_e = mu_e, 
                                mu_c = mu_c, lower = lower, 
                                upper = upper)
  
#Calculate overall MLE 
  mle_ovr <- mu_e[stage] - mu_c[stage]
  
#Calculate and return the ubc-MLE   
  cbc_mle = tryCatch(cbc_MLE(mle = mle_ovr,
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
                      return(paste0("Error when computing estimator please check 
                                    input/control parameters ", e$message))
                    })
  return(cbc_mle)
}
