#Normal Estimators input functions ---------------------------------------------

-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
##                           paired trial design estimators                                  ##
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------

#CODE FOR THE OVERALL MLE FUNCTION ---------------------------------------------
n_ovr_mle_single_arm <- function(mu, mu_null, n_e, var, stage, max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
                                         max_stage = max_stage)

#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }


  #Calculate and return the overall MLE
  mle_ovr = (mu[stage] - mu_null)
  return(mle_ovr)
}





#CODE FOR THE STAGE 1 MLE FUNCTION ---------------------------------------------
n_stage1_mle_single_arm <- function(mu, mu_null, n_e, var, stage, max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
                                         max_stage = max_stage)

#Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

#Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

#Calculate and return the stage 1 MLE
  mle_stg1 = (mu[1] - mu_null)
  return(mle_stg1)
}





#CODE FOR THE MUE FUNCTION -----------------------------------------------------
n_mue_single_arm <- function(mu, mu_null, n_e, var, lower, upper, search, stage,
                           tolerance, max_iter, tail_type, n, max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                       var = var, stage = stage,
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
  data <- data_norm_tran_single_arm(n = n_e, variance = var, mu = mu,
                                    mu_null = mu_null, lower = lower,
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
                                input/control parameters ",e$message))
                })
  return(mue_val)
}





#CODE FOR THE UMVUE FUNCTION ---------------------------------------------------
n_umvue_single_arm <- function(mu, mu_null, n_e, var, lower, upper,
                             eps, abs_tol, rel_tol, stage, subdivisions,
                             max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
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
  data <- data_norm_tran_single_arm(n = n_e, variance = var, mu = mu,
                                    mu_null = mu_null, lower = lower,
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
                                input/control parameters ",e$message))
                })
  return(umvue)
}





#CODE FOR THE ubc-MLE FUNCTION -------------------------------------------------
n_ubc_mle_single_arm <- function(mu, mu_null, n_e, var,
                               lower, upper, eps, abs_tol, rel_tol,
                               search, stage, max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
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
  data <- data_norm_tran_single_arm(n = n_e, variance = var, mu = mu,
                                    mu_null = mu_null, lower = lower,
                                    upper = upper)

#Calculate overall MLE
  mle_ovr <- mu[stage] - mu_null

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
                                    input/control parameters ",e$message))
                    })
  return(ubc_mle)
}





#CODE FOR THE cMLE FUNCTION ----------------------------------------------------
n_cmle_single_arm <- function(mu, mu_null, n_e, var, stage,
                            stage_conditional, max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <- vector_input(var)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
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
  cMLE = (((mu[stage_conditional] * n_e[stage_conditional]) -
            (mu[stage_conditional - 1] * n_e[stage_conditional - 1])) /
            (n_e[stage_conditional] - n_e[stage_conditional - 1])) - mu_null
  return(cMLE)
}





#CODE FOR THE cMUE FUNCTION ----------------------------------------------------
n_cmue_single_arm <- function(mu, mu_null, n_e, var, lower, upper,
                            search, stage, tolerance, max_iter, tail_type, n,
                            max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
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
  data <- data_norm_tran_single_arm(n = n_e, variance = var, mu = mu,
                                    mu_null = mu_null, lower = lower,
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
                                  input/control parameters ",e$message))
                  })
  return(cMUE_val)
}





#CODE FOR THE UMVCUE FUNCTION --------------------------------------------------
n_umvcue_single_arm <- function(mu, mu_null, n_e, var, lower, upper,
                              eps, abs_tol, rel_tol, stage, subdivisions,
                              max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
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
  data <- data_norm_tran_single_arm(n = n_e, variance = var, mu = mu,
                                    mu_null = mu_null, lower = lower,
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
                                  input/control parameters ",e$message))
                  })
  return(UMVCUE)
}





#CODE FOR THE cbc-MLE FUNCTION -------------------------------------------------
n_cbc_mle_single_arm <- function(mu, mu_null, n_e, var, lower,
                               upper, eps, abs_tol, rel_tol, search, stage,
                               max_stage){

#Vectorize arguments
  mu <- vector_input(mu)
  mu_null <- vector_input(mu_null)
  n_e <- vector_input(n_e)
  var <-vector_input(var)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

#Check for any errors in input data
  errors <- normal_single_arm_data_check(mu = mu, mu_null = mu_null, n_e = n_e,
                                         var = var, stage = stage,
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
  data <- data_norm_tran_single_arm(n = n_e, variance = var, mu = mu,
                                    mu_null = mu_null, lower = lower,
                                    upper = upper)

#Calculate overall MLE
  mle_ovr <- mu[stage] - mu_null

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
                                    input/control parameters ",e$message))
                    })
  return(cbc_mle)
}

-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
##                           parallel trial design estimators                                  ##
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------

  #CODE FOR THE OVERALL MLE FUNCTION ---------------------------------------------
n_ovr_mle_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, stage,
                               max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
n_stage1_mle_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, stage,
                                  max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
n_mue_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e,
                           lower, upper, search, stage, tolerance,
                           max_iter, tail_type, n, max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_parallel(n_e = n_e, n_c = n_c, var_e = var_e,
                                  var_c = var_c, mu_e = mu_e, mu_c = mu_c,
                                  lower = lower, upper = upper)

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
n_umvue_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, lower, upper,
                             eps, abs_tol, rel_tol, stage, subdivisions,
                             max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)
  lower <- vector_input(lower)
  upper <- vector_input(upper)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

  #Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))

  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  #Transform data into correct format
  data <- data_norm_tran_parallel(n_e = n_e, n_c = n_c, var_e = var_e,
                                  var_c = var_c, mu_e = mu_e, mu_c = mu_c,
                                  lower = lower, upper = upper)

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
n_ubc_mle_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e,
                               lower, upper, eps, abs_tol, rel_tol,
                               search, stage, max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_parallel(n_e = n_e, n_c = n_c, var_e = var_e,
                                  var_c = var_c, mu_e = mu_e, mu_c = mu_c,
                                  lower = lower, upper = upper)

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
n_cmle_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, stage,
                            stage_conditional, max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
  cMLE = (((mu_e[stage_conditional] * n_e[stage_conditional]) -
             (mu_e[stage_conditional - 1] * n_e[stage_conditional - 1])) /
            (n_e[stage_conditional] - n_e[stage_conditional - 1])) -
    (((mu_c[stage_conditional] * n_c[stage_conditional]) -
        (mu_c[stage_conditional - 1] * n_c[stage_conditional - 1]))/
       (n_c[stage_conditional] - n_c[stage_conditional - 1]))
  return(cMLE)
}





#CODE FOR THE cMUE FUNCTION ----------------------------------------------------
n_cmue_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, lower, upper,
                            search, stage, tolerance, max_iter, tail_type, n,
                            max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_parallel(n_e = n_e, n_c = n_c, var_e = var_e,
                                  var_c = var_c, mu_e = mu_e, mu_c = mu_c,
                                  lower = lower, upper = upper)

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
n_umvcue_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, lower, upper,
                              eps, abs_tol, rel_tol, stage, subdivisions,
                              max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)
  lower <- vector_input(lower)
  upper <- vector_input(upper)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

  #Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))

  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  #Transform data into correct format
  data <- data_norm_tran_parallel(n_e = n_e, n_c = n_c, var_e = var_e,
                                  var_c = var_c, mu_e = mu_e, mu_c = mu_c,
                                  lower = lower, upper = upper)

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
n_cbc_mle_parallel <- function(mu_c, mu_e, n_c, n_e, var_c, var_e, lower,
                               upper, eps, abs_tol, rel_tol, search, stage,
                               max_stage){

  #Vectorize arguments
  mu_e <- vector_input(mu_e)
  mu_c <- vector_input(mu_c)
  n_e <- vector_input(n_e)
  n_c <- vector_input(n_c)
  var_e <-vector_input(var_e)
  var_c <-vector_input(var_c)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_parallel_data_check(mu_e = mu_e, mu_c = mu_c, n_e = n_e,
                                       n_c = n_c, var_e = var_e, var_c = var_c,
                                       stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_parallel(n_e = n_e, n_c = n_c, var_e = var_e,
                                  var_c = var_c, mu_e = mu_e, mu_c = mu_c,
                                  lower = lower, upper = upper)

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
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------
##                           crossover trial design estimators                                 ##
-------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------

  #Normal Estimators input functions ---------------------------------------------

#CODE FOR THE OVERALL MLE FUNCTION ---------------------------------------------
n_ovr_mle_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, stage,
                                max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  #Calculate and return the overall MLE
  mle_ovr = 0.5 * (mu_ab[stage] - mu_ba[stage])
  return(mle_ovr)
}





#CODE FOR THE STAGE 1 MLE FUNCTION ---------------------------------------------
n_stage1_mle_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba,
                                   stage, max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  #Calculate and return the stage 1 MLE
  mle_stg1 = 0.5 * (mu_ab[1] - mu_ba[1])
  return(mle_stg1)
}





#CODE FOR THE MUE FUNCTION -----------------------------------------------------
n_mue_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, lower,
                            upper, search, stage, tolerance, max_iter,
                            tail_type, n, max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_crossover(n_ab = n_ab, n_ba = n_ba, var_ab = var_ab,
                                   var_ba = var_ba, mu_ab = mu_ab,
                                   mu_ba = mu_ba, lower = lower, upper = upper)

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
n_umvue_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, lower,
                              upper, eps, abs_tol, rel_tol, stage, subdivisions,
                              max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)
  lower <- vector_input(lower)
  upper <- vector_input(upper)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

  #Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))

  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  #Transform data into correct format
  data <- data_norm_tran_crossover(n_ab = n_ab, n_ba = n_ba, var_ab = var_ab,
                                   var_ba = var_ba, mu_ab = mu_ab,
                                   mu_ba = mu_ba, lower = lower, upper = upper)

  #Calculate and return the UMVUE
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
n_ubc_mle_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba,
                                lower, upper, eps, abs_tol, rel_tol,
                                search, stage, max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_crossover(n_ab = n_ab, n_ba = n_ba, var_ab = var_ab,
                                   var_ba = var_ba, mu_ab = mu_ab,
                                   mu_ba = mu_ba, lower = lower, upper = upper)

  #Calculate overall MLE
  mle_ovr <- 0.5*(mu_ab[stage] - mu_ba[stage])

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
n_cmle_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, stage,
                             stage_conditional, max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

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
  cMLE = 0.5*((((mu_ab[stage_conditional] * n_ab[stage_conditional]) -
                  (mu_ab[stage_conditional - 1] * n_ab[stage_conditional - 1])) /
                 (n_ab[stage_conditional] - n_ab[stage_conditional - 1])) -
                (((mu_ba[stage_conditional] * n_ba[stage_conditional]) -
                    (mu_ba[stage_conditional - 1] * n_ba[stage_conditional - 1])) /
                   (n_ba[stage_conditional] - n_ba[stage_conditional - 1])))
  return(cMLE)
}





#CODE FOR THE cMUE FUNCTION ----------------------------------------------------
n_cmue_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, lower,
                             upper, search, stage, tolerance, max_iter,
                             tail_type, n,max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_crossover(n_ab = n_ab, n_ba = n_ba, var_ab = var_ab,
                                   var_ba = var_ba, mu_ab = mu_ab,
                                   mu_ba = mu_ba, lower = lower, upper = upper)

  #Calculate and return the cMUE
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
n_umvcue_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, lower,
                               upper, eps, abs_tol, rel_tol, stage,
                               subdivisions, max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)
  lower <- vector_input(lower)
  upper <- vector_input(upper)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

  #Check for any errors in the stage inputs
  errors <- c(errors, stage_input_checks(stage = stage, max_stage = max_stage))

  #Check for any errors in the limits inputs
  errors <- c(errors, check_limits(lower, upper, stage, max_stage))

  #Check and return any errors
  if(is.null(errors) == FALSE) {
    return(errors)
  }

  #Transform data into correct format
  data <- data_norm_tran_crossover(n_ab = n_ab, n_ba = n_ba, var_ab = var_ab,
                                   var_ba = var_ba, mu_ab = mu_ab,
                                   mu_ba = mu_ba, lower = lower, upper = upper)

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
n_cbc_mle_crossover <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, n_ba, lower,
                                upper, eps, abs_tol, rel_tol, search, stage,
                                max_stage){

  #Vectorize arguments
  mu_ab <- vector_input(mu_ab)
  mu_ba <- vector_input(mu_ba)
  n_ab <- vector_input(n_ab)
  n_ba <- vector_input(n_ba)
  var_ab <- vector_input(var_ab)
  var_ba <- vector_input(var_ba)
  lower <- vector_input(lower)
  upper <- vector_input(upper)
  search <- vector_input(search)

  #Check for any errors in input data
  errors <- normal_crossover_data_check(mu_ab = mu_ab, mu_ba = mu_ba,
                                        n_ab = n_ab, n_ba = n_ba,
                                        var_ab = var_ab, var_ba = var_ba,
                                        stage = stage, max_stage = max_stage)

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
  data <- data_norm_tran_crossover(n_ab = n_ab, n_ba = n_ba, var_ab = var_ab,
                                   var_ba = var_ba, mu_ab = mu_ab,
                                   mu_ba = mu_ba, lower = lower, upper = upper)

  #Calculate overall MLE
  mle_ovr <- 0.5 * (mu_ab[stage] - mu_ba[stage])

  #Calculate and return the cbc-MLE
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
