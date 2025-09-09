#Function to calculate ubc-MLE -------------------------------------------------
ubc_MLE <- function(mle, stage, max_stage = max_stage, info_vector, c_region, 
                    theta_lower = -0.5, theta_upper = 0.5, tolerance = 1e-6,
                    eps = 1e-5, rel_tol = .Machine$double.eps^0.25,
                    abs_tol = .Machine$double.eps^0.25, subdivisions = 100L){
  
  #Function to calculate the overall bias
  bias_fun <- Vectorize(function(theta_tilde){
    
    #Set up vector of conditional bias
    cond_exp <- rep(0, stage)
    
    #For loop to calculate conditional expectations from stage 1 to k
    for (i in 1:stage) {
      
      #Check if stage 1
      if(i == 1){
        
        #Function to integrate over for stage 1
        integrate_fn <- Vectorize(function(z){
          z/sqrt(info_vector[i]) * dsnorm(x = z, stage = i,
                                          max_stage = max_stage,
                                          c_region = c_region, 
                                          info_vector = info_vector,
                                          theta = theta_tilde, eps = eps)
        },vectorize.args = "z")
        
        #Add lower and upper integral for conditional expectation at stage i
        cond_exp[i] = integrate(integrate_fn, lower = -Inf, 
                                upper = Inf, rel.tol = rel_tol,  
                                abs.tol = abs_tol)$value  
        next
        
      }
      
      #Check if stage k
      if(i == stage){
        #Function to integrate over for stage k
        integration_k <- Vectorize(integrate_fn <- function(z){
          z/sqrt(info_vector[i]) * dsnorm(x = z, stage = i,
                                          max_stage = max_stage,
                                          c_region = c_region,
                                          info_vector = info_vector,
                                          theta = theta_tilde, eps = eps)
        }, vectorize.args = "z")
        
        #Integrate over function to get conditional expectation at stage k
        cond_exp[i] = integrate(integration_k, lower = -Inf, 
                                upper = Inf, rel.tol = rel_tol,  
                                abs.tol = abs_tol)$value
        next
      }
      
      #Check if not stage 1 or k
      else{
        #Function to integrate over for stage i
        integration_mid <- Vectorize(integrate_mid <- function(z){
          z/sqrt(info_vector[i]) * dsnorm(x = z, stage = i, 
                                          max_stage = max_stage,
                                          c_region = c_region,
                                          info_vector = info_vector,
                                          theta = theta_tilde, eps = eps)
        }, vectorize.args = "z")
        
        #Add lower and upper integral for conditional expectation at stage i
        cond_exp[i] = integrate(integration_mid, lower = -Inf, 
                                upper = Inf, rel.tol = rel_tol,  
                                abs.tol = abs_tol)$value  
        next
      }
    }
    
    #Calculate overall bias by summing expectations minus theta tilde
    ovr_bias = sum(cond_exp)
    
    #Return value to optimize over
    return((theta_tilde + ovr_bias - mle - theta_tilde)^2)
  },vectorize.args = "theta_tilde")
  
  #optimize overall bias function 
  ubc_mle_est = optimize(bias_fun, interval = c(theta_lower, theta_upper))$minimum
  
  #Return estimate of ubc-MLE
  return(ubc_mle_est)
}
