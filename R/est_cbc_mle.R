#Function to calculate the cbc-MLE ---------------------------------------------
cbc_MLE <- function(mle, stage, max_stage, info_vector, c_region, 
                    theta_lower = -0.5, theta_upper = 0.5, tolerance = 1e-6,
                    eps = 1e-5, rel_tol = .Machine$double.eps^0.25,
                    abs_tol = .Machine$double.eps^0.25, subdivisions = 100L){
  
  #Set up function to calculate conditional bias
  cond_bias_fn <- Vectorize(function(theta_tilde){
  #Calculate integrating constant (probability of stopping at stage k)
    int_constant  = psnorm(lower = -Inf, upper = Inf, stage = stage, max_stage = max_stage, c_region = c_region, 
                           info_vector = info_vector, theta = theta_tilde)
  #Calculate scaling value for conditional density
    scaling = 1 / int_constant
    
    
    wrap <- Vectorize(function(x) {
      scaling * dsnorm(x = x, stage = stage,
                       max_stage = max_stage,
                       c_region = c_region,
                       info_vector = info_vector,
                       theta = theta_tilde, eps = eps)
    },vectorize.args = "x")
    
    int_check <- tryCatch(integrate(wrap, lower = -Inf, upper = Inf)$value, error = function(e) { return( 1e6)})
    if(int_check == 1e6) {
      return(1e6)
    }
    if( abs(int_check - 1) > 1e-3) {
      return(1e6)
    }
    

  #Set up function for calculating the integral
    integration_cond = Vectorize(integrate_cond_fn <- function(z){
      z/sqrt(info_vector[stage]) * dsnorm(x = z, stage = stage,
                                      max_stage = max_stage,
                              c_region = c_region,
                              info_vector = info_vector,
                              theta = theta_tilde, eps = eps) * scaling
    }, vectorize.args = "z")
  
  #Integrate over function to get expectation of MLE conditional on stage k
    cond_exp = integrate(integration_cond, lower = -Inf, upper = Inf, rel.tol = rel_tol,  abs.tol = abs_tol, subdivisions = subdivisions)$value
  
  #Calculate conditional bias as expectation minus theta tilde
    cond_bias = cond_exp - theta_tilde
  
  #Return the optimization value
    return((theta_tilde + cond_bias - mle)^2)
    },vectorize.args = "theta_tilde")
  #Optimize the conditional bias function over interval
  cbc_mle_est = optimize(cond_bias_fn, interval = c(theta_lower,theta_upper))$minimum
  
  #Return estimate of cbc-MLE
  return(cbc_mle_est)
}

#cbc_MLE(mle = MLE, info_vector = c(Ihat_1,Ihat_2, 450), c_region = c_region, k = 2)
#cbc_MLE(mle = 0.16,stage = 1, max_stage = 2, info_vector = c(Ihat_1,Ihat_2), c_region = matrix(c(-Inf, e1), ncol = 2, byrow = TRUE), theta_lower = -1, theta_upper= 1)
