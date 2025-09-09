#Function to calculate the UMVUE -----------------------------------------------
umvue <- function(stage, max_stage, c_region, info_vector, z, theta = 0, eps = 1e-5, 
                    rel_tol = .Machine$double.eps^0.25,
                    abs_tol = .Machine$double.eps^0.25, subdivisions = 100L) {
  
  error_indicator <- 0
  
  if(stage == 1){
    return(z/sqrt(info_vector[1]))
  }
  
  conditional_density <- Vectorize(function(z1) {
    z1_conditional_dsnorm(z1 = z1, z = z, stage = stage, max_stage = max_stage, c_region = c_region, 
                          info_vector = info_vector, theta = theta, eps = eps)
    
  }, vectorize.args = "z1")
  
  
  norm_const <- tryCatch(integrate(conditional_density, lower = c_region[1, 1], 
                          upper = c_region[1, 2], rel.tol = rel_tol,  
                          abs.tol = abs_tol)$value, error = function(e) {
                            error_indicator <- 1
                            return("Error in calculating the normalising constant: check trial data inputs and control inputs")
                          })
  if(error_indicator == 1) {
    return(norm_const)
  }
  scaling <- 1 / (norm_const * sqrt(info_vector[1])) 
  
  umvue <- tryCatch(integrate(function(x) { conditional_density(x) * scaling * x }, 
                     lower = c_region[1, 1],  upper = c_region[1, 2], 
                     rel.tol = rel_tol, abs.tol = abs_tol, subdivisions = 
                       subdivisions)$value,
  error = function(e) {
    error_indicator <- 1
    return("Error in integrating the umvue: check trial data inputs and control inputs")
  })
  
  if(error_indicator == 1) {
    return(umvue)
  }
  
  return(umvue)
}



