
#c_region = matrix(c(-0.6933,3.7103, 1.0042, 2.5114), ncol = 2, byrow = TRUE)
#n_e <- c(188, 376, 564)
# n_c <- c(188, 376, 564)
# s_e <- rbinom(1,c(188, 376, 564 ),0.25)
# s_c <- rbinom(1,c(188, 376, 564 ),0.1)
# 
# data <- data_bnry_tran(n_c,n_e,s_c,s_e, lower = c_region[, 1], upper = c_region[, 2])
# #Function to calculate the UMVCUE ----------------------------------------------
# stage = 1
# max_stage = 3
# info_vector = data$info_vector
# z = data$z[1]
# theta = 0
# eps = 1e-5
# rel_tol = .Machine$double.eps^0.25
# abs_tol = .Machine$double.eps^0.25
# subdivisions = 100L

umvcue <- function(stage, max_stage, c_region, info_vector, z = NULL, theta = 0, eps = 1e-5, 
                     rel_tol = .Machine$double.eps^0.25, 
                     abs_tol = .Machine$double.eps^0.25, subdivisions = 100L) {
  
#Calculate density conditional on zk1 ------------------------------------------
  conditional_density <- Vectorize(function(zk) {
    zk_conditional_dsnorm(zk = zk, z = z, stage = stage, max_stage = max_stage, c_region = c_region, 
                           info_vector = info_vector, theta = theta, eps = eps)
    
  }, vectorize.args = "zk")
  
#Calculate normalizing constant for integrating over expectation ---------------
  norm_const <- integrate(conditional_density, lower = c_region[stage, 1], 
                          upper = c_region[stage, 2], rel.tol = rel_tol,  
                          abs.tol = abs_tol)$value
  
#Calculate scaling constant for integration ------------------------------------
  scaling <- 1 / (norm_const) 
  
#Calculate UMVCUE estimate by integrating over conditional MLE -----------------
  if(stage == 1) {
    
    umvcue_estimate <- integrate(function(x) { conditional_density(x) * 
        scaling * z/sqrt(info_vector[stage]) }, lower = c_region[stage, 1],  
        upper = c_region[stage, 2], rel.tol = rel_tol, abs.tol = abs_tol, 
        subdivisions = subdivisions)$value
    
    #Return UMVCUE -----------------------------------------------------------------
    return(umvcue_estimate)
    
  }  
  
  umvcue_estimate <- integrate(function(x) { conditional_density(x) * 
      scaling * (z *sqrt(info_vector[stage]) - x * sqrt(info_vector[stage-1]))/
      (info_vector[stage] - info_vector[stage-1]) }, lower = c_region[stage, 1],  
      upper = c_region[stage, 2], rel.tol = rel_tol, abs.tol = abs_tol, 
      subdivisions = subdivisions)$value
  
#Return UMVCUE -----------------------------------------------------------------
  return(umvcue_estimate)
}
