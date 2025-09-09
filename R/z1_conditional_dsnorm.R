
#this assumes z is nto in the continuation region
# add error check for z
z1_conditional_dsnorm <- Vectorize(function(z1, z, stage, max_stage, c_region = NULL, info_vector, theta, eps = 1e-5) {
  
  corr_mat <- create_correlation_matrix(info_vector[1:stage])
  mean_vector <- theta * sqrt(info_vector[1:stage])
  
  if (stage == 2) {
    lower_m1 = c(z1 - eps/2, z - eps/2)
    upper_m1 = c(z1 + eps/2, z + eps/2)
    
    lower_m2 = c(c_region[1:(stage - 1), 1], z - eps/2)
    upper_m2 = c(c_region[1:(stage - 1), 2], z + eps/2)
  }
  if(stage > 2) {
    lower_m1 = c(z1 - eps/2,c_region[2:(stage - 1), 1], z - eps/2)
    upper_m1 = c(z1 + eps/2,c_region[2:(stage - 1), 2], z + eps/2)
    
    lower_m2 = c(c_region[1:(stage - 1), 1], z - eps/2)
    upper_m2 = c(c_region[1:(stage - 1), 2], z + eps/2)
  }
  
  m1 <- pmvnorm(lower = lower_m1, upper = upper_m1, mean = mean_vector, corr = corr_mat)[[1]]
  
  m2 <- pmvnorm(lower = lower_m2, upper = upper_m2, mean = mean_vector, corr = corr_mat)[[1]]
  
  return(m1/(m2 * eps))

},vectorize.args = "z1")
