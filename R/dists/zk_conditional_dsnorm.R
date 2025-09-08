

##### need to fix
zk_conditional_dsnorm <- Vectorize(function(zk, z = NULL, stage, max_stage, c_region = NULL, info_vector, theta, eps = 1e-5) {
  
  corr_mat <- create_correlation_matrix(info_vector[1:stage])
  mean_vector <- theta * sqrt(info_vector[1:stage])
  if (stage == 1) {
    lower_m1 = c(zk - eps/2)
    upper_m1 = c(zk + eps/2)
    
    m1 <- pmvnorm(lower = lower_m1, upper = upper_m1, mean = mean_vector, sigma = corr_mat[1, 1])[[1]]
    return(m1/eps)
  }
  
  if (stage == 2) {
    lower_m1 = c(zk - eps/2, z - eps/2)
    upper_m1 = c(zk + eps/2, z + eps/2)
    
    lower_m2 = c(c_region[1:(stage - 1), 1], z - eps/2)
    upper_m2 = c(c_region[1:(stage - 1), 2], z + eps/2)
  }
  if(stage > 2) {
    lower_m1 = c(c_region[1:(stage-2), 1],zk - eps/2, z - eps/2)
    upper_m1 = c(c_region[1:(stage-2), 2],zk + eps/2, z + eps/2)
    
    lower_m2 = c(c_region[1:(stage - 1), 1], z - eps/2)
    upper_m2 = c(c_region[1:(stage - 1), 2], z + eps/2)
  }
  
  if(stage == max_stage) {
    lower_m1 = c(c_region[ , 1],zk - eps/2)
    upper_m1 = c(c_region[ , 2],zk + eps/2)
    m1 <- pmvnorm(lower = lower_m1, upper = upper_m1, mean = mean_vector, corr = corr_mat)[[1]]
    return(m1/eps)
  }
  
  m1 <- pmvnorm(lower = lower_m1, upper = upper_m1, mean = mean_vector, corr = corr_mat)[[1]]
  
  m2 <- pmvnorm(lower = lower_m2, upper = upper_m2, mean = mean_vector, corr = corr_mat)[[1]]
  
  return(m1/(m2 * eps))
  
},vectorize.args = "zk")
