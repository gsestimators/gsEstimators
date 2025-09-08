#Transformation function for binary input data ---------------------------------

#Function to split string inputs -----------------------------------------------
vector_input <- function(str){
  
  str_vec <- trimws(strsplit(str, ",")[[1]])
  
  as.numeric(str_vec)
}





#Function to calculate operating characteristics for Two-arm binary data -------
data_bnry_tran <- function(n_c, n_e, s_c, s_e, lower, upper){
  
#Calculate pooled probabilities 
  p_tilde = (s_e + s_c)/(n_e + n_c)
  
#Calculate information vector 
  info_vector = 1/(p_tilde*(1-p_tilde)*(1/n_e + 1/n_c))
  
#Construct matrix of continuation regions 
  c_region = matrix(c(lower, upper), ncol = 2, byrow = FALSE)
  
#Calculate vector of test statistics 
  z = (s_e/n_e - s_c/n_c)*sqrt(info_vector)
  
#Return all transformed data 
  return(list(info_vector = info_vector, c_region = c_region, z = z))
}





#Function to calculate operating characteristics for single arm binary data ----
data_bnry_tran_single_arm <- function(n, s, p_0, lower, upper){
  
#Calculate information vector 
  info_vector = n / (p_0 * (1 - p_0))

#Construct matrix of continuation regions   
  c_region = matrix(c(lower, upper), ncol = 2, byrow = FALSE)
  
#Calculate vector of test statistics 
  z = ((s / n) - p_0)*sqrt(info_vector)
  
#Return all transformed data 
  return(list(info_vector = info_vector, c_region = c_region, z = z))
}

data_test_2 = data_bnry_tran_single_arm(n = c(101, 143), s = c(27,42), p_0 = 0.1, lower = -Inf, upper = 2.797)
