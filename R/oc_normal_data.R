#Transformation function for normal input data ---------------------------------

vector_input <- function(str){
  
  str_vec <- trimws(strsplit(str, ",")[[1]])
  
  as.numeric(str_vec)
}





#Function for a Parallel Two-Treatment comparison------------------------------
data_norm_tran_parallel <- function(n_e, n_c, var_e, var_c, mu_e, mu_c, lower, 
                                    upper){
  
#Calculate information vector --------------------------------------------------
  info_vector = (var_e/n_e + var_c/n_c)^(-1)
  
#Construct matrix of continuation regions --------------------------------------
  c_region = matrix(c(lower, upper), ncol = 2, byrow = FALSE)
  
#Calculate vector of test statistics -------------------------------------------
  z = (mu_e - mu_c)*sqrt(info_vector) 
  
#Return all transformed data ---------------------------------------------------
  return(list(info_vector = info_vector, c_region = c_region, z = z))
}




#Function for a single study population ----------------------------------------
data_norm_tran_single_arm <- function(n, mu, variance, mu_null, lower, upper){
  
#Calculate information vector --------------------------------------------------
  info_vector = n / variance
  
#Construct matrix of continuation regions --------------------------------------
  c_region = matrix(c(lower, upper), ncol = 2, byrow = FALSE)
  
#Calculate vector of test statistics -------------------------------------------
  z = (mu - mu_null)*sqrt(info_vector)
  
#Return all transformed data ---------------------------------------------------
  return(list(info_vector = info_vector, c_region = c_region, z = z))
}




#Function for a paired two treatment design ------------------------------------
data_norm_tran_paired <- function(n, variance, mu_e, mu_c, lower, upper){
  
#Calculate information vector --------------------------------------------------
  info_vector = n / variance
  
#Construct matrix of continuation regions --------------------------------------
  c_region = matrix(c(lower, upper), ncol = 2, byrow = FALSE)
  
#Calculate vector of test statistics -------------------------------------------
  z = (mu_e - mu_c) * sqrt(info_vector)
  
#Return all transformed data ---------------------------------------------------
  return(list(info_vector = info_vector, c_region = c_region, z = z))
}





#Function for a two period crossover -------------------------------------------
data_norm_tran_crossover <- function(n_ab, n_ba, var_ab, var_ba, mu_ab, mu_ba, 
                                     lower, upper){
  
#Calculate information vector --------------------------------------------------
  info_vector = 4 * ((var_ab / n_ab) + (var_ba / n_ba))^(-1)
  
#Construct matrix of continuation regions --------------------------------------
  c_region = matrix(c(lower, upper), ncol = 2, byrow = FALSE)
  
#Calculate vector of test statistics -------------------------------------------
  z = 0.5 * (mu_ab + mu_ba) * sqrt(info_vector)
  
#Return all transformed data ---------------------------------------------------
  return(list(info_vector = info_vector, c_region = c_region, z = z))
}


