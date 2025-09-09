#Input error checking for normal (Parallel) data -------------------------------
normal_single_arm_data_check <- function(mu, mu_null, n_e, var, stage, 
                                         max_stage){
  
  #Set error message list to an empty vector 
  error_message_list <- c()
  
  #Check length of probability under null hypothesis is equal to 1
  if(length(mu_null) != 1){
    error_message_list <- c(error_message_list, "number of inputted 
                            means under null hypothesis must equal 1")
  }
  
#Check that the variance is positive
  if(var <= 0){
    error_message_list <- c(error_message_list, 
                            "Variance must be strictly positive")
  } 
  
#Check that inputted variance is of length 1
  if(length(var) != 1){
    error_message_list <- c(error_message_list, 
                            "Length of input of variances must equal 1")
  }
  
#Ensure integer values inputted for control arm sample sizes
  if(all(floor(n_e) == n_e) == FALSE) { 
    error_message_list <- c(error_message_list, 
                            "List of experimental arm sample sizes must take 
                            integer values")
  }
  
#Ensure length of cumulative experimental sample sizes equals stage reached
  if((length(n_e) == stage) == FALSE) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            experimental arm sample sizes must equal the 
                            realised trial stage")
  }
  
#Ensure cumulative experimental sample size is strictly increasing 
  if(vector_increasing_strict(n_e) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative experimental  
                            arm sample sizes must be strictly increasing")
  }
  
#Ensure positivity of cumulative experimental sample sizes 
  if(all(n_e > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative experimental arm sample sizes must be 
                            positive")
  } 
  
#Ensure length of list of experimental arm means is the same as stage reached -1
  if((length(mu) == stage) == FALSE){
    error_message_list <- c(error_message_list, "Length of list of  
                            experimental arm means must equal the realised 
                            trial stage")
  } 
  
  return(error_message_list)
}