#Input error checking for normal (Parallel) data -------------------------------
normal_parallel_data_check <- function(mu_e, mu_c, n_e, n_c, var_e, var_c, 
                                       stage, max_stage){
  
#Set error message list to an empty vector 
  error_message_list <- c()
  
#Check that the variance in the control arm is positive
  if(var_c <= 0){
    error_message_list <- c(error_message_list, 
                            "Variance in the control arm must be strictly 
                            positive")
  }
  
#Check that the variance in the control arm is positive
  if(var_e <= 0){
    error_message_list <- c(error_message_list, 
                            "Variance in the experimental arm must be strictly 
                            positive")
  }
  
#Check that variance in control arm is length 1
  if(length(var_c) != 1){
    error_message_list <- c(error_message_list, 
                            "Length of input of control arm variances must be 
                            length 1")
  }
  
#Check that variance in experimental arm is length 1
  if(length(var_e) != 1){
    error_message_list <- c(error_message_list, 
                            "Length of input of experimental arm variances must 
                            be length 1")
  }

#Ensure integer values inputted for control arm sample sizes
  if(all(floor(n_c) == n_c) == FALSE) { 
    error_message_list <- c(error_message_list, 
                            "List of control arm sample sizes must take 
                            integer values")
  }
  
#Ensure integer values inputted for experimental arm sample sizes
  if(all(floor(n_e) == n_e) == FALSE) { 
    error_message_list <- c(error_message_list, 
                            "List of experimental arm sample sizes must take 
                            integer values")
  }
  
#Ensure length of cumulative control sample sizes equals stage reached 
  if((length(n_c) == stage) == FALSE) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            control arm sample sizes must equal the realised 
                            trial stage")
  }
  
#Ensure length of cumulative experimental sample sizes equals stage reached
  if((length(n_e) == stage) == FALSE) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            experimental arm sample sizes must equal the 
                            realised trial stage")
  }
  
#Ensure cumulative control sample size is strictly increasing 
  if(vector_increasing_strict(n_c) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative control arm 
                            sample sizes must be strictly increasing")
  }
  
#Ensure cumulative experimental sample size is strictly increasing 
  if(vector_increasing_strict(n_e) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative experimental  
                            arm sample sizes must be strictly increasing")
  }
  
#Ensure positivity of cumulative control sample sizes 
  if(all(n_c > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative control arm sample sizes must be 
                            positive")
  }
  
#Ensure positivity of cumulative experimental sample sizes 
  if(all(n_e > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative experimental arm sample sizes must be 
                            positive")
  }
  
#Ensure length of list of control arm means is the same as stage reached -1
  if((length(mu_c) == stage) == FALSE){
    error_message_list <- c(error_message_list, "Length of list of  
                            control arm means must equal the realised 
                            trial stage")
  } 
  
#Ensure length of list of experimental arm means is the same as stage reached -1
  if((length(mu_e) == stage) == FALSE){
    error_message_list <- c(error_message_list, "Length of list of  
                            experimental arm means must equal the realised 
                            trial stage")
  } 
  
  return(error_message_list)
}
