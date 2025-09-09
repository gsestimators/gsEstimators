#Input error checking for normal (Paired) data -------------------------------
normal_paired_data_check <- function(mu_e, mu_c, var, n_p, stage, max_stage){
  
#Set error message list to an empty vector 
  error_message_list <- c()
  
#Check that the variance is positive
  if(var <= 0){
    error_message_list <- c(error_message_list, 
                            "Variance must be be strictly positive")
  }
  
#Check that variance inputted is length 1
  if(length(var) != 1){
    error_message_list <- c(error_message_list, 
                            "Length of input of variance must be length 1")
  }
  
#Ensure integer values inputted paired samples at each stage
  if(all(floor(n_p) == n_p) == FALSE) { 
    error_message_list <- c(error_message_list, 
                            "List of paired samples at each stage must take 
                            integer values")
  }
  
#Ensure length of paired samples sizes at each stage equals stage reached
  if((length(n_p) == stage) == FALSE) { 
    error_message_list <- c(error_message_list, "Length of list of paired 
                            samples at each stage must equal the realised trial 
                            stage")
  }
  
#Ensure list of paired samples sizes at each stage is strictly increasing 
  if(vector_increasing_strict(n_p) == FALSE) {
    error_message_list <- c(error_message_list, "List of paired samples at each 
                            stage must be strictly increasing")
  }
  
#Ensure positivity of list of paired samples sizes at each stage 
  if(all(n_p > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            paired samples sizes at each stage must be 
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