#Input error checking for normal (Crossover) data ------------------------------
normal_crossover_data_check <- function(mu_ab, mu_ba, var_ab, var_ba, n_ab, 
                                        n_ba, stage, max_stage){
  
  #Set error message list to an empty vector 
  error_message_list <- c()
  
#Check that the variances in the AB period arm inputted is positive
  if(var_ab <= 0){
    error_message_list <- c(error_message_list, 
                            "Variance in the AB period arm must be strictly 
                            positive")
  }
  
#Check that the length of list of variances in the AB period arm is length 1
  if(length(var_ab) != 1){
    error_message_list <- c(error_message_list, 
                            "Length of list of inputted variances in the AB 
                            period arm must be length 1")
  }
  
#Check that the variances in the BA period arm inputted is positive
  if(var_ba <= 0){
    error_message_list <- c(error_message_list, 
                            "Variance in the BA period arm must be strictly 
                            positive")
  }
  
#Check that the length of list of variances in the BA period arm is length 1
  if(length(var_ba) != 1){
    error_message_list <- c(error_message_list, 
                            "Length of list of inputted variances in the BA 
                            period arm must be length 1")
  }
  
#Ensure integer values inputted for AB arm sample sizes
  if(all(floor(n_ab) == n_ab) == FALSE) { 
    error_message_list <- c(error_message_list, 
                            "List of AB period arm sample sizes must take 
                            integer values")
  }
  
#Ensure integer values inputted for BA arm sample sizes
  if(all(floor(n_ba) == n_ba) == FALSE) { 
    error_message_list <- c(error_message_list, 
                            "List of BA period arm sample sizes must take 
                            integer values")
  }
  
#Ensure length of AB period arm sample sizes equals stage reached 
  if((length(n_ab) == stage) == FALSE) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            AB period arm sample sizes must equal the realised 
                            trial stage")
  }
  
#Ensure length of BA period arm sample sizes equals stage reached
  if((length(n_ba) == stage) == FALSE) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            BA period arm sample sizes must equal the 
                            realised trial stage")
  }
  
#Ensure cumulative AB period arm sample size is strictly increasing 
  if(vector_increasing_strict(n_ab) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative AB period 
                            arm sample sizes must be strictly increasing")
  }
  
#Ensure cumulative BA period arm sample size is strictly increasing 
  if(vector_increasing_strict(n_ba) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative BA period 
                            arm sample sizes must be strictly increasing")
  }
  
#Ensure positivity of cumulative AB period arm sample sizes 
  if(all(n_ab > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative AB period arm sample size must be 
                            positive")
  }
  
#Ensure positivity of cumulative BA period arm sample sizes 
  if(all(n_ba > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative BA period arm sample size must be 
                            positive")
  }
  
#Ensure length of list of AB period arm means is the same as stage reached -1
  if((length(mu_ab) == stage) == FALSE){
    error_message_list <- c(error_message_list, "Length of list of  
                            AB period arm means must equal the realised 
                            trial stage")
  } 
  
#Ensure length of list of BA period arm means is the same as stage reached -1
  if((length(mu_ba) == stage) == FALSE){
    error_message_list <- c(error_message_list, "Length of list of  
                            BA period arm means must equal the realised 
                            trial stage")
  }
  
  return(error_message_list)
}