#Input error checking for binary data two-arm-----------------------------------
binary_single_arm_check_sample_sizes <- function(succ, n_e, p_0, stage, 
                                                 max_stage) {
  
#Set error message list to an empty vector 
  error_message_list <- c()
  
#Check probability under null is between -1 and 1
  if(p_0< -1 | p_0 >1){
    error_message_list <- c(error_message_list, "Probability under null 
                            hypothesis cannot be greater than 1 or less than 
                            -1")
  }
  
#Check length of probability under null hypothesis is equal to 1
  if(length(p_0) != 1){
    error_message_list <- c(error_message_list, "number of inputted 
                            probabilities under null hypothesis must equal 1")
  }
  
#Ensure integer values inputted for experimental arm sample sizes
  if(all(floor(n_e) == n_e) == FALSE) { 
    error_message_list <- c(error_message_list, "List of experimental arm sample
                            sizes must take integer values")
  }
  
#Ensure length of cumulative experimental sample sizes equals stage reached 
  if(length(n_e) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            control arm sample sizes must equal the realised 
                            trial stage")
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
  
#Ensure integer values inputted for experimental arm events
  if(all(floor(succ) == succ) == FALSE) { 
    error_message_list <- c(error_message_list, "List of experimental arm events 
                            must take integer values")
  } 
  
#Ensure length of cumulative experimental total events equals max stages 
  if(length(succ) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative   
                            number of events in the experimental arm must equal 
                            the number of trial stages")
  } 

#Ensure cumulative experimental total events is strictly increasing 
  if(vector_increasing(succ) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative number of     
                            events in experimental arm must be strictly 
                            increasing")
  } 

#Ensure positivity of cumulative experimental total events 
  if(all(succ > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            number of events in the experimental arm must 
                            be positive")
  }
  
#Ensure that total experimental events never exceeds sample size at each stage 
  if(any(succ > n_e)) {
    error_message_list <- c(error_message_list, "Number of events in the  
                            experimental arm must be less than the sample size 
                            at each stage")
  }

  return(error_message_list)
}