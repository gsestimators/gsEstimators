#Input error checking for binary data two-arm-----------------------------------
binary_2_arm_check_sample_sizes <- function(n_c, n_e, s_c, s_e, stage, 
                                            max_stage) {
  
#Set error message list to an empty vector 
  error_message_list <- c()
  
#Ensure integer values inputted for control arm sample sizes
  if(all(floor(n_c) == n_c) == FALSE) { 
    error_message_list <- c(error_message_list, "List of control arm sample 
                            sizes must take integer values")
  }
  
#Ensure integer values inputted for experimental arm sample sizes
  if(all(floor(n_e) == n_e) == FALSE) { 
    error_message_list <- c(error_message_list, "List of experimental arm sample
                            sizes must take integer values")
  }
  
#Ensure integer values inputted for control arm events
  if(all(floor(s_c) == s_c) == FALSE) { 
    error_message_list <- c(error_message_list, "List of control arm events must
                            take integer values")
  }
  
#Ensure integer values inputted for experimental arm events
  if(all(floor(s_e) == s_e) == FALSE) { 
    error_message_list <- c(error_message_list, "List of experimental arm events 
                            must take integer values")
  }
  
#Ensure length of cumulative control sample sizes equals max stages 
  if(length(n_c) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            control arm sample sizes must equal the number
                            of trial stages")
  }
  
#Ensure length of cumulative experimental sample sizes equals max stages 
  if(length(n_e) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            control arm sample sizes must equal the realised 
                            trial stage")
  }
  
#Ensure length of cumulative control total events equals max stages 
  if(length(s_c) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            experimental arm sample sizes must equal the 
                            realised trial stage")
  }
  
#Ensure length of cumulative experimental total events equals max stages 
  if(length(s_e) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative   
                            number of events in the experimental arm must equal 
                            the number of trial stages")
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

#Ensure cumulative control total events is strictly increasing 
  if(vector_increasing(s_c) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative number of     
                            events in control arm must be strictly increasing")
  }
  
#Ensure cumulative experimental total events is strictly increasing 
  if(vector_increasing(s_e) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative number of     
                            events in experimental arm must be strictly 
                            increasing")
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

#Ensure positivity of cumulative control total events 
  if(all(s_c > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            number of events in the control arm must be 
                            positive")
  }

#Ensure positivity of cumulative experimental total events 
  if(all(s_e > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            number of events in the experimental arm must 
                            be positive")
  }
  
#Ensure that total control events never exceeds sample size at each stage 
  if(any(s_c > n_c)) {
    error_message_list <- c(error_message_list, "Number of events in the control 
                            arm must be less than the sample size at each 
                            stage")
  }

#Ensure that total experimental events never exceeds sample size at each stage 
  if(any(s_e > n_e)) {
    error_message_list <- c(error_message_list, "Number of events in the  
                            experimental arm must be less than the sample size 
                            at each stage")
  }
  return(error_message_list)
}