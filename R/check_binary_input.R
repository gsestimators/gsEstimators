#Input error checking for binary data ------------------------------------------

#Function to check whether a vector is strictly increasing ---------------------
vector_increasing_strict <- function(input) {
  if(length(input) == 1){
    return(TRUE)
  }
  all(diff(input) > 0)
}



#Function to check whether a vector is increasing ------------------------------
vector_increasing <- function(input) {
  if(length(input) == 1){
    return(TRUE)
  }
  if(all(input == -Inf)){
    return(TRUE)
  }
  all(diff(input) >= 0)
}



#Function to check whether a vector is strictly decreasing ---------------------
vector_decreasing_strict <- function(input) {
  if(length(input) == 1){
    return(TRUE)
  }
  all(diff(input) < 0)
}



#Function to check whether a vector is decreasing ------------------------------
vector_decreasing <- function(input) {
  if(length(input) == 1){
    return(TRUE)
  }
  if(all(input == Inf)){
    return(TRUE)
  }
  all(diff(input) <= 0)
}





#Error checking functions for the Trial Data -----------------------------------
binary_2_arm_check_sample_sizes <- function(n_c, n_e, s_c, s_e, stage, max_stage) {
  
#Set error message list to an empty vector -------------------------------------
  error_message_list <- c()
  
  if(is.integer(n_c) == FALSE) { 
    error_message_list <- c(error_message_list, "List of control arm sample sizes must take integer values")
  }
  
  #Ensure length of cumulative experimental sample sizes equals max stages -------
  if(is.integer(n_e) == FALSE) { 
    error_message_list <- c(error_message_list, "List of experimental arm sample sizes must take integer values")
  }
  
  #Ensure length of cumulative control total events equals max stages ------------
  if(is.integer(s_c) == FALSE) { 
    error_message_list <- c(error_message_list, "List of control arm events must take integer values")
  }
  if(is.integer(s_e) == FALSE) { 
    
    #Ensure length of cumulative experimental total events equals max stages -------
    error_message_list <- c(error_message_list, "List of experimental arm events must take integer values")
  }
#Ensure length of cumulative control sample sizes equals max stages ------------
  if(length(n_c) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            control arm sample sizes must equal the number
                            of trial stages")
  }
  
#Ensure length of cumulative experimental sample sizes equals max stages -------
  if(length(n_e) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative  
                            experimental arm sample sizes must equal the 
                            number of trial stages")
  }
  
#Ensure length of cumulative control total events equals max stages ------------
  if(length(s_c) != stage) { 
    error_message_list <- c(error_message_list, "Length of list of cumulative   
                            number of events in the control arm must equal the
                            number of trial stages")
  }
  if(length(s_e) != stage) { 

#Ensure length of cumulative experimental total events equals max stages -------
    error_message_list <- c(error_message_list, "Length of list of cumulative   
                            number of events in the experimental arm must equal 
                            the number of trial stages")
  }
  
#Ensure cumulative control sample size is strictly increasing ------------------
  if(vector_increasing_strict(n_c) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative control arm 
                            sample sizes must be strictly increasing")
  }

#Ensure cumulative experimental sample size is strictly increasing -------------
  if(vector_increasing_strict(n_e) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative experimental  
                            arm sample sizes must be strictly increasing")
  }

#Ensure cumulative control total events is strictly increasing -----------------
  if(vector_increasing(s_c) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative number of     
                            events in control arm must be strictly increasing")
  }
  
#Ensure cumulative experimental total events is strictly increasing ------------
  if(vector_increasing(s_e) == FALSE) {
    error_message_list <- c(error_message_list, "List of cumulative number of     
                            events in experimental arm must be strictly 
                            increasing")
  }
  
#Ensure positivity of cumulative control sample sizes --------------------------
  if(all(n_c > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative control arm sample sizes must be 
                            positive")
  }
  
#Ensure positivity of cumulative experimental sample sizes ---------------------
  if(all(n_e > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            cumulative experimental arm sample sizes must be 
                            positive")
  }

#Ensure positivity of cumulative control total events --------------------------
  if(all(s_c > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            number of events in the control arm must be 
                            positive")
  }

#Ensure positivity of cumulative experimental total events ---------------------
  if(all(s_e > 0) == FALSE) {
    error_message_list <- c(error_message_list, "All values in list of 
                            number of events in the experimental arm must 
                            be positive")
  }
  
#Ensure that total control events never exceeds sample size at each stage ------
  if(any(s_c > n_c)) {
    error_message_list <- c(error_message_list, "Number of events in the control 
                            arm must be less than the sample size at each 
                            stage")
  }

#Ensure that total experimental events never exceeds sample size at each stage -
  if(any(s_e > n_e)) {
    error_message_list <- c(error_message_list, "Number of events in the  
                            experimental arm must be less than the sample size 
                            at each stage")
  }
  return(error_message_list)
}





#Error checking functions for the futility and superiority bounds --------------
check_limits <- function(lower, upper, stage, max_stage){
  
#Set error list to empty vector ------------------------------------------------
  error_message_list <- c()
  
#Check to see if either superiority or futility bound inputs are empty ---------
  if(is.null(lower) | is.null(upper)){
    error_message_list <- c(error_message_list, "Lists of Superiority and  
                            Futility bounds cannot be empty")
    return(error_message_list)
  }
  
#Check that superiority and futility bound inputs are the same length ----------
  if(length(lower) != length(upper)) {
    error_message_list <- c(error_message_list, "List of Superiority bounds must  
                            be same length as the list Futility bounds")
    return(error_message_list)
  }
  
#Ensure that futility bounds are lower than superiorty at each stage -----------
  if(all(lower < upper) == FALSE) {
    error_message_list <- c(error_message_list, "All values in the list of  
                            Superiority bounds must be greater than Futility 
                            bounds")
  }
  
#Ensure that futility bounds are increasing ------------------------------------
  if(vector_increasing(lower) == FALSE) {
    error_message_list <- c(error_message_list, "List of Futility bounds must be 
                            non-decreasing")
  }
  
#Ensure that superiority bounds are decreasing ---------------------------------
  if(vector_decreasing(upper) == FALSE) {
    error_message_list <- c(error_message_list, "List of Superiority bounds must  
                            be decreasing")
  }
  
#Ensure length of futility bound input is one less than the max stage ----------
  if(length(lower) != max_stage - 1) {
    error_message_list <- c(error_message_list, "Length of list of Futility  
                            bounds should be one less than the maximum number of stages")
  }
  
#Ensure length of superiority bound input is one less than max stage -----------
  if(length(upper) != stage - 1) {
    error_message_list <- c(error_message_list, "Length of list of Superiority  
                            bounds should be one less than the maximum number of 
                            stages")
  }
  
#Return list of errors ---------------------------------------------------------
  return(error_message_list)
}

