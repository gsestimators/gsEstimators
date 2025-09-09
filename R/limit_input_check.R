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