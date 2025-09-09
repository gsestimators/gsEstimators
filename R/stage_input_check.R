#Input error checking for stage inputs -----------------------------------------

stage_input_checks <- function(stage, max_stage){
  
  error_message_list <- c()
  
#Check that stage reached is an integer
  if((floor(stage) == stage) == FALSE){
    error_message_list <- c(error_message_list, "Stage reached must be an 
                            integer value")
  }
  
#Check that max stage is an integer
  if((floor(max_stage) == max_stage) == FALSE){
    error_message_list <- c(error_message_list, "Maximum stages must be an 
                            integer value")
  }
  
#Check that stage is less that or equal to max stage 
  if(stage > max_stage){
    error_message_list <- c(error_message_list, "Stage reached cannot be greater
                            than maximum stages")
  }
#Check positivity of stage reached
  if(stage <= 0){
    error_message_list <- c(error_message_list, "Stage reached cannot be less 
                            than zero")
  }
  
#Check positivity of maximum stages
  if(max_stage <= 0){
    error_message_list <- c(error_message_list, "Maximum stages cannot be less 
                            than zero")
  }
  
  return(error_message_list)
}