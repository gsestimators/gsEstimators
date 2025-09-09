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