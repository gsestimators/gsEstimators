#' Cumulative Distribution Function for Group Sequential Test Statistics
#'
#' Computes the cumulative distribution function for group sequential test 
#' statistics with user-specified integration limits. This allows computation
#' of P(lower < Z_k < upper | continuation to stage k).
#'
#' @param lower Numeric, lower integration limit for Z_k (-Inf for left tail)
#' @param upper Numeric, upper integration limit for Z_k (+Inf for right tail)
#' @param k Integer, current stage number (1, 2, 3, ...)
#' @param c_region Matrix with 2 columns, continuation region bounds for stages 1 to k-1.
#'                 Column 1: lower bounds, Column 2: upper bounds. NULL for k=1.
#' @param info_vector Numeric vector of information levels at each stage
#' @param theta Numeric, true parameter value (effect size)
#'
#' @return Numeric, probability P(lower < Z_k < upper | continuation to stage k)
#'
#' @details
#' The function computes P(lower < Z_k < upper | (Z_1, ..., Z_{k-1}) ∈ C_{k-1})
#' where C_{k-1} is the continuation region up to stage k-1.
#' 
#' For k=1: Returns P(lower < Z_1 < upper) using standard normal CDF
#' For k>1: Uses multivariate normal CDF with continuation region constraints
#' 
#' @references
#' Jennison, C. and Turnbull, B. W. (2000). 
#' \emph{Group Sequential Methods with Applications to Clinical Trials}. 
#' Chapman & Hall/CRC. Chapter 2.
#'
#' @note
#' To get standard CDF behavior:
#' - Use lower = -Inf, upper = x for F(x)
#' - Use lower = x, upper = Inf for 1 - F(x)
#' - Use finite lower and upper for P(lower < X < upper)
#'
#' @examples
#' # Stage 1 CDF
#' psnorm(lower = -Inf, upper = 1.5, k = 1, info_vector = c(0.5, 1.0), theta = 0.2)
#' 
#' # Stage 2 probability in interval, given continuation
#' c_region <- matrix(c(-1, 1), nrow = 1, ncol = 2)
#' psnorm(lower = 0, upper = 2, k = 2, c_region = c_region, 
#'        info_vector = c(0.5, 1.0), theta = 0.2)
#'
psnorm <- function(lower, upper, stage, max_stage, c_region = NULL, info_vector, theta) {
  
  if (lower == upper){
    return(0)
  }
  
  if(stage == max_stage) {
  if (stage == 1) {
    mean_z <- theta * sqrt(info_vector[1])
    
    # P(lower < Z_1 < upper) = P(Z_1 < upper) - P(Z_1 < lower)
    p_lower <- pnorm(lower, mean = mean_z, sd = 1, lower.tail = TRUE)
    p_upper <- pnorm(upper, mean = mean_z, sd = 1, lower.tail = FALSE)
    
    return(1 - (p_lower + p_upper))
  }
  
  # Stage k > 1: Multivariate normal CDF with continuation constraints
  
  # Create correlation matrix for stages 1 to k
  cor_mat <- create_correlation_matrix(info_vector)
  
  # Mean vector under alternative hypothesis θ
  mean_z <- theta * sqrt(info_vector)
  
  # Integration limits: continuation region bounds + [lower, upper] for Z_k
  lower_lim <- c(c_region[, 1], lower)
  upper_lim <- c(c_region[, 2], upper)
  # Compute P(Z_1 ∈ C_1, ..., Z_{k-1} ∈ C_{k-1}, lower < Z_k < upper)
  p1 <- pmvnorm(lower = lower_lim, 
                upper = upper_lim, 
                mean = mean_z, 
                corr = cor_mat)[[1]]
  
  return(p1)
  }
  else {
    
    futility <- c_region[stage, 1]
    efficacy <- c_region[stage, 2]
    
    if(upper >= efficacy & lower >= efficacy) {
      if (stage == 1) {
        mean_z <- theta * sqrt(info_vector[1])
        p_lower <- pnorm(lower, mean = mean_z, sd = 1, lower.tail = TRUE)
        p_upper <- pnorm(upper, mean = mean_z, sd = 1, lower.tail = FALSE)
        return(1 - (p_lower + p_upper))
      }
      cor_mat <- create_correlation_matrix(info_vector[1:stage])
      mean_z <- theta * sqrt(info_vector[1:stage])
      lower_lim <- c(c_region[1:(stage - 1), 1], lower)
      upper_lim <- c(c_region[1:(stage - 1), 2], upper)
      p1 <- pmvnorm(lower = lower_lim, upper = upper_lim,  mean = mean_z, sigma = cor_mat)[[1]]
      return(p1)
    }
    if(lower <= futility & upper <= futility) {
      if (stage == 1) {
        mean_z <- theta * sqrt(info_vector[1])
        p_lower <- pnorm(lower, mean = mean_z, sd = 1, lower.tail = TRUE)
        p_upper <- pnorm(upper, mean = mean_z, sd = 1, lower.tail = FALSE)
        return(1 - (p_lower + p_upper))
      }
      cor_mat <- create_correlation_matrix(info_vector[1:stage])
      mean_z <- theta * sqrt(info_vector[1:stage])
      lower_lim <- c(c_region[1:(stage - 1), 1], lower)
      upper_lim <- c(c_region[1:(stage - 1), 2], upper)
      p1 <- pmvnorm(lower = lower_lim, upper = upper_lim,  mean = mean_z, sigma = cor_mat)[[1]]
      return(p1)
    }
    if(lower <= futility & upper >= futility & upper <= efficacy) {
      if (stage == 1) {
        mean_z <- theta * sqrt(info_vector[1])
        p_lower <- pnorm(lower, mean = mean_z, sd = 1, lower.tail = TRUE)
        p_upper <- pnorm(futility, mean = mean_z, sd = 1, lower.tail = FALSE)
        return(1 - (p_lower + p_upper))
      }
      cor_mat <- create_correlation_matrix(info_vector[1:stage])
      mean_z <- theta * sqrt(info_vector[1:stage])
      lower_lim <- c(c_region[1:(stage - 1), 1], lower)
      upper_lim <- c(c_region[1:(stage - 1), 2], futility)
      p1 <- pmvnorm(lower = lower_lim, upper = upper_lim,  mean = mean_z, sigma = cor_mat)[[1]]
      return(p1)
    }
    if(lower >= futility & lower <= efficacy & upper >= efficacy) {
      if (stage == 1) {
        mean_z <- theta * sqrt(info_vector[1])
        p_lower <- pnorm(efficacy, mean = mean_z, sd = 1, lower.tail = TRUE)
        p_upper <- pnorm(upper, mean = mean_z, sd = 1, lower.tail = FALSE)
        return(1 - (p_lower + p_upper))
      }
      cor_mat <- create_correlation_matrix(info_vector[1:stage])
      mean_z <- theta * sqrt(info_vector[1:stage])
      lower_lim <- c(c_region[1:(stage - 1), 1], efficacy)
      upper_lim <- c(c_region[1:(stage - 1), 2], upper)
      p1 <- pmvnorm(lower = lower_lim, upper = upper_lim,  mean = mean_z, sigma = cor_mat)[[1]]
      return(p1)
    }
    if(lower <= futility & upper >= efficacy) {
      if (stage == 1) {
        mean_z <- theta * sqrt(info_vector[1])
        
        p_lower_1 <- pnorm(lower, mean = mean_z, sd = 1, lower.tail = TRUE)
        p_upper_1 <- pnorm(futility, mean = mean_z, sd = 1, lower.tail = FALSE)
        
        p_lower_2 <- pnorm(efficacy, mean = mean_z, sd = 1, lower.tail = TRUE)
        p_upper_2 <- pnorm(upper, mean = mean_z, sd = 1, lower.tail = FALSE)
        
        part_1 = 1 - (p_lower_1 + p_upper_1)
        part_2 = 1 - (p_lower_2 + p_upper_2)
        
        return(part_1 + part_2)
      }
      cor_mat <- create_correlation_matrix(info_vector[1:stage])
      mean_z <- theta * sqrt(info_vector[1:stage])
      
      lower_lim_1 <- c(c_region[1:(stage - 1), 1], efficacy)
      upper_lim_1 <- c(c_region[1:(stage - 1), 2], upper)
      
      lower_lim_2 <- c(c_region[1:(stage - 1), 1], lower)
      upper_lim_2 <- c(c_region[1:(stage - 1), 2], futility)
      
      part_1 <- pmvnorm(lower = lower_lim_1, upper = upper_lim_1,  mean = mean_z, sigma = cor_mat)[[1]]
      part_2 <- pmvnorm(lower = lower_lim_2, upper = upper_lim_2,  mean = mean_z, sigma = cor_mat)[[1]]
      return(part_1 + part_2)
    }
    else{
      return(0)
    }
  } 
   
}


#psnorm(lower = Inf, upper = Inf, stage = 3, max_stage =3, info_vector = c(Ihat_1, Ihat_2, 450), c_region = matrix(c(-Inf, 2.797, -Inf, 2.101), ncol = 2, byrow = TRUE), theta = 0.1 )
