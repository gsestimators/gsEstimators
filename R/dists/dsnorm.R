#' Density Function for Group Sequential Normal Test Statistics
#'
#' Computes the density function p(k,z;θ) for group sequential test statistics as show in ref(the refference to the book).
#' For k=1, this reduces to the standard normal density. For k>1, it computes
#' the conditional density given the continuation region constraints from 
#' previous stages.
#'
#' @param x Numeric value at which to evaluate the density
#' @param k Integer, current end stage number (1, 2, 3, ...)
#' @param c_region Matrix with 2 columns, continuation region bounds for stages 1 to k-1.
#'                 Column 1: lower bounds, Column 2: upper bounds. NULL for k=1.
#' @param info_vector Numeric vector of information levels at each stage
#' @param theta Numeric, true parameter value (effect size)
#' @param eps Numeric, small value for numerical differentiation by finite differences (default: 1e-5)
#'
#' @return Numeric, density value at x
#'
#' @details
#' The function implements the density of Z_k | (Z_1, ..., Z_{k-1}) ∈ C_{k-1}
#' where C_{k-1} is the continuation region up to stage k-1.
#' 
#' For k=1: Returns dnorm(x - θ√I_1)
#' For k>1: Uses numerical differentiation of the multivariate normal CDF
#'
#' @references
#' Jennison, C. and Turnbull, B. W. (2000). 
#' \emph{Group Sequential Methods with Applications to Clinical Trials}. 
#' Chapman & Hall/CRC. Chapter 2.
#'
#' @examples
#' # Stage 1 density
#' dsnorm(x = 1.5, k = 1, info_vector = c(0.5, 1.0), theta = 0.2)
#' 
#' # Stage 2 density with continuation region
#' c_region <- matrix(c(-1, 1), nrow = 1, ncol = 2)
#' dsnorm(x = 1.2, k = 2, c_region = c_region, info_vector = c(0.5, 1.0), theta = 0.2)
#'
dsnorm <- function(x, stage, max_stage, c_region = NULL, info_vector, theta, eps = 1e-5) {
  
  if(stage == max_stage) {
    if (stage == 1) {
      return(dnorm(x - theta * sqrt(info_vector[1])))
    }
    # Stage k > 1: Conditional density via numerical differentiation
    # Create correlation matrix for stages 1 to k
    cor_mat <- create_correlation_matrix(info_vector)
    
    # Mean vector under alternative hypothesis θ
    mean_z <- theta * sqrt(info_vector)
    
    # Integration limits: continuation region bounds + small interval around x
    lower_lim <- c(c_region[, 1], x - (eps/2))
    upper_lim <- c(c_region[, 2], x + (eps/2))
    
    # Compute P(Z_1 ∈ C_1, ..., Z_{k-1} ∈ C_{k-1}, Z_k ∈ [x-ε/2, x+ε/2])
    integral_value <- pmvnorm(lower = lower_lim,
                              upper = upper_lim,
                              mean = mean_z,
                              corr = cor_mat)[[1]]
    
    # Numerical derivative: f(x) ≈ P(x-ε/2 < X < x+ε/2) / ε
    return(integral_value / eps)
  }
  else {
    
    futility <- c_region[stage, 1]
    efficacy <- c_region[stage, 2]
    
    if(x > futility & x < efficacy) {
      return(0)
    }

    if (stage == 1) {
      return(dnorm(x - theta * sqrt(info_vector[1])))
    }
    # Stage k > 1: Conditional density via numerical differentiation
    # Create correlation matrix for stages 1 to k
    cor_mat <- create_correlation_matrix(info_vector[1:stage])
    
    # Mean vector under alternative hypothesis θ
    mean_z <- theta * sqrt(info_vector[1:stage])
    
    # Integration limits: continuation region bounds + small interval around x
    lower_lim <- c(c_region[1:(stage - 1), 1], x - (eps/2))
    upper_lim <- c(c_region[1:(stage - 1), 2], x + (eps/2))
    
    # Compute P(Z_1 ∈ C_1, ..., Z_{k-1} ∈ C_{k-1}, Z_k ∈ [x-ε/2, x+ε/2])
    integral_value <- pmvnorm(lower = lower_lim,
                              upper = upper_lim,
                              mean = mean_z,
                              corr = cor_mat)[[1]]
    
    # Numerical derivative: f(x) ≈ P(x-ε/2 < X < x+ε/2) / ε
    return(integral_value / eps)
    
  }
}



#x = seq(from = -6, to = 6, length = 1000)
#hold <- c()
#for( i in 1:1000) {
 # hold[i] = dsnorm(x[i], stage = 1, max_stage = 2, c_region = matrix(c(-2.797, 2.797),ncol = 2, byrow = TRUE),
  #                 info_vector = c(313,394),theta =0)
#}
#plot(x,hold)
