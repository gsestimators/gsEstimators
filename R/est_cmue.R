#' Conditional Median Unbiased Estimator
#'
#' Computes the conditional median unbiased estimator by finding the parameter value θ such 
#' that the rejection probability equals 0.5 under the alternative 
#' hypothesis.
#'
#' @param k Integer, total number of stages in the design
#' @param c_region Matrix with 2 columns and k - 1 rows, containing continuation region bounds.
#'                 Column 1: lower bounds, Column 2: upper bounds for each stage
#' @param info_vector Numeric vector of information levels at each stage (length k)
#' @param z Numeric, observed test statistic value at the final stage
#' 
#' Control parameters for the root finding:
#' @param tolerance Numeric, tolerance for root finding (default: 1e-6)
#' @param theta_lower Numeric, lower bound for θ search (default: 0)
#' @param theta_upper Numeric, upper bound for θ search (default: 5)
#' See ?uniroot for more details
#'
#' @return Numeric value representing the conditional median unbiased estimate of theta
#' 
# stage = 3
# max_stage = 3
# c_region = data$c_region
# info_vector = data$info_vector
# z = data$z[3]
# theta_lower = -0.2
# theta_upper = 1
# tolerance = 1e-6

cmue <- function(stage, max_stage, c_region, info_vector, z, 
                 theta_lower = 0, theta_upper = 5, tolerance = .Machine$double.eps^0.25) {
  median_objective <- Vectorize(function(theta) {
    normalising_constant = psnorm(lower = -Inf, upper = Inf, stage = stage, max_stage = max_stage, c_region = c_region, 
                                  info_vector = info_vector, theta = theta)
    #print(normalising_constant)
    cdf_value <- 1/normalising_constant * psnorm(lower = z, upper = Inf, stage = stage, max_stage = max_stage, c_region = c_region, 
                        info_vector = info_vector, theta = theta)
    
    if(is.infinite(cdf_value) == TRUE){
      return(1e6)
    }
    if(is.na(cdf_value) == TRUE) {
      return(1e6)
    }
    return((cdf_value - 0.5)^2)
  },vectorize.args = "theta")

  tryCatch({
    root_result <- optimize(median_objective, interval = c(theta_lower, theta_upper), tol = tolerance)$minimum
    return(root_result)
    
  }, error = function(e) {
    return(paste("Root-finding failed:", e$message,
                 "\nConsider adjusting the bounds [theta_lower, theta_upper] or tolerance."))
  })
}

# cmue(stage = 1, max_stage = 2, c_region = matrix(c(-Inf, 2.787), ncol = 2, byrow = TRUE), info_vector = c(313, 394), z =2.9, tolerance = 1e-6, 
#   theta_lower = -1, theta_upper = 1)
# 
# plot(x, median_objective(x), ylim = c(0,0.61))
# 
# median_objective(0.4)

#data_test <- data_norm_tran_parallel(n_e = c(34,78), n_c = c(34,78), var_e = 25, var_c = 25, mu_e = c(5.2,6), mu_c = c(3,2.8), lower = -Inf, upper = 2.797)
#cmue(stage = 2, max_stage = 2, c_region = data_test$c_region, info_vector = data_test$info_vector, z = data_test$z[2], theta_lower = -1, theta_upper = 1, tolerance = 1e-5)
