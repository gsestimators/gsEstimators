#' Conditional Maximum Likelihood Estimator
#'
#' Computes the conditional maximum likelihood estimator (CMLE) of the effect size 
#' \eqn{\theta} given that the trial has continued to stage \eqn{k}, and the final 
#' test statistic \eqn{Z_k = z} was observed.
#'
#' @param k Integer, final stage of the trial where it stopped
#' @param c_region Matrix with 2 columns and \eqn{k - 1} rows, specifying the continuation 
#'        region bounds for stages 1 to \eqn{k-1}. Column 1: lower bounds, Column 2: upper bounds
#' @param info_vector Numeric vector of information levels at each stage (length \eqn{k})
#' @param z Numeric, observed value of the test statistic \eqn{Z_k}
#' @param theta_lower Numeric, lower bound for \eqn{\theta} in optimization (default: 0)
#' @param theta_upper Numeric, upper bound for \eqn{\theta} in optimization (default: 5)
#' @param theta_initial Numeric, initial value for \eqn{\theta} to start optimization
#' @param eps Numeric, small value used for numerical differentiation by finite differences in `dsnorm()` (default: 1e-5)
#'
#' @return Numeric value representing the conditional MLE of \eqn{\theta}
#'
#' @details
#' This function finds the value of \eqn{\theta} that maximizes the conditional 
#' likelihood of the observed final-stage test statistic \eqn{Z_k = z}, given that the 
#' trial reached stage \eqn{k}. It does so by maximizing the conditional density computed 
#' by the `dsnorm()` function.
#'
#' For \eqn{k = 1}, this reduces to standard MLE under a single-stage normal model.
#' For \eqn{k > 1}, the function numerically maximizes the conditional likelihood 
#' using the Brent optimization method.
#'
#' @references
#' Jennison, C. and Turnbull, B. W. (2000). 
#' \emph{Group Sequential Methods with Applications to Clinical Trials}. 
#' Chapman & Hall/CRC. Chapter 2 and 3.
#'
#' @examples
#' # Example for stage 1 (no continuation region)
#' cmle(k = 1, c_region = NULL, info_vector = c(0.5, 1.0), z = 1.5, 
#'      theta_initial = 0.5)
#'
#' # Example for stage 2 with continuation region
#' c_region <- matrix(c(-1, 1), nrow = 1, ncol = 2)
#' cmle(k = 2, c_region = c_region, info_vector = c(0.5, 1.0), z = 1.2,
#'      theta_lower = -2, theta_upper = 2, theta_initial = 0)
#'
cmle <- function(k, c_region, info_vector, z, 
                 theta_lower = 0, theta_upper = 5, theta_initial, eps = 1e-5) {
  
  theta_objective <- function(theta) {
    return(dsnorm(x = z, k = k, c_region = c_region, info_vector = info_vector, theta = theta, eps = eps))
  }

  tryCatch({

    max_result <- optim(theta_initial, theta_objective, method = "Brent", lower = theta_lower, upper = theta_upper, control = list(fnscale = -1))
    return(max_result$par)
    
  }, error = function(e) {
    stop("Optimisation failed: ", e$message, 
         "\nConsider adjusting the bounds [theta_lower, theta_upper], theta_initial or tolerance.")
  })
}
