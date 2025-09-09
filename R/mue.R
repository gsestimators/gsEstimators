#' Calculate MUE (Median Unbiased Estimator) for Group Sequential Design
#'
#' Computes the median unbiased estimator by finding the parameter value θ such
#' that the overall rejection probability equals 0.5 under the alternative
#' hypothesis. This provides an unbiased estimate of the treatment effect that
#' adjusts for the bias introduced by sequential monitoring.
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
#' @return Numeric, the median unbiased estimate of θ
#'
#' @details
#' The median unbiased estimator solves for θ such that P_θ(reject H_0) = 0.5,
#' where the rejection probability is computed across all k stages:
#'
#' \itemize{
#'   \item Stage 1: P(Z_1 ∉ C_1) = P(Z_1 < c_{1,1}) + P(Z_1 > c_{1,2})
#'   \item Stages 2 to k-1: P(continue to stage i, then reject at stage i)
#'   \item Stage k: P(continue to stage k, then |Z_k| > z)
#' }
#'
#' The median unbiased estimator has the property that P(θ̂ > θ) = P(θ̂ < θ) = 0.5
#' for the true parameter value θ, providing an unbiased estimate that accounts
#' for the sequential nature of the design.
#'
#' @note
#' This function assumes the final stage uses a two-sided test with critical
#' value z. The continuation regions should be properly defined for stages 1 to k.
#'
#' @references
#' Jennison, C. and Turnbull, B. W. (2000).
#' \emph{Group Sequential Methods with Applications to Clinical Trials}.
#' Chapman & Hall/CRC. Chapter 7.
#'
#' @examples
#' # Example with 3-stage design
#' k <- 3
#' c_region <- matrix(c(-1, 1, -1.5, 1.5, -2, 2), nrow = (k - 1), ncol = 2)
#' info_vector <- c(0.33, 0.67, 1.0)
#' z <- 1.96
#'
#' mue_estimate <- mue(k = k, c_region = c_region,
#'                     info_vector = info_vector, z = z)
#'@export

mue <- function(stage, max_stage, c_region, info_vector, z, tolerance = 1e-6,
                theta_lower = 0, theta_upper = 5, tail_type = "upper", max_iter = 1000, n = 100) {


  if(theta_lower == theta_upper) {
    return("Lower and upper search regions must be different see control inputs")
  }


  # Inner function: compute rejection probability for given theta
  p_reject <- Vectorize( function(theta) {

    # Initialize probability parts for each stage
    p_parts <- rep(0, times = stage)
    p_lower_tail <- rep(0, times = stage)
    p_upper_tail <- rep(0, times = stage)

    # Stage 1: Direct rejection probability
    for (i in 1:stage) {
      if (i == 1) {
        if(tail_type == "upper") {
          up <- c_region[1, 2]

          p_upper <- psnorm(lower = up, upper = Inf, stage = 1,
                            max_stage = 1,
                            info_vector = info_vector[1], theta = theta)

          p_parts[i] <- p_upper

          next
        }
        if(tail_type == "lower") {
          lo <- c_region[1, 1]

          p_lower <- psnorm(lower = -Inf, upper = lo, stage = 1,
                            max_stage = 1,
                            info_vector = info_vector[1], theta = theta)

          p_parts[i] <- p_lower

          next
        }
        if(tail_type == "two") {

          up <- c_region[1, 2]
          lo <- c_region[1, 1]

          p_upper <- psnorm(lower = up, upper = Inf, stage = 1,
                            max_stage = 1,
                            info_vector = info_vector[1], theta = theta)

          p_lower <- psnorm(lower = -Inf, upper = lo, stage = 1,
                            max_stage = 1,
                            info_vector = info_vector[1], theta = theta)

          p_upper_tail[i] <- p_upper
          p_lower_tail[i] <- p_lower

          next
        }
      }

      # Final stage k: Two-sided test with observed value z
      if (i == stage) {

        if(stage == max_stage) {
          if(tail_type == "lower") {

            p_lower <- psnorm(lower = -Inf, upper = z, stage = i, max_stage = max_stage,
                         c_region = c_region[ , , drop = FALSE],
                         info_vector = info_vector[1:i], theta = theta)

            p_parts[i] <- p_lower

            next
          }
          if(tail_type == "upper") {

            p_upper <- psnorm(lower = z, upper = Inf, stage = i, max_stage = max_stage,
                         c_region = c_region[ , , drop = FALSE],
                         info_vector = info_vector[1:i], theta = theta)

            p_parts[i] <- p_upper

            next
          }
          if(tail_type == "two") {

            p_upper <- psnorm(lower = z, upper = Inf, stage = i, max_stage = max_stage,
                         c_region = c_region[ , , drop = FALSE],
                         info_vector = info_vector[1:i], theta = theta)

            p_lower <- psnorm(lower = -Inf, upper = z, stage = i, max_stage = max_stage,
                         c_region = c_region[ , , drop = FALSE],
                         info_vector = info_vector[1:i], theta = theta)

            p_upper_tail[i] <- p_upper
            p_lower_tail[i] <- p_lower

            next

          }
        }

        if(tail_type == "lower") {

          p_lower <- psnorm(lower = -Inf, upper = z, stage = i, max_stage = max_stage,
                            c_region = c_region[1:i , , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_parts[i] <- p_lower

          next
        }
        if(tail_type == "upper") {

          p_upper <- psnorm(lower = z, upper = Inf, stage = i, max_stage = max_stage,
                            c_region = c_region[ 1:i, , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_parts[i] <- p_upper

          next
        }
        if(tail_type == "two") {

          p_upper <- psnorm(lower = z, upper = Inf, stage = i, max_stage = max_stage,
                            c_region = c_region[ 1:i, , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_lower <- psnorm(lower = -Inf, upper = z, stage = i, max_stage = max_stage,
                            c_region = c_region[1:i , , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_upper_tail[i] <- p_upper
          p_lower_tail[i] <- p_lower

          next

        }
      }

      # Intermediate stages 2 to k-1: Rejection outside continuation region
      else {
        if(tail_type == "upper") {

          up <- c_region[i, 2]

          p_upper <- psnorm(lower = up, upper = Inf, stage = i,
                            max_stage = max_stage,c_region = c_region[1:i, , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_parts[i] <- p_upper

          next
        }
        if(tail_type == "lower") {

          lo <- c_region[i, 1]

          p_lower <- psnorm(lower = -Inf, upper = lo, stage = i,
                            max_stage = max_stage,c_region = c_region[1:i, , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_parts[i] <- p_lower

          next
        }
        if(tail_type == "two") {

          up <- c_region[i, 2]
          lo <- c_region[i, 1]

          p_upper <- psnorm(lower = up, upper = Inf, stage = i,
                            max_stage = max_stage,c_region = c_region[1:i, , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_lower <- psnorm(lower = -Inf, upper = lo, stage = i,
                            max_stage = max_stage,c_region = c_region[1:i, , drop = FALSE],
                            info_vector = info_vector[1:i], theta = theta)

          p_upper_tail[i] <- p_upper
          p_lower_tail[i] <- p_lower

          next
        }
      }
    }

    if(tail_type == "two"){
      return((2 * min(sum(p_lower_tail), sum(p_upper_tail))) - 0.5)
    }

    else{
    return(sum(p_parts) - 0.5)
    }

  }, vectorize.args = "theta")

  mue_result <- uniroot.all(p_reject,
                            lower = theta_lower,
                            upper = theta_upper,
                            tol = tolerance,
                            maxiter = max_iter,
                            n = n)

  return(mue_result)
}
