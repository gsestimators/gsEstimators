create_correlation_matrix <- function(info_vector, enforce_monotone = TRUE) {
  i_len <- length(info_vector)
  
  if (enforce_monotone) {
    info_vector <- cummax(info_vector)
  }
  
  corr_mat <- matrix(NA_real_, ncol = i_len, nrow = i_len)
  
  for (i in 1:i_len) {
    for (j in 1:i_len) {
      if (i <= j) {
        corr_mat[i, j] <- sqrt(info_vector[i] / info_vector[j])
      } else {
        corr_mat[i, j] <- sqrt(info_vector[j] / info_vector[i])
      }
    }
  }
  
  return(corr_mat)
}