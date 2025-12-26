# Function to generate covariance matrix
cov_matrix <- function(times, stds, type, rho) {
  
  # Check input parameters
  n_times <- length(times)
  n_std <- length(stds)
  if (n_times != n_std) {
    stop(paste("ERROR: Number of measurement times (", n_times, ") does not equal
               number of standard deviation values (", n_std, ").", 
               sep = ""))
  }
  else if (toupper(type) != "AR" & toupper(type) != "CS") {
    stop(paste("ERROR: Covariance structure (", type, ") must be AR or CS.",
               sep = ""))
  }
  else if (rho < 0 | rho > 1) {
    stop("ERROR: rho must be between 0 and 1.")
  }
  
  # Initialize covariance matrix
  cov <- diag(stds * stds)
  n <- length(times)
  
  # Option 1: compound symmetry
  if (toupper(type) == "CS") {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {cov[i,j] <- stds[i] * stds[j] * rho}
      }
    }
  }
  
  # Option 2: first-order autoregressive 
  else if (toupper(type) == "AR") {
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {cov[i,j] <- stds[i] * stds[j] * rho ^ abs(i - j)}
      }
    }
  }
  
  # Output covariance matrix + times in last column
  return(cbind(cov, times)) 
  
}