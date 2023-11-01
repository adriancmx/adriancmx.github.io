heteroPCA <- function(data, rank, max_iterations, tolerance) {
  # Step 1: Initialize
  n <- nrow(data)
  p <- ncol(data)
  Sigma_hat <- cov(data)  # Sample covariance matrix
  N <- Sigma_hat
  diag(N) <- 0
  converged <- FALSE
  iteration <- 0
  
  while (!converged && iteration < max_iterations) {
    iteration <- iteration + 1
    
    # Step 2.1: Perform SVD on N
    svd_result <- svd(N)
    U <- svd_result$u
    V <- svd_result$v
    D <- svd_result$d
    
    # Rank-r approximation
    U_r <- U[, 1:rank]
    V_r <- V[, 1:rank]
    D_r <- diag(D[1:rank])
    delta_N <- U_r %*% D_r %*% t(V_r)
    
    # Step 2.2: Update N
    DelN <- N
    diag(DelN) <- 0
    N_new <- diag(diag(delta_N)) + DelN
    
    # Check for convergence
    if (max(abs(N - N_new)) < tolerance || iteration == max_iterations) {
      convergence <- TRUE
      break
    }
    
    N <- N_new
  }
  
  # Step 3: Output U
  U_hat <- U[, 1:rank]
  return(list(U = U_hat, convergence = convergence))
}

ordinaryPCA <- function(data, rank, max_iterations, tolerance) {
  n <- nrow(data)
  p <- ncol(data)
  cov_matrix <- cov(data)  # Step 1: Calculate the sample covariance matrix
  
  U <- NULL
  convergence <- FALSE
  for (t in 1:max_iterations) {
    # Step 2.1: Perform SVD on the covariance matrix
    svd_result <- svd(cov_matrix)
    lambda <- svd_result$d
    U <- svd_result$u
    
    # Get the best rank-r approximation
    Ur <- U[, 1:rank]
    lambdar <- diag(lambda[1:rank])
    delta_cov <- Ur %*% lambdar %*% t(Ur)
    
    # Step 2.2: Update the covariance matrix
    cov_matrix_new <- delta_cov
    
    # Step 2.3: Increment iteration count
    t <- t + 1
    
    # Check for convergence
    if (max(abs(cov_matrix - cov_matrix_new)) < tolerance || t == max_iterations) {
      convergence <- TRUE
      break
    }
    
    cov_matrix <- cov_matrix_new
  }
  
  # Step 3: Output the matrix U and convergence status
  U <- U[, 1:rank]
  U <- t(U)  # Transpose to have dimensions n x r
  return(list(U = U, convergence = convergence))
}

# Sample data for testing
set.seed(123)
Y <- matrix(rnorm(100), nrow = 10)  # 10 observations, 10 variables

# Apply heteroPCA on the sample data
rank <- 3  # Desired rank
max_iter <- 100  # Maximum number of iterations
tolerance <- 1e-6  # Replace with the tolerance for convergence

U <- heteroPCA(Y, rank, max_iter, tolerance)

# Display the estimated principal components
print(U)

# Perform ordinary PCA
result <- ordinaryPCA(data, rank, max_iter, tolerance)

# Access the results
print(result$U)
print(result$convergence)