#You can change the window here
window_start <- 2
window_end <- 6

alpha_matrix_window2 <- matrix(0, nrow = num_genes, ncol = num_genes)

for (gi in 1:num_genes) {
  X_matrix <- distance_matrix[window_start:(window_end - 2), ]
  Y_vector <- distance_matrix[(window_start + 1):(window_end - 1), gi]

  prev_alpha <- alpha_matrix_window1[, gi]  

  lambda1 <- 0.01  # check both lambdas. Best way to do cross-validation? 
  lambda2 <- 0.001  # smaller, better? 
  
  # Initial guess for alpha: can choose alpha from previous window
  #initial_alpha <- rep(0, ncol(X_matrix)) 
  
  objective_function <- function(alpha, X_matrix, Y_vector, lambda1, lambda2, prev_alpha) {
    residual <- Y_vector - X_matrix %*% alpha  
    sparsity_penalty <- lambda1 * sum(abs(alpha))  
    smoothness_penalty <- lambda2 * sum(abs(alpha - prev_alpha)) 
    return(sum(residual^2) + sparsity_penalty + smoothness_penalty)
  }
  
  optim_result <- optim(
    #par = initial_alpha,
    par=prev_alpha,
    fn = objective_function,
    X_matrix = X_matrix,
    Y_vector = Y_vector,
    lambda1 = lambda1,
    lambda2 = lambda2,
    prev_alpha = prev_alpha,
    method = "Nelder-Mead"
    #method = "L-BFGS-B"
    #method = "Nelder-Mead"
    #method="SANN"
    #method="CG" #good for sparse problems
    #method = "BFGS" #which method is a btter fit here?
  )
  alpha_matrix_window2[, gi] <- optim_result$par
}

print(alpha_matrix_window2)
alpha_matrix_window1
